/* 
 * Copyright (c) 2026 Nick Xenias
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU Lesser General Public License as   
 * published by the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * Lesser General Lesser Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
//==============================================================================
// Name:     evm.h
//
// Purpose:  Implement efficient vector math functions using SIMD instructions
//
// Created:  2026/03/28
//
// Description:
//  Implements some common efficient vector math functions for x86 processors
//  using SIMD (single instruction multiple data) intrinsics.  These functions
//  use runtime dispatch to automatically choose the most advanced instruction
//  set based on the architecture being run on.  For example, if a function
//  supports both an AVX and an AVX512 branch, it will take the AVX512 branch
//  if executed on an Icelake cpu or the AVX branch if executed on a Haswell
//  cpu, regardless of compiler optimization flags.  Note that use of this
//  header requires you to compile with c++17.  In some cases (unit/regression
//  testing, etc.) it may be desireable to disable certain branches.  This can
//  be done by specifying one or more of the following macros at compile time:
//    -DDISABLE_AVX512
//    -DDISABLE_AVX2
//    -DDISABLE_AVX
//    -DDISABLE_SSE4 (disables both 4.1 and 4.2)
//    -DDISABLE_SSSE3
//    -DDISABLE_SSE3
//    -DDISABLE_SSE2
//    -DDISABLE_SSE
//    -DDISABLE_SIMD (disables all branches, uses single-element "for" loops)
//
//  Note that currently there are only implementations for AVX, AVX2, and AVX512
//  branches along with a default (non-SIMD) branch, so disabling the SSE
//  branches will not have any practical effect.
//
//  One of the hidden costs of SIMD instructions is loading data from memory
//  into a SIMD register and storing it back into memory after a computation.
//  Loads and stores are more efficient if memory is aligned on address
//  boundaries that are a multiple of the SIMD register size.  To simplify this
//  a memory-aligned allocator is provided which is compatible with c++ stl
//  containers (vector, etc.) using aligned_allocator.h.  Furthermore, if the
//  vector itself is a multiple of the SIMD register size in bytes, many
//  functions can skip costly edge-condition checking to ensure memory is not
//  over-indexed on the final load/store operation of a vector.  For example,
//  if you have a <float> type vector with 7 elements, using SSE
//  instructions will need to do 2 loads, the first loading elements 0-3 and
//  the second loading 4-7.  However since the vector is only 7 elements long
//  there is no element 7 and it will be over-indexed by the second load.  This
//  requires an edge case check at the end of each buffer which typically
//  involves a mem copy to a location that has sufficient size for the load or
//  store.  Newer x86 instructions sets have masked load/store operations which
//  avoid this, however these are not available on older architectures.
//
//  To help mitigate the cost of load operations, this code base includes an
//  aligned allocator for stl vectors (aligned_allocator.h) which both aligns
//  memory on boundaries that are a multiple of the SIMD register size and
//  reserves additional space in the buffer equivalent to one full SIMD
//  register.  There is also a helper function to return the number of elements
//  per SIMD register based on data type (float, complex<float>, etc.).  Below
//  is an example of creating a SIMD aligned vector of complex single precision
//  floating point elements, sized to a multiple of the SIMD register size:
//    ALIGNED_VECTOR(complex<float> ) myvec(EVM::getSize<complex<float> >(29));
//  Note that the EVM::getSize<complex<float> >(29) function call will return
//  the next higher number (>=29) of 64-bit (complex float) elements that is a
//  multiple of the SIMD register size on the target platform. 
//
//  For maximum speed/efficiency, the function EVM::enableFastDenormals()
//  should be called at the start of any program using evm.h.  This will avoid
//  large penalties when denormal (very small floating point) numbers are
//  produced by setting them to 0.  Setting this at the start of the program
//  will cause any threads spawned by the program to inherit the setting as
//  well.
//  
//  The following mathematical operations are currently implemented:
//
//  Multiply/Divide:
//    scale = multiply an input buffer by a constant scale factor
//    mult = point-by-point vector multiply
//    multc = point-by-point complex conjugate multiply, second input conjugated
//    div = point-by-point vector divide, true division (this is the preferred
//        implementation for complex divisors, achieving the most accuracy with
//        speeds as fast as divr on most AVX-enabled CPUs)
//    divnr = point-by-point vector divide using a reciprocal approximation with
//        1 Newton-Raphson stage for refinement, this is the preferred
//        implementation for moderate buffer sizes (fit in cache) with a
//        reasonable tradeoff between accuracy (within 2 ulp of true division
//        with AVX512, 3 ulp with AVX2/AVX, and within 1 ulp ~97-99% of the
//        time) and speed
//    divr = point-by-point vector divide using fast reciprocal approximation
//        without Newton-Raphson refinement, this favors speed more than
//        accuracy and can result in maximum error of ~512 ulp near 1.0 for
//        AVX512 (6.1e-5 normalized absolute error) or 3072 ulp near 1.0 for
//        AVX2/AVX (3.7e-4 normalized absolute error)
//    divc1 = point-by-point complex conjugate divide, first input conjugated,
//        true division
//    divc1nr = point-by-point complex conjugate divide, first input
//        conjugated, with 1 Newton-Raphson stage for refinement, this is the
//        preferred implementation for moderate buffer sizes (fit in cache)
//        with real divisors for a reasonable tradeoff between accuracy and
//        speed (for complex divisors div is more accurate and as fast or
//        faster on certain CPUs)
//    divc1r = point-by-point complex conjugate divide, first input conjugated,
//        using fast reciprocal approximation without Newton-Raphson refinement
//    divc2 = point-by-point complex conjugate divide, second input conjugated,
//        true division
//    divc2nr = point-by-point complex conjugate divide, second input
//        conjugated, with 1 Newton-Raphson stage for refinement, this is the
//        preferred implementation for moderate buffer sizes (fit in cache)
//        with a reasonable tradeoff between accuracy and speed
//    divc2r = point-by-point complex conjugate divide, second input conjugated,
//        using fast reciprocal approximation without Newton-Raphson refinement
//
//  Dot-product:
//    dotp = dot-product of two input buffers
//    dotpc = complex-complex conjugate dot-product, second input conjugated
//
//  Power/Magnitude:
//    pow2 = squares each element of an input buffer
//    pow4 = raises each element of an input buffer to the 4th power
//    pow8 = raises each element of an input buffer to the 8th power
//    mag = magnitude of input buffer
//    magr = magnitude of a complex input buffer using a reciprocal square root
//        root approximation (trades accuracy for speed)
//    magSq = magnitude squared of a complex input buffer, faster than mag as
//        it doesn't have to compute a square root
//
//  Specialized:
//    tune = tunes complex data with a numerically controlled oscillator (NCO)
//    fs2Shift = Fs/2 frequency shift for complex data
//
//  Important notes on speed/accuracy:
//    1. Some division methods compute 2 registers worth of divides at once and
//       others do not (this may even differ between the AVX512 and AVX/AVX2
//       implementations).  This is done because the divisor for 2 dividend
//       registers can be packed into a single register (i.e. with a complex
//       divisor) which can reduce the number of div or rcp calls.  For AVX512
//       this is typically faster but for AVX/AVX2 the extra shuffle and
//       multiply instructions can actually use more cycles than the rcp call.
//    2. On some CPUs the div function may be just as fast or slightly faster
//       than the reciprocal approximation versions (divr, divnr).  This is
//       because the division operation itself can be pipelined while the other
//       instructions are running, whereas the reciprocal approximation
//       versions are using shuffles and multiplies which use the same ports
//       as the other instructions (which are also shuffles and multiplies).
//    3. The "true" division methods involving complex dividend, divisor, or
//       both actually compute 1/divisor then multiply the dividend by this.
//       This greatly improves the speed (see note 1) but comes at the cost of
//       slightly lower accuracy (up to 1.5 ULP vs. .5 ULP for dividend/divisor)
//       and range restrictions (1/divisor can overflow or underflow if the
//       magnitude of the divisor is < 2.9e-39 or > 8.5e37).
//  
//  Important Compiler-dependent notes (most users will not care about this,
//  but here it is anyway):
//    1. FP exceptions: for speed, some methods (i.e. division) load the
//       final partial register with zeros in the unused elements and divide
//       them anyway (0/0, or 0*rcp(0)), which raises FE_INVALID (SIGFPE if FP
//       exceptions are unmasked) and/or FE_DIVBYZERO even on valid input.
//       Those elements are never stored, so results are unaffected; this only
//       matters to code that tests or traps FP exceptions.  Note that even
//       when masked instructions are used, icpx's default floating-point model
//       computes them anyway, so they can raise FE_INVALID as well (icpx's
//       strict model, which would prevent this, cannot be combined with the
//       target-attribute dispatch used here).  Applying setUnusedToOne() to
//       the divisor avoids it at the cost of 2 extra instructions per call
//       with a partial register.
//    2. Complex conj(),*,/ operators: icpx 2026.1 does not emit the
//       std::complex<float>(re, im) constructor when it is only used from
//       target()-versioned functions, which leaves an undefined reference at
//       link time.  conj() and arithmetic between a float and a complex both
//       use that constructor, so the non-SIMD ("default") versions that would
//       need it either expand the computation out with .real()/.imag() or
//       copy a constexpr complex (built at compile time) and set its real()
//       and imag() parts from the inputs.  The division methods do not use
//       expanded computations as this would result in the same range
//       restrictions as the SIMD paths (i.e. very large or very small inputs
//       yielding NaN), and the intent of the default methods is to give exact
//       performance to the native / operator at the expense of speed.
//
//==============================================================================
#ifndef EFFICIENT_VECTOR_MATH_H
#define EFFICIENT_VECTOR_MATH_H

#include <math.h>
#include <cstdint>   // for int32_t
#include <string>
#include <cstring>   // for memcpy/memset
#include <algorithm> // for std::min
#include <complex>
#include <type_traits> // for static_assert type checking
#include <x86intrin.h>
#include <vector>
#include "aligned_allocator.h"

// check if user wishes to disable all SIMD and explicittly set the macros
#if defined(DISABLE_SIMD)
#define DISABLE_AVX512
#define DISABLE_AVX2
#define DISABLE_AVX
#define DISABLE_SSE4
#define DISABLE_SSSE3
#define DISABLE_SSE3
#define DISABLE_SSE2
#define DISABLE_SSE
#endif // DISABLE_SIMD

// for memory allocation and register size, start with largest and work down
// to smallest, as long as an instruction set is not explicitly disabled,
// assume it can/will be used (i.e. account for worst-case alignment and size)
#if !defined(DISABLE_AVX512)
  typedef evm_64 evm_align;
  #define MAX_SIMD_REG 64 // zmm 64-byte
#elif !defined(DISABLE_AVX2) || !defined(DISABLE_AVX)
  typedef evm_32 evm_align;
  #define MAX_SIMD_REG 32  // ymm 32-byte
#elif !defined(DISABLE_SSE4) || !defined(DISABLE_SSSE3) || !defined(DISABLE_SSE3) || !defined(DISABLE_SSE2) || !defined(DISABLE_SSE)
  typedef evm_16 evm_align;
  #define MAX_SIMD_REG 16  // xmm 16-byte
#else
  // default case, no extra bytes but align start by 16 to maybe get some extra
  // speedup from automatic vectorization by compiler
  typedef evm_16 evm_align;
  #define MAX_SIMD_REG 0  // no SIMD registers used
#endif

// define some aliases for using memory aligned allocators based on the SIMD
// registers being used
#ifndef ALIGNED_VECTOR
#define ALIGNED_VECTOR(type) std::vector<type, aligned_allocator<type, evm_align> >
#endif
#ifndef ALIGNED_ITERATOR
#define ALIGNED_ITERATOR(type) std::vector<type, aligned_allocator<type, evm_align> >::iterator
#endif

namespace EVM
{
  using std::complex;
  
  // mask table for the final partial register of the AVX/AVX2 kernels, rows
  // can be cast as __m256i* directly and row n selects the first n floats
  alignas(32) inline constexpr int32_t masks[9][8] = {
    { 0, 0, 0, 0, 0, 0, 0, 0},   // remainder 0 (unused)
    {-1, 0, 0, 0, 0, 0, 0, 0},   // remainder 1
    {-1,-1, 0, 0, 0, 0, 0, 0},   // remainder 2
    {-1,-1,-1, 0, 0, 0, 0, 0},   // remainder 3
    {-1,-1,-1,-1, 0, 0, 0, 0},   // remainder 4
    {-1,-1,-1,-1,-1, 0, 0, 0},   // remainder 5
    {-1,-1,-1,-1,-1,-1, 0, 0},   // remainder 6
    {-1,-1,-1,-1,-1,-1,-1, 0},   // remainder 7
    {-1,-1,-1,-1,-1,-1,-1,-1}    // remainder 4 complex elements (masks[rem<<1])
  };

  // Returns the mask for AVX512 masked loads/stores of the final partial
  // register: the low n bits are set (0 <= n <= 16), so floats 0..n-1 are
  // loaded/stored and the rest are skipped (maskz loads set them to 0)
  __attribute__((__always_inline__))
  constexpr __mmask16 mask16(const int n)
  {
    return static_cast<__mmask16>(0xffff >> (16 - n));
  }

  //==================//
  // Helper Functions //
  //==================//
  // Returns a buffer size in elements that is greater than or equal to the
  // desired size but would result in the buffer being a multiple of
  // EVM_ALIGNMENT bytes.  This can be used to compute required buffer sizes
  // for using the more efficient internal methods which require the length
  // argument to be a multiple of the SSE/AVX register size.  Be careful to
  // check the specific methods to determine whether the length pertains to
  // Complex8 or float data types, as some methods have mixed inputs.
  //    desired = desired buffer size in elements
  //    EVM::getSize = buffer size closest to desired which is greater than
  //        or equal to desired and is a multiple of the AVX/SSE register size
  template<typename T>
  inline size_t getSize(const size_t desired);
  
  // Returns the highest optimization level as a string, note that in order to
  // output a string, this must be a wrapper around a const char* method as
  // the intel compiler (based on LLVM) doesn't instantiate the std::string
  // constructor when multiversioning is used
  //    EVM::getOptLevel = string representation of highest optimization level
  inline const char* getOptLevelRaw();
  inline std::string getOptLevel();

  // Enables fast denormals (floating point numbers very close to 0) by
  // flushing produced denormals to 0 and interpreting them as 0.  This can
  // result in large savings (on the order of 80 cycles) but must be set before
  // any threads are spawned in order for threads to inherit the setting.  This
  // method should be called at the start of any program using evm.h unless
  // expressly undesired by other functionality in the program.  This same
  // functionality will usually be enabled automatically with -ffast-math but
  // best not to rely on that.
  inline void enableFastDenormals()
  {
    #if defined(__SSE__) || defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    #endif
  };

  // Changes the lanes of an AVX register that msk does not select from 0.0f to
  // 1.0f to avoid things like divide by 0 tripping FE_INVALID on unused
  // elements (_mm256_maskload_ps sets unused elements to 0, this is intended
  // to set the unused elements on the output of maskload to 1).  While
  // division by 0 in the unused elements is not an error per-se, it would
  // raise FE_INVALID (SIGFPE if FP exceptions are unmasked) even though values
  // are never stored or used for output.  If a user or debugger is trapping
  // that error it could potentially mask a true FE_INVALID error later in the
  // program by triggering on these unused elements.  Note that this
  // implementation deliberately does not use _mm256_blendv_ps as certain
  // compilers (GCC 12) can expand that in AVX code to a branch per element
  // resulting in significantly slower code.  While the current implementation
  // avoids this branching, it only works for 0-valued elements.  Currently
  // unused: the division methods skip it for speed (see note 1 at the top of
  // this file), it is kept for code where FE_INVALID-clean tails matter.
  //    in = loaded data to mask, replacing 0's with 1's
  //    msk = mask to select which values are replaced with 1's, see
  //         _mm256_maskload_ps for details
  //    setUnusedToOne = input register with any 0-valued elements not selected
  //        by msk set to 1.0f
  __attribute__((__target__("avx"), __always_inline__))
  inline __m256 setUnusedToOne(const __m256 in, const __m256i msk)
  {
    return _mm256_or_ps(in, _mm256_andnot_ps(_mm256_castsi256_ps(msk),
        _mm256_set1_ps(1.0f)));
  }

  //==============//
  // Vector Scale //
  //==============//
  // Multiply each element of a vector by a constant scale factor.  The output
  // buffer may be the same as the input buffer.
  //    in1 = buffer to apply scale factor to
  //    scaleFactor = scale factor to apply to in1
  //    len = number of elements to scale
  //    out = output buffer, product of in1*scaleFactor
  inline void scale(float const * const in1, const float scaleFactor,
      const int len, float * const out);
  inline void scale(complex<float> const * const in1, const float scaleFactor,
      const int len, complex<float> * const out);

  //=================//
  // Vector Multiply //
  //=================//
  // Point by point multiply of two input buffers.  The output buffer may be
  // the same as one of the input buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply
  //    len = number of elements to multiply
  //    out = output buffer to recieve result of in1*in2
  inline void mult(float const * const in1, float const * const in2,
      const int len, float * const out);
  inline void mult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  inline void mult(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  
  // Point by point complex-conjugate multiply of two input buffers, second
  // input conjugated.  The output buffer may be the same as one of the input
  // buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply, conjugated
  //    len = number of elements to multiply
  //    out = output buffer to recieve result of complex multiply
  inline void multc(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  inline void multc(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  
  //============================//
  // Vectory Multiply and Scale //
  //============================//
  // Point by point multiply of two input buffers with a constant scale factor.
  // The output buffer may be the same as one of the input buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply
  //    len = number of elements to multiply
  //    scale = scale factor
  //    out = output buffer to recieve result of complex multiply
  inline void mults(float const * const in1, float const * const in2,
      const int len, const float scale, float * const out);
  inline void mults(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out);
  inline void mults(float const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out);
  // Point by point complex conjugate multiply of two input buffers with a
  // constant scale factor.  The output buffer may be the same as one of the
  // input buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply, conjugated
  //    len = number of elements to multiply
  //    scale = scale factor
  //    out = output buffer to recieve result of complex multiply
  inline void multcs(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out);

  //===============//
  // Vector Divide //
  //===============//
  // Point by point divide of two input buffers.  The output buffer may be
  // the same as one of the input buffers.  This version uses a true divide.
  // Note that for complex divisors, AVX results are only reliable while |in2|
  // is roughly between 1e-19 and 1.8e19 and |in1|·|in2| stays below 3.4e38.
  // Outside that range the quotient can be 0, INF or NaN, especially when
  // enableFastDenormals() has been called.
  //    in1 = first input buffer for divide (dividend)
  //    in2 = second input buffer for divide (divisor)
  //    len = number of elements to divide
  //    out = output buffer to recieve result of complex division
  inline void div(float const * const in1, float const * const in2,
      const int len, float * const out);
  inline void div(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  inline void div(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  inline void div(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out);
  
  // Point by point divide of two input buffers.  The output buffer may be
  // the same as one of the input buffers.  This version uses a reciprocal
  // approximation with 1 Newton-Raphson stage for refinement to achieve an
  // answer comparable to true division for single precision float elements but
  // with increased efficiency.  Note that due to the reciprocal approximation
  // and Newton-Raphson stage, these methods may produce different results than
  // normal division when the divisor is a special floating point value (0
  // yields NaN rather than INF, INF yields NaN rather than 0, and subnormals
  // may yield -+INF with the wrong sign).  Also note that for complex
  // divisors, AVX results are only reliable while |in2| is roughly between
  // 1e-19 and 9e18 and |in1|·|in2| stays below 3.4e38.  Outside that range
  // the quotient can be 0, INF or NaN, especially when enableFastDenormals()
  // has been called.
  //    in1 = first input buffer for divide (dividend)
  //    in2 = second input buffer for divide (divisor)
  //    len = number of elements to divide
  //    out = output buffer to recieve result of complex division
  inline void divnr(float const * const in1, float const * const in2,
      const int len, float * const out);
  inline void divnr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  inline void divnr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  inline void divnr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out);

  // Point by point divide of two input buffers using a reciprocal
  // approximation without Newton-Raphson refinement.  The output buffer may be
  // the same as one of the input buffers.  This version is more efficient than
  // the normal div() and divnr() functions but can have larger error from a
  // true divide.   Note that for complex divisors, AVX results are only
  // reliable while |in2| is roughly between 1e-19 and 9e18 and |in1|·|in2|
  // stays below 3.4e38.  Outside that range the quotient can be 0, INF or NaN,
  // especially when enableFastDenormals() has been called.
  //    in1 = first input buffer for divide (dividend)
  //    in2 = second input buffer for divide (divisor)
  //    len = number of elements to divide
  //    out = output buffer to recieve result of complex division
  inline void divr(float const * const in1, float const * const in2,
      const int len, float * const out);
  inline void divr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  inline void divr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  inline void divr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out);
  
  // Point by point complex-conjugate divide of two input buffers, first
  // input conjugated.  The output buffer may be the same as one of the input
  // buffers.  This version uses a true divide.  Note that for complex
  // divisors, AVX results are only reliable while |in2| is roughly between
  // 1e-19 and 1.8e19 and |in1|·|in2| stays below 3.4e38.  Outside that range
  // the quotient can be 0, INF or NaN, especially when enableFastDenormals()
  // has been called.
  //    in1 = first input buffer for divide (dividend), conjugated
  //    in2 = second input buffer for divide (divisor)
  //    len = number of elements to divide
  //    out = output buffer to receive result of complex division
  inline void divc1(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  inline void divc1(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out);
  // Point by point complex-conjugate divide of two input buffers, second
  // input conjugated.  The output buffer may be the same as one of the input
  // buffers. This version uses a true divide.  Note that for complex divisors,
  // AVX results are only reliable while |in2| is roughly between 1e-19 and
  // 1.8e19 and |in1|·|in2| stays below 3.4e38.  Outside that range the
  // quotient can be 0, INF or NaN, especially when enableFastDenormals() has
  // been called.
  //    in1 = first input buffer for divide (dividend)
  //    in2 = second input buffer for divide (divisor), conjugated
  //    len = number of elements to divide
  //    out = output buffer to receive result of complex division
  inline void divc2(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  inline void divc2(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  
  // Point by point complex-conjugate divide of two input buffers, first
  // input conjugated. This version uses a reciprocal approximation with 1
  // Newton-Raphson stage for refinement to achieve an answer comparable to
  // true division for single precision float elements but with increased
  // efficiency.  Note that due to the reciprocal approximation and
  // Newton-Raphson stage, these methods may produce different results than
  // normal division when the divisor is a special floating point value (0
  // yields NaN rather than INF, INF yields NaN rather than 0, and subnormals
  // may yield -+INF with the wrong sign).  Also Note that for complex
  // divisors, AVX results are only reliable while |in2| is roughly between
  // 1e-19 and 9e18 and |in1|·|in2| stays below 3.4e38.  Outside that range
  // the quotient can be 0, INF or NaN, especially when enableFastDenormals()
  // has been called.
  //    in1 = first input buffer for divide (dividend), conjugated
  //    in2 = second input buffer for divide (divisor)
  //    len = number of elements to divide
  //    out = output buffer to receive result of complex division
  inline void divc1nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  inline void divc1nr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out);
  // Point by point complex-conjugate divide of two input buffers, second
  // input conjugated.  This version uses a reciprocal approximation with 1
  // Newton-Raphson stage for refinement to achieve an answer comparable to
  // true division for single precision float elements but with increased
  // efficiency.  Note that due to the reciprocal approximation and
  // Newton-Raphson stage, these methods may produce different results than
  // normal division when the divisor is a special floating point value (0
  // yields NaN rather than INF, INF yields NaN rather than 0, and subnormals
  // may yield -+INF with the wrong sign).  Also note that for complex
  // divisors, AVX results are only reliable while |in2| is roughly between
  // 1e-19 and 9e18 and |in1|·|in2| stays below 3.4e38.  Outside that range
  // the quotient can be 0, INF or NaN, especially when enableFastDenormals()
  // has been called.
  //    in1 = first input buffer for divide (dividend)
  //    in2 = second input buffer for divide (divisor), conjugated
  //    len = number of elements to divide
  //    out = output buffer to receive result of complex division
  inline void divc2nr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  inline void divc2nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);

  // Point by point complex-conjugate divide of two input buffers, first
  // input conjugated, using a reciprocal approximation without Newton-Raphson
  // refinement.  Note that for complex divisors, AVX results are only reliable
  // while |in2| is roughly between 1e-19 and 9e18 and |in1|·|in2| stays below
  // 3.4e38.  Outside that range the quotient can be 0, INF or NaN, especially
  // when enableFastDenormals() has been called.
  //    in1 = first input buffer for divide (dividend), conjugated
  //    in2 = second input buffer for divide (divisor)
  //    len = number of elements to divide
  //    out = output buffer to receive result of complex division
  inline void divc1r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  inline void divc1r(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out);
  // Point by point complex-conjugate divide of two input buffers, second
  // input conjugated, using a reciprocal approximation without Newton-Raphson
  // refinement.  Note that for complex divisors, AVX results are only reliable
  // while |in2| is roughly between 1e-19 and 9e18 and |in1|·|in2| stays below
  // 3.4e38.  Outside that range the quotient can be 0, INF or NaN, especially
  // when enableFastDenormals() has been called.
  //    in1 = first input buffer for divide (dividend)
  //    in2 = second input buffer for divide (divisor), conjugated
  //    len = number of elements to divide
  //    out = output buffer to receive result of complex division
  inline void divc2r(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  inline void divc2r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  
  //=============//
  // Dot-Product //
  //=============//
  // Dot-product of two input buffers
  //    in1 = first input buffer for dot-product
  //    in2 = second input buffer for dot-product
  //    len = number of elements to compute dot-product over
  //    EVM::dotp = dot-product of in1 and in2
  inline float dotp(float const * const in1, float const * const in2,
      const int len);
  inline complex<float> dotp(complex<float> const * const in1,
      complex<float> const * const in2, const int len);
  inline complex<float> dotp(float const * const in1,
      complex<float> const * const in2, const int len);
  
  // Point by point complex-conjugate dot-product of two input buffers, second
  // input conjugated
  //    in1 = first input buffer for dot-product
  //    in2 = second input buffer for dot-product, conjugated
  //    len = number of elements to compute dot-product over
  //    EVM::dotpc = dot-product of in1 and in2
  inline complex<float> dotpc(complex<float> const * const in1,
      complex<float> const * const in2, const int len);
  inline complex<float> dotpc(float const * const in1,
      complex<float> const * const in2, const int len);


  //================//
  // Raise to Power //
  //================//

  //=============================//
  // Magnitude/Magnitude Squared //
  //=============================//

  //=======================//
  // Specialized Functions //
  //=======================//

  //================================//
  // Helper Function Implementation //
  //================================//
  template<typename T>
  inline size_t getSize(const size_t desired)
  {
    static_assert(
        std::is_same<T, float>::value ||
        std::is_same<T, complex<float> >::value,
        "[EVM::getSize] Cannot use evm.h with types other than float or "
        "complex<float>");
    if (MAX_SIMD_REG == 0) return desired; // all SIMD disabled
    int regSize; // zmm, ymm, or xmm register size in bytes
    if (__builtin_cpu_supports("avx512f"))
    {
      regSize = 64;
    }
    else if (__builtin_cpu_supports("avx2") || __builtin_cpu_supports("avx"))
    {
      regSize = 32;
    }
    else // everything else use 16 as pretty much every x86 supports SSE
    {
      regSize = 16;
    }
    // decrease size if certain instruction sets were explicitly disabled
    regSize = std::min(regSize, MAX_SIMD_REG);
    // if we're already an integer multiple of the simd register size, return
    if((desired*sizeof(T))%regSize == 0) return desired;
    // note: this logic requires that regSize is an integer multiple of
    // sizeof(T), i.e. that an integer number of T objects fit in a simd
    // register (this is true for float, double, complex<float>, etc.),
    // this is guaranteed by the static assert above
    return desired+((regSize-((desired*sizeof(T))%regSize))/sizeof(T));
  }
  
  // Returns the highest optimization level as a string, only provide
  // implementation for each target if it wasn't explicitly disabled,
  // note this must return a char* vs. a string due to a quirk of the intel
  // compiler which will not instantiate string inside of multiversioning
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline const char* getOptLevelRaw() {return "AVX512";}
  #endif
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2")))
  inline const char* getOptLevelRaw() {return "AVX2";}
  #endif
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline const char* getOptLevelRaw() {return "AVX";}
  #endif
  #if !defined(DISABLE_SSE4)
  __attribute__((__target__("sse4.2")))
  inline const char* getOptLevelRaw() {return "SSE4.2";}
  __attribute__((__target__("sse4.1")))
  inline const char* getOptLevelRaw() {return "SSE4.1";}
  #endif
  #if !defined(DISABLE_SSSE3)
  __attribute__((__target__("ssse3")))
  inline const char* getOptLevelRaw() {return "SSSE3";}
  #endif
  #if !defined(DISABLE_SSE3)
  __attribute__((__target__("sse3")))
  inline const char* getOptLevelRaw() {return "SSE3";}
  #endif
  #if !defined(DISABLE_SSE2)
  __attribute__((__target__("sse2")))
  inline const char* getOptLevelRaw() {return "SSE2";}
  #endif
  #if !defined(DISABLE_SSE)
  __attribute__((__target__("sse")))
  inline const char* getOptLevelRaw() {return "SSE";}
  #endif
  __attribute__((__target__("default")))
  inline const char* getOptLevelRaw() {return "SIMD Disabled";}

  inline std::string getOptLevel() {return std::string(getOptLevelRaw());}
  
  //======================================//
  // Vector Scale Function Implementation //
  //======================================//
  // real data, real scale factor
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void scale(float const * const in1, const float scaleFactor,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld;
    __m512 sc = _mm512_set1_ps(scaleFactor);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld = _mm512_loadu_ps(&in1[i]);
      ld = _mm512_mul_ps(ld, sc);
      _mm512_storeu_ps(&out[i], ld);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld = _mm512_mul_ps(ld, sc);
      _mm512_mask_storeu_ps(&out[i], mk, ld);
    }
    return;
  }
  #endif // AVX512 real x real scale
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void scale(float const * const in1, const float scaleFactor,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld;
    __m256 sc = _mm256_set1_ps(scaleFactor);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld = _mm256_loadu_ps(&in1[i]);
      ld = _mm256_mul_ps(ld, sc);
      _mm256_storeu_ps(&out[i], ld);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld = _mm256_maskload_ps(&in1[i], msk);
      ld = _mm256_mul_ps(ld, sc);
      _mm256_maskstore_ps(&out[i], msk, ld);
    }
    return;
  }
  #endif // AVX real x real scale
  __attribute__((__target__("default")))
  inline void scale(float const * const in1, const float scaleFactor,
      const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = scaleFactor*in1[i];
    return;
  }

  // Complex data, real scale factor
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void scale(complex<float> const * const in1, const float scaleFactor,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld;
    __m512 sc = _mm512_set1_ps(scaleFactor);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld = _mm512_mul_ps(ld, sc);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = mask16(rem<<1);
      ld = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i]));
      ld = _mm512_mul_ps(ld, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld);
    }
    return;
  }
  #endif // end AVX512 complex x real scale
  #if !defined(DISABLE_AVX) // AVX complex x real scale
  __attribute__((__target__("avx")))
  inline void scale(complex<float> const * const in1, const float scaleFactor,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld;
    __m256 sc = _mm256_set1_ps(scaleFactor);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld = _mm256_mul_ps(ld, sc);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld = _mm256_mul_ps(ld, sc);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld);
    }
    return;
  }
  #endif // AVX512 complex x real scale
  __attribute__((__target__("default")))
  inline void scale(complex<float> const * const in1, const float scaleFactor,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*scaleFactor;
  }

  // real data, complex scale factor
  #if !defined(DISABLE_AVX512) // AVX512 real data x complex scale
  __attribute__((__target__("avx512f")))
  inline void scale(float const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3;
    const __m512 re = _mm512_set1_ps(scaleFactor.real());
    const __m512 im = _mm512_set1_ps(scaleFactor.imag());
    const __m512i p1 = _mm512_setr_epi32(0,16,1,17,2,18,3,19,4,20,5,21,6,22,
        7,23);
    const __m512i p2 = _mm512_setr_epi32(8,24,9,25,10,26,11,27,12,28,13,29,
        14,30,15,31);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_mul_ps(ld1, re);
      ld3 = _mm512_mul_ps(ld1, im);
      ld1 = _mm512_permutex2var_ps(ld2, p1, ld3);
      ld2 = _mm512_permutex2var_ps(ld2, p2, ld3);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld2);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld2 = _mm512_mul_ps(ld1, re);
      ld3 = _mm512_mul_ps(ld1, im);
      ld1 = _mm512_permutex2var_ps(ld2, p1, ld3);
      ld2 = _mm512_permutex2var_ps(ld2, p2, ld3);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]),
          mask16((rem-8)<<1), ld2);
    }
    else if(rem)
    {
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld2 = _mm512_mul_ps(ld1, re);
      ld3 = _mm512_mul_ps(ld1, im);
      ld1 = _mm512_permutex2var_ps(ld2, p1, ld3);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]),
          mask16(rem<<1), ld1);
    }
    return;
  }
  #endif // end AVX512 real x complex scale
  #if !defined(DISABLE_AVX) // AVX real x complex scale
  __attribute__((__target__("avx")))
  inline void scale(float const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3;
    __m256 re = _mm256_set1_ps(scaleFactor.real());
    __m256 im = _mm256_set1_ps(scaleFactor.imag());
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_mul_ps(ld1, re); // [a0,a1,a2,a3,a4,a5,a6,a7]
      ld3 = _mm256_mul_ps(ld1, im); // [b0,b1,b2,b3,b4.b5,b6,b7]
      ld1 = _mm256_unpacklo_ps(ld2, ld3); // [a0,b0,a1,b1,a4,b4,a5,b5]
      ld2 = _mm256_unpackhi_ps(ld2, ld3); // [a2,b2,a3,b3,a6,b6,a7,b7]
      ld3 = _mm256_permute2f128_ps(ld1, ld2, 0x20); // 
      ld1 = _mm256_permute2f128_ps(ld1, ld2, 0x31); // 
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld3);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem>4) // if remainder is > 4, need 2 registers worth
    {
      // 2 floats per complex, so double rem for msk2
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      const __m256i msk2 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      ld2 = _mm256_mul_ps(ld1, re); // [r0,r1,r2,r3,r4,r5,r6,r7]
      ld3 = _mm256_mul_ps(ld1, im); // [i0,i1,i2,i3,i4.i5,i6,i7]
      ld1 = _mm256_unpacklo_ps(ld2, ld3); // [r0,i0,r1,i1,r4,i4,r5,i5]
      ld2 = _mm256_unpackhi_ps(ld2, ld3); // [r2,i2,r3,i3,r6,i6,r7,i7]
      ld3 = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [r0,i0,r1,i1,r2,i2,r3,i3]
      ld1 = _mm256_permute2f128_ps(ld1, ld2, 0x31); // [r4,i4,r5,i5,r6,i6,r7,i7]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld3);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
          ld1);
    }
    else if(rem)
    {
      // 2 floats per complex, so double rem for msk2
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      const __m256i msk2 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      ld2 = _mm256_mul_ps(ld1, re); // [r0,r1,r2,r3,r4,r5,r6,r7]
      ld3 = _mm256_mul_ps(ld1, im); // [i0,i1,i2,i3,i4.i5,i6,i7]
      ld1 = _mm256_unpacklo_ps(ld2, ld3); // [r0,i0,r1,i1,r4,i4,r5,i5]
      ld2 = _mm256_unpackhi_ps(ld2, ld3); // [r2,i2,r3,i3,r6,i6,r7,i7]
      ld3 = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [r0,i0,r1,i1,r2,i2,r3,i3]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
          ld3);
    }
    return;
  }
  #endif // end AVX real x complex scale
  __attribute__((__target__("default"))) // default real x complex scale
  inline void scale(float const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = scaleFactor*in1[i];
    return;
  }

  // Complex data, complex scale factor
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void scale(complex<float> const * const in1,
      const complex<float> scaleFactor, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    __m512 re = _mm512_set1_ps(scaleFactor.real());
    __m512 im = _mm512_set1_ps(scaleFactor.imag());
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      ld2 = _mm512_mul_ps(ld2, im); // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = mask16(rem<<1);
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i]));
      ld2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      ld2 = _mm512_mul_ps(ld2, im); // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x complex scale
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void scale(complex<float> const * const in1,
      const complex<float> scaleFactor, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2;
    __m256 re = _mm256_set1_ps(scaleFactor.real());
    __m256 im = _mm256_set1_ps(scaleFactor.imag());
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1);
      ld1 = _mm256_mul_ps(ld1, re);
      ld2 = _mm256_mul_ps(ld2, im);
      ld1 = _mm256_addsub_ps(ld1, ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1);
      ld1 = _mm256_mul_ps(ld1, re);
      ld2 = _mm256_mul_ps(ld2, im);
      ld1 = _mm256_addsub_ps(ld1, ld2);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX complex x complex scale
  __attribute__((__target__("default")))
  inline void scale(complex<float> const * const in1,
      const complex<float> scaleFactor, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*scaleFactor;
    return;
  }

  //=========================================//
  // Vector Multiply Function Implementation //
  //=========================================//
  // real x real
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void mult(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld1 = _mm512_mul_ps(ld1, ld2);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld1 = _mm512_mul_ps(ld1, ld2);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 real x real
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void mult(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_mul_ps(ld1, ld2);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld2);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX real x real
  __attribute__((__target__("default")))
  inline void mult(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i];
  }

  // complex x complex
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void mult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = mask16(rem<<1);
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in1[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in2[i])); // B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x complex
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void mult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX complex x complex
  __attribute__((__target__("default")))
  inline void mult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i];
    return;
  }
  
  // real x complex
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))  // AVX512 real x complex
  inline void mult(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 real x complex
  #if !defined(DISABLE_AVX2) // AVX2 real x complex 
  __attribute__((__target__("avx2")))
  inline void mult(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX2 real x complex
  #if !defined(DISABLE_AVX) // AVX real x complex
  __attribute__((__target__("avx")))
  inline void mult(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX real x complex
  __attribute__((__target__("default"))) // default real x complex
  inline void mult(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i];
  } // end default real x complex

  // complex x conj(complex)
  #if !defined(DISABLE_AVX512) // AVX512 complex x conj(complex)
  __attribute__((__target__("avx512f")))
  inline void multc(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = mask16(rem<<1);
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in1[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in2[i])); // B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x conj(complex)
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void multc(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex x conj(complex)
  #if !defined(DISABLE_AVX) // AVX complex x conj(complex)
  __attribute__((__target__("avx")))
  inline void multc(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex x conj(complex)
  __attribute__((__target__("default")))
  inline void multc(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    // written out rather than in1[i]*conj(in2[i]), see the icpx note above
    for(int i=0; i<len; ++i)
    {
      const float re = in1[i].real()*in2[i].real() +
          in1[i].imag()*in2[i].imag();
      const float im = in1[i].imag()*in2[i].real() -
          in1[i].real()*in2[i].imag();
      out[i].real(re);
      out[i].imag(im);
    }
    return;
  }
  
  // real x conj(complex)
  #if !defined(DISABLE_AVX512) // AVX512 real x conj(complex)
  __attribute__((__target__("avx512f")))
  inline void multc(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      // negate every other element, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_xor_ps (casting should be free)
      sc = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc),neg));
      ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      sc = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc),neg));
      ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      sc = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc),neg));
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 real x conj(complex)
  #if !defined(DISABLE_AVX2) // AVX2 real x conj(complex)
  __attribute__((__target__("avx2")))
  inline void multc(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      sc = _mm256_xor_ps(sc, neg);    // negate every other element
      ld1 = _mm256_xor_ps(ld1, neg);  // negate every other element
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        sc = _mm256_xor_ps(sc, neg);    // negate every other element
        ld1 = _mm256_xor_ps(ld1, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        sc = _mm256_xor_ps(sc, neg);    // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // AVX2 real x conj(complex)
  #if !defined(DISABLE_AVX) // AVX real x conj(complex)
  __attribute__((__target__("avx")))
  inline void multc(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      sc1 = _mm256_xor_ps(sc1, neg);      // negate every other element
      sc2 = _mm256_xor_ps(sc2, neg);      // negate every other element
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        sc1 = _mm256_xor_ps(sc1, neg);      // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);      // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc1 = _mm256_xor_ps(sc1, neg);        // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX real x conj(complex)
  __attribute__((__target__("default"))) // default real x conj(complex)
  inline void multc(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    // written out rather than in1[i]*conj(in2[i]), see the icpx note above
    for(int i=0; i<len; ++i)
    {
      const float re = in1[i]*in2[i].real();
      const float im = -in1[i]*in2[i].imag();
      out[i].real(re);
      out[i].imag(im);
    }
    return;
  }
  
  //====================================================//
  // Vectory Multiply and Scale Function Implementation //
  //====================================================//
  // real x real x scale
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void mults(float const * const in1, float const * const in2,
      const int len, const float scale, float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    const __m512 sc = _mm512_set1_ps(scale);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld1 = _mm512_mul_ps(ld1, ld2);
      ld1 = _mm512_mul_ps(ld1, sc);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld1 = _mm512_mul_ps(ld1, ld2);
      ld1 = _mm512_mul_ps(ld1, sc);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 real x real x scale
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void mults(float const * const in1, float const * const in2,
      const int len, const float scale, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2;
    const __m256 sc = _mm256_set1_ps(scale);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_mul_ps(ld1, ld2);
      ld1 = _mm256_mul_ps(ld1, sc);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld2);
      ld1 = _mm256_mul_ps(ld1, sc);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX real x real x scale
  __attribute__((__target__("default"))) // default real x real x scale
  inline void mults(float const * const in1, float const * const in2,
      const int len, const float scale, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i]*scale;
  }

  // complex x complex x scale
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void mults(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, sh, re, im;
    const __m512 sc = _mm512_set1_ps(scale);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld2 = _mm512_mul_ps(ld2, sc);
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = mask16(rem<<1);
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in1[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in2[i])); // B
      ld2 = _mm512_mul_ps(ld2, sc);
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x complex x scale
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void mults(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    const __m256 sc = _mm256_set1_ps(scale);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld2 = _mm256_mul_ps(ld2, sc);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      ld2 = _mm256_mul_ps(ld2, sc);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX complex x complex x scale
  __attribute__((__target__("default"))) // default complex x complex x scale
  inline void mults(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i]*scale;
  }

  // real x complex x scale
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void mults(float const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    const __m512 sca = _mm512_set1_ps(scale);
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld1 = _mm512_mul_ps(ld1, sca);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld1 = _mm512_mul_ps(ld1, sca);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld1 = _mm512_mul_ps(ld1, sca);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 real x complex x scale
  #if !defined(DISABLE_AVX2) // AVX2 real x complex x scale
  __attribute__((__target__("avx2")))
  inline void mults(float const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    const __m256 sca = _mm256_set1_ps(scale);
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld1 = _mm256_mul_ps(ld1, sca);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      ld1 = _mm256_mul_ps(ld1, sca);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX2 real x complex x scale
  #if !defined(DISABLE_AVX) // AVX real x complex x scale
  __attribute__((__target__("avx")))
  inline void mults(float const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2;
    const __m256 sca = _mm256_set1_ps(scale);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld1 = _mm256_mul_ps(ld1, sca);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      ld1 = _mm256_mul_ps(ld1, sca);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX real x complex x scale
  __attribute__((__target__("default"))) // default real x complex x scale
  inline void mults(float const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i]*scale;
  }

  // complex x conj(complex) x scale
  #if !defined(DISABLE_AVX512) // AVX512 complex x conj(complex) x scale
  __attribute__((__target__("avx512f")))
  inline void multcs(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, sh, re, im;
    const __m512 sc = _mm512_set1_ps(scale);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld2 = _mm512_mul_ps(ld2, sc);
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = mask16(rem<<1);
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in1[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in2[i])); // B
      ld2 = _mm512_mul_ps(ld2, sc);
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x conj(complex) x scale
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void multcs(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    const __m256 sc = _mm256_set1_ps(scale);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld2 = _mm256_mul_ps(ld2, sc);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      ld2 = _mm256_mul_ps(ld2, sc);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex x conj(complex) x scale
  #if !defined(DISABLE_AVX) // AVX complex x conj(complex) x scale
  __attribute__((__target__("avx")))
  inline void multcs(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    const __m256 sc = _mm256_setr_ps(scale, -scale, scale, -scale, scale,
        -scale, scale, -scale);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_mul_ps(ld2, sc);// conj(B) x scale
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_mul_ps(ld2, sc);// conj(B) x scale
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex x conj(complex) x scale
  // default complex x conj(complex) x scale
  __attribute__((__target__("default")))
  inline void multcs(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out)
  {
    // written out rather than in1[i]*conj(in2[i])*scale, see the icpx note
    // above
    for(int i=0; i<len; ++i)
    {
      const float re = (in1[i].real()*in2[i].real() +
          in1[i].imag()*in2[i].imag())*scale;
      const float im = (in1[i].imag()*in2[i].real() -
          in1[i].real()*in2[i].imag())*scale;
      out[i].real(re);
      out[i].imag(im);
    }
  }

  //=======================================//
  // Vector Divide Function Implementation //
  //=======================================//
  // real / real
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void div(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld1 = _mm512_div_ps(ld1, ld2);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld1 = _mm512_maskz_div_ps(mk, ld1, ld2);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 real / real
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void div(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_div_ps(ld1, ld2);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      ld1 = _mm256_div_ps(ld1, ld2);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX real / real
  __attribute__((__target__("default"))) // default real / real
  inline void div(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }

  // complex / complex
  #if !defined(DISABLE_AVX512) // AVX512 complex / complex
  __attribute__((__target__("avx512f")))
  inline void div(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, ld4, sh, re, im, rc;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 8)
    {
      const __m512 one = _mm512_set1_ps(1.0f);
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld4 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm512_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in1[i+8]));
        ld4 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i]+8));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load A and B and compute A/B
    ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i])); // A
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i])); // B
    sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
    im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
    re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
    ld2 = _mm512_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br7)^2,(Bi7)^2]
    im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
    sh = _mm512_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
    ld2 = _mm512_add_ps(ld2, sh); // |B|^2
    ld1 = _mm512_maskz_div_ps(mk, ld1, ld2);// (A*conj(B))/|B|^2 == A/B
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    return;
  }
  #endif // AVX512 complex / complex
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void div(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm256_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm256_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute A/B
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    ld1 = _mm256_div_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX2 complex / complex
  #if !defined(DISABLE_AVX) // AVX complex / complex
  __attribute__((__target__("avx")))
  inline void div(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // conjugate B and D
        ld2 = _mm256_xor_ps(ld2, neg); // conj(B) = [Br0,-Bi0,...,Br3,-Bi3]
        ld4 = _mm256_xor_ps(ld4, neg); // conj(D) = [Dr0,-Di0,...,Dr3,-Di3]
        // compute complex numerator A*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [-Bi0,-Bi0,-Bi1,-Bi1,...,-Bi3,-Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,-Ar0*Bi0,...,-Ai3*Bi3,-Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [-Di0,-Di0,-Di1,-Di1,...,-Di3,-Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [-Ci0*Di0,-Cr0*Di0,...,-Ci3*Di3,-Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm256_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // conjugate B and D
        ld2 = _mm256_xor_ps(ld2, neg); // conj(B) = [Br0,-Bi0,...,Br3,-Bi3]
        ld4 = _mm256_xor_ps(ld4, neg); // conj(D) = [Dr0,-Di0,...,Dr3,-Di3]
        // compute complex numerator A*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [-Bi0,-Bi0,-Bi1,-Bi1,...,-Bi3,-Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,-Ar0*Bi0,...,-Ai3*Bi3,-Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [-Di0,-Di0,-Di1,-Di1,...,-Di3,-Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [-Ci0*Di0,-Cr0*Di0,...,-Ci3*Di3,-Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm256_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute A/B
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    ld2 = _mm256_xor_ps(ld2, neg); // conj(B) = [Br0,-Bi0,...,Br3,-Bi3]
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [-Bi0,-Bi0,-Bi1,-Bi1,...,-Bi3,-Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,-Ar0*Bi0,...,-Ai3*Bi3,-Ar3*Bi3]
    ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
    ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    ld1 = _mm256_div_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX complex / complex
  __attribute__((__target__("default"))) // default complex / complex
  inline void div(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }

  // real / complex
  #if !defined(DISABLE_AVX512) // AVX512 real / complex
  __attribute__((__target__("avx512f")))
  inline void div(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    // permutation constants to copy each real element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 8)
    {
      const __m512 one = _mm512_set1_ps(1.0f);
      for(; i<len-15; i+=16) // process 16 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_loadu_ps(&in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        // negate every other element, note the casting is done to avoid
        // requiring the AVX512DQ instruction _mm512_xor_ps (casting should be
        // free)
        sc1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc1),neg));
        sc2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc2),neg));
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * conj(in2[i+8])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm512_mul_ps(ld3, rc); // (in1*conj(in2[i+8]))/|in2[i+8]|^2
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        // negate every other element, note the casting is done to avoid
        // requiring the AVX512DQ instruction _mm512_xor_ps (casting should be
        // free)
        sc1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc1),neg));
        sc2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc2),neg));
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * conj(in2[i+8])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        // (in1*conj(in2[i+8]))/|in2[i+8]|^2
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load in1 and in2 and compute in1/in2
    ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i]));
    sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
    // negate every other element, note the casting is done to avoid requiring
    // the AVX512DQ instruction _mm512_xor_ps (casting should be free)
    sc1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc1),neg));
    ld1 = _mm512_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
    sc2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm512_add_ps(sc2, ld1);  // |in2[i]|^2
    ld2 = _mm512_maskz_div_ps(mk, ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    return;
  }
  #endif // end AVX512 real / complex
  #if !defined(DISABLE_AVX2) // AVX2 real / complex
  __attribute__((__target__("avx2,fma")))
  inline void div(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 8 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        // note msk2 accounts for 2 reals per element for the complex buffer
        const __m256i msk1 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem]));
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_maskload_ps(&in1[i], msk1);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // note msk2 accounts for 2 reals per element for the complex buffer
    const __m256i msk1 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem]));
    const __m256i msk2 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load in1 and in2 and compute in1/in2
    ld1 = _mm256_maskload_ps(&in1[i], msk1);
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
        msk2);
    sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
    sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
    ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
    sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
    ld2 = _mm256_div_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
    return;
  }
  #endif // AVX2 real / complex
  #if !defined(DISABLE_AVX) // AVX real / complex
  __attribute__((__target__("avx")))
  inline void div(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 8 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        // note msk2 accounts for 2 reals per element for the complex buffer
        const __m256i msk1 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem]));
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_maskload_ps(&in1[i], msk1);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // note msk2 accounts for 2 reals per element for the complex buffer
    const __m256i msk1 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem]));
    const __m256i msk2 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load in1 and in2 and compute in1/in2
    ld1 = _mm256_maskload_ps(&in1[i], msk1);
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
        msk2);
    sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
    ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
    sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
    sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
    ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
    sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
    ld2 = _mm256_div_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
    return;
  }
  #endif // end AVX real / complex
  __attribute__((__target__("default"))) // default real / complex
  inline void div(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    // complex(in1[i], 0) without the complex<float>(re, im) constructor, see
    // the icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> a = zero; // imaginary part stays 0
    for(int i=0; i<len; ++i)
    {
      a.real(in1[i]);
      out[i] = a/in2[i];
    }
  }

  // complex / real
  #if !defined(DISABLE_AVX512) // AVX512 complex / real
  __attribute__((__target__("avx512f")))
  inline void div(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    const __m512 one = _mm512_set1_ps(1.0f);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      // compute 1/in2 with a single divide, then spread it to each complex
      ld1 = _mm512_div_ps(one, ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      // compute 1/in2 with a single divide, then spread it to each complex
      ld1 = _mm512_maskz_div_ps(mask16(rem), one, ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_maskz_mul_ps(mk, ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_maskz_div_ps(mk, ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 complex / real
  #if !defined(DISABLE_AVX2) // AVX2 complex / real
  __attribute__((__target__("avx2")))
  inline void div(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256 one = _mm256_set1_ps(1.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // compute 1/in2 with a single divide, then spread it to each complex
      ld1 = _mm256_div_ps(one, ld1);
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        // compute 1/in2 with a single divide, then spread it to each complex
        ld1 = _mm256_div_ps(one, ld1);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_div_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // AVX2 complex / real
  #if !defined(DISABLE_AVX) // AVX complex / real
  __attribute__((__target__("avx")))
  inline void div(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2;
    const __m256 one = _mm256_set1_ps(1.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // compute 1/in2 with a single divide, then spread it to each complex
      ld1 = _mm256_div_ps(one, ld1);
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        // compute 1/in2 with a single divide, then spread it to each complex
        ld1 = _mm256_div_ps(one, ld1);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_div_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX complex / real
  __attribute__((__target__("default"))) // default complex / real
  inline void div(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  } // end default real x complex

  // real / real with reciprocal approx and 1 stage Newton-Raphson refinement
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void divnr(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, inv;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m512 two = _mm512_set1_ps(2.0f);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld1 = _mm512_loadu_ps(&in1[i]);
      inv = _mm512_rcp14_ps(ld2);   // old
      // use one newton-raphson stage (new = old*(2-orig*old)) where 'new' is
      // updated reciprocal approx, 'old' is the current reciprocal approx,
      // and 'orig' is the original value we're trying to invert
      ld2 = _mm512_fnmadd_ps(ld2, inv, two);  // 2-orig*old
      ld2 = _mm512_mul_ps(ld2, inv);          // old*(2-orig*old)
      ld1 = _mm512_mul_ps(ld1, ld2);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      inv = _mm512_maskz_rcp14_ps(mk, ld2);   // old
      ld2 = _mm512_fnmadd_ps(ld2, inv, two);  // 2-orig*old
      ld2 = _mm512_mul_ps(inv, ld2);          // old*(2-orig*old)
      ld1 = _mm512_mul_ps(ld1, ld2);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // end AVX512 real * recip(real) with NR
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void divnr(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, inv;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      inv = _mm256_rcp_ps(ld2); // old
      // use one newton-raphson stage (new = old*(2-orig*old)) where 'new' is
      // updated reciprocal approx, 'old' is the current reciprocal approx,
      // and 'orig' is the original value we're trying to invert
      ld2 = _mm256_fnmadd_ps(ld2, inv, two);  // 2-orig*old
      ld2 = _mm256_mul_ps(inv, ld2);          // old*(2-orig*old)
      ld1 = _mm256_mul_ps(ld1, ld2);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      inv = _mm256_rcp_ps(ld2);             // old
      ld2 = _mm256_fnmadd_ps(ld2, inv, two);  // 2-orig*old
      ld2 = _mm256_mul_ps(inv, ld2);          // old*(2-orig*old)
      ld1 = _mm256_mul_ps(ld1, ld2);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // end AVX2 real * recip(real) with NR
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void divnr(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, inv;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      inv = _mm256_rcp_ps(ld2);   // old
      // use one newton-raphson stage (new = old*(2-orig*old)) where 'new' is
      // updated reciprocal approx, 'old' is the current reciprocal approx,
      // and 'orig' is the original value we're trying to invert
      ld2 = _mm256_mul_ps(ld2, inv); // orig*old
      ld2 = _mm256_sub_ps(two, ld2); // 2-orig*old
      ld2 = _mm256_mul_ps(inv, ld2); // old*(2-orig*old)
      ld1 = _mm256_mul_ps(ld1, ld2);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      inv = _mm256_rcp_ps(ld2);    // old
      ld2 = _mm256_mul_ps(ld2, inv); // orig*old
      ld2 = _mm256_sub_ps(two, ld2); // 2-orig*old
      ld2 = _mm256_mul_ps(inv, ld2); // old*(2-orig*old)
      ld1 = _mm256_mul_ps(ld1, ld2);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // end AVX real * recip(real) with NR
  __attribute__((__target__("default"))) // default real * recip(real) with NR
  inline void divnr(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }

  // complex / complex with reciprocal approximation and 1 stage Newton-Raphson
  #if !defined(DISABLE_AVX512) // AVX512 complex * recip(complex) with NR
  __attribute__((__target__("avx512f")))
  inline void divnr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, ld4, sh, re, im, rc;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m512 two = _mm512_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld4 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm512_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in1[i+8]));
        ld4 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i]+8));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load A and B and compute A/B
    ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i])); // A
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i])); // B
    sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
    im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
    re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
    ld2 = _mm512_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br7)^2,(Bi7)^2]
    im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
    sh = _mm512_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
    ld2 = _mm512_add_ps(ld2, sh); // |B|^2
    // compute 1/|B|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    re = _mm512_maskz_rcp14_ps(mk, ld2);    // inv
    ld2 = _mm512_fnmadd_ps(ld2, re, two); // 2-orig*inv
    ld2 = _mm512_mul_ps(re, ld2); // inv*(2-orig*inv) == 1/|B|^2
    ld1 = _mm512_mul_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    return;
  }
  #endif // AVX512 complex * recip(complex) with NR
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void divnr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm256_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm256_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute A/B
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    // compute 1/|B|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    re = _mm256_rcp_ps(ld2);    // inv
    ld2 = _mm256_fnmadd_ps(ld2, re, two); // 2-orig*inv
    ld2 = _mm256_mul_ps(re, ld2); // inv*(2-orig*inv) == 1/|B|^2
    ld1 = _mm256_mul_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX2 complex * recip(complex) with NR
  #if !defined(DISABLE_AVX) // AVX complex * recip(complex) with NR
  __attribute__((__target__("avx")))
  inline void divnr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // conjugate B and D
        ld2 = _mm256_xor_ps(ld2, neg); // conj(B) = [Br0,-Bi0,...,Br3,-Bi3]
        ld4 = _mm256_xor_ps(ld4, neg); // conj(D) = [Dr0,-Di0,...,Dr3,-Di3]
        // compute complex numerator A*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [-Bi0,-Bi0,-Bi1,-Bi1,...,-Bi3,-Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,-Ar0*Bi0,...,-Ai3*Bi3,-Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [-Di0,-Di0,-Di1,-Di1,...,-Di3,-Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [-Ci0*Di0,-Cr0*Di0,...,-Ci3*Di3,-Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm256_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // conjugate B and D
        ld2 = _mm256_xor_ps(ld2, neg); // conj(B) = [Br0,-Bi0,...,Br3,-Bi3]
        ld4 = _mm256_xor_ps(ld4, neg); // conj(D) = [Dr0,-Di0,...,Dr3,-Di3]
        // compute complex numerator A*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [-Bi0,-Bi0,-Bi1,-Bi1,...,-Bi3,-Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,-Ar0*Bi0,...,-Ai3*Bi3,-Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [-Di0,-Di0,-Di1,-Di1,...,-Di3,-Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [-Ci0*Di0,-Cr0*Di0,...,-Ci3*Di3,-Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm256_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute A/B
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    ld2 = _mm256_xor_ps(ld2, neg); // conj(B) = [Br0,-Bi0,...,Br3,-Bi3]
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [-Bi0,-Bi0,-Bi1,-Bi1,...,-Bi3,-Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,-Ar0*Bi0,...,-Ai3*Bi3,-Ar3*Bi3]
    ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
    ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    // compute 1/|B|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    re = _mm256_rcp_ps(ld2);    // inv
    ld2 = _mm256_mul_ps(ld2, re); // orig*inv
    ld2 = _mm256_sub_ps(two, ld2);// 2-orig*inv
    ld2 = _mm256_mul_ps(re, ld2); // inv*(2-orig*inv) == 1/|B|^2
    ld1 = _mm256_mul_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX complex * recip(complex) with NR
  __attribute__((__target__("default"))) // default complex * recip(complex) NR
  inline void divnr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }

  // real / complex using reciprocal approximation and 1 stage Newton-Raphson
  #if !defined(DISABLE_AVX512) // AVX512 real * recip(complex) with NR
  __attribute__((__target__("avx512f")))
  inline void divnr(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    // permutation constants to copy each real element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m512 two = _mm512_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 16 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_loadu_ps(&in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        // negate every other element, note the casting is done to avoid
        // requiring the AVX512DQ instruction _mm512_xor_ps (casting should be
        // free)
        sc1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc1),neg));
        sc2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc2),neg));
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * conj(in2[i+8])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm512_mul_ps(ld3, rc); // (in1*conj(in2[i+8]))/|in2[i+8]|^2
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        // negate every other element, note the casting is done to avoid
        // requiring the AVX512DQ instruction _mm512_xor_ps (casting should be
        // free)
        sc1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc1),neg));
        sc2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc2),neg));
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * conj(in2[i+8])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        // (in1*conj(in2[i+8]))/|in2[i+8]|^2
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load in1 and in2 and compute in1/in2
    ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i]));
    sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
    // negate every other element, note the casting is done to avoid requiring
    // the AVX512DQ instruction _mm512_xor_ps (casting should be free)
    sc1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc1),neg));
    ld1 = _mm512_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
    sc2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm512_add_ps(sc2, ld1);  // |in2[i]|^2
    // compute 1/|in2[i]|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    sh = _mm512_maskz_rcp14_ps(mk, sc2);    // inv
    sc2 = _mm512_fnmadd_ps(sc2, sh, two); // 2-orig*inv
    sc2 = _mm512_mul_ps(sh, sc2); // inv*(2-orig*inv) == 1/|in2[i]|^2
    ld2 = _mm512_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    return;
  }
  #endif // end AVX512 real * recip(complex) with NR
  #if !defined(DISABLE_AVX2) // AVX2 real * recip(complex) with NR
  __attribute__((__target__("avx2,fma")))
  inline void divnr(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 8 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        // note msk2 accounts for 2 reals per element for the complex buffer
        const __m256i msk1 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem]));
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_maskload_ps(&in1[i], msk1);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // note msk2 accounts for 2 reals per element for the complex buffer
    const __m256i msk1 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem]));
    const __m256i msk2 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load in1 and in2 and compute in1/in2
    ld1 = _mm256_maskload_ps(&in1[i], msk1);
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
        msk2);
    sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
    sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
    ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
    sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
    // compute 1/|in2[i]|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    sh = _mm256_rcp_ps(sc2);    // inv
    sc2 = _mm256_fnmadd_ps(sc2, sh, two); // 2-orig*inv
    sc2 = _mm256_mul_ps(sh, sc2); // inv*(2-orig*inv) == 1/|in2[i]|^2
    ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
    return;
  }
  #endif // AVX2 real * recip(complex) with NR
  #if !defined(DISABLE_AVX) // AVX real * recip(complex) with NR
  __attribute__((__target__("avx")))
  inline void divnr(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 8 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        // note msk2 accounts for 2 reals per element for the complex buffer
        const __m256i msk1 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem]));
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_maskload_ps(&in1[i], msk1);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // note msk2 accounts for 2 reals per element for the complex buffer
    const __m256i msk1 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem]));
    const __m256i msk2 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load in1 and in2 and compute in1/in2
    ld1 = _mm256_maskload_ps(&in1[i], msk1);
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
        msk2);
    sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
    ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
    sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
    sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
    ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
    sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
    // compute 1/|in2[i]|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    sh = _mm256_rcp_ps(sc2);    // inv
    sc2 = _mm256_mul_ps(sc2, sh); // orig*inv
    sc2 = _mm256_sub_ps(two, sc2);// 2-orig*inv
    sc2 = _mm256_mul_ps(sh, sc2); // inv*(2-orig*inv) == 1/|in2[i]|^2
    ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
    return;
  }
  #endif // end AVX real * recip(complex) with NR
  __attribute__((__target__("default"))) // default real*recip(complex) with NR
  inline void divnr(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    // complex(in1[i], 0) without the complex<float>(re, im) constructor, see
    // the icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> a = zero; // imaginary part stays 0
    for(int i=0; i<len; ++i)
    {
      a.real(in1[i]);
      out[i] = a/in2[i];
    }
  }

  // complex / real using reciprocal approximation and 1 stage Newton-Raphson
  #if !defined(DISABLE_AVX512) // AVX512 complex * recip(real) with NR
  __attribute__((__target__("avx512f")))
  inline void divnr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m512 two = _mm512_set1_ps(2.0f);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc = _mm512_rcp14_ps(ld1);
      ld1 = _mm512_fnmadd_ps(ld1, sc, two); // 2-orig*inv
      ld1 = _mm512_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc = _mm512_maskz_rcp14_ps(mask16(rem), ld1);
      ld1 = _mm512_fnmadd_ps(ld1, sc, two); // 2-orig*inv
      ld1 = _mm512_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc = _mm512_maskz_rcp14_ps(mask16(rem), ld1);
      ld1 = _mm512_fnmadd_ps(ld1, sc, two); // 2-orig*inv
      ld1 = _mm512_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 complex * recip(real) with NR
  #if !defined(DISABLE_AVX2) // AVX2 complex * recip(real) with NR
  __attribute__((__target__("avx2,fma")))
  inline void divnr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc = _mm256_rcp_ps(ld1);
      ld1 = _mm256_fnmadd_ps(ld1, sc, two); // 2-orig*inv
      ld1 = _mm256_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sc = _mm256_rcp_ps(ld1);
        ld1 = _mm256_fnmadd_ps(ld1, sc, two); // 2-orig*inv
        ld1 = _mm256_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sc = _mm256_rcp_ps(ld1);
        ld1 = _mm256_fnmadd_ps(ld1, sc, two); // 2-orig*inv
        ld1 = _mm256_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX2 complex * recip(real) with NR
  #if !defined(DISABLE_AVX) // AVX complex * recip(real) with NR
  __attribute__((__target__("avx")))
  inline void divnr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc1 = _mm256_rcp_ps(ld1);    // inv
      ld1 = _mm256_mul_ps(ld1, sc1); // orig*inv
      ld1 = _mm256_sub_ps(two, ld1); // 2-orig*inv
      ld1 = _mm256_mul_ps(sc1, ld1); // inv*(2-orig*inv) == 1/in2
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sc1 = _mm256_rcp_ps(ld1);    // inv
        ld1 = _mm256_mul_ps(ld1, sc1); // orig*inv
        ld1 = _mm256_sub_ps(two, ld1); // 2-orig*inv
        ld1 = _mm256_mul_ps(sc1, ld1); // inv*(2-orig*inv) == 1/in2
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sc1 = _mm256_rcp_ps(ld1);    // inv
        ld1 = _mm256_mul_ps(ld1, sc1); // orig*inv
        ld1 = _mm256_sub_ps(two, ld1); // 2-orig*inv
        ld1 = _mm256_mul_ps(sc1, ld1); // inv*(2-orig*inv) == 1/in2
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX complex * recip(real) with NR
  __attribute__((__target__("default"))) // default complex*recip(real) with NR
  inline void divnr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }

  // real / real using reciprocal approximation
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f"))) // AVX512 real * recip(real)
  inline void divr(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_rcp14_ps(ld2);
      ld1 = _mm512_mul_ps(ld1, ld2);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_rcp14_ps(mk, ld2);
      ld1 = _mm512_mul_ps(ld1, ld2);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 real * recip(real)
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx"))) // AVX real * recip(real)
  inline void divr(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_rcp_ps(ld2);
      ld1 = _mm256_mul_ps(ld1, ld2);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_rcp_ps(ld2);
      ld1 = _mm256_mul_ps(ld1, ld2);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX real * recip(real)
  __attribute__((__target__("default"))) // default real * recip(real)
  inline void divr(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }

  // complex / complex using reciprocal approximation
  #if !defined(DISABLE_AVX512) // AVX512 complex * recip(complex)
  __attribute__((__target__("avx512f")))
  inline void divr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, ld4, sh, re, im, rc;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld4 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm512_mul_ps(ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load A, B, C, and D, and compute A/B and C/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in1[i+8]));
        ld4 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i]+8));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*conj(B)
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator C*conj(D)
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0+Di0Ci0,Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*conj(B))/|B|^2 == A/B
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc); // (C*conj(D))/|D|^2 == C/D
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load A and B and compute A/B
    ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i])); // A
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i])); // B
    sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
    im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
    re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
    ld2 = _mm512_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br7)^2,(Bi7)^2]
    im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
    sh = _mm512_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
    ld2 = _mm512_add_ps(ld2, sh); // |B|^2
    ld2 = _mm512_maskz_rcp14_ps(mk, ld2);  // 1/|B|^2
    ld1 = _mm512_mul_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    return;
  }
  #endif // AVX512 complex * recip(complex)
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma"))) // AVX2 complex * recip(complex)
  inline void divr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0+Bi0Ai0,Br0*Ai0-Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex * recip(complex)
  #if !defined(DISABLE_AVX) // AVX complex * recip(complex)
  __attribute__((__target__("avx")))
  inline void divr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_addsub_ps(ld1, im);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_addsub_ps(ld1, im);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (A*conj(B))/|B|^2 == A/B
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex * recip(complex)
  __attribute__((__target__("default"))) // default complex * recip(complex)
  inline void divr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }
  
  // real / complex using reciprocal approximation
  #if !defined(DISABLE_AVX512) // AVX512 real * recip(complex)
  __attribute__((__target__("avx512f")))
  inline void divr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    // permutation constants to copy each real element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 16 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_loadu_ps(&in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        // negate every other element, note the casting is done to avoid
        // requiring the AVX512DQ instruction _mm512_xor_ps (casting should be
        // free)
        sc1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc1),neg));
        sc2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc2),neg));
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * conj(in2[i+8])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm512_mul_ps(ld3, rc); // (in1*conj(in2[i+8]))/|in2[i+8]|^2
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*conj(B) and in1*conj(D)
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        // negate every other element, note the casting is done to avoid
        // requiring the AVX512DQ instruction _mm512_xor_ps (casting should be
        // free)
        sc1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc1),neg));
        sc2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc2),neg));
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * conj(in2[i+8])
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*conj(in2[i]))/|in2[i]|^2
        // (in1*conj(in2[i+8]))/|in2[i+8]|^2
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load in1 and in2 and compute in1/in2
    ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i]));
    sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
    // negate every other element, note the casting is done to avoid requiring
    // the AVX512DQ instruction _mm512_xor_ps (casting should be free)
    sc1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(sc1),neg));
    ld1 = _mm512_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
    sc2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm512_add_ps(sc2, ld1);  // |in2[i]|^2
    sc2 = _mm512_maskz_rcp14_ps(mk, sc2);  // 1/|in2[i]|^2
    ld2 = _mm512_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    return;
  }
  #endif // end AVX512 real * recip(complex)
  #if !defined(DISABLE_AVX2) // AVX2 real * recip(complex)
  __attribute__((__target__("avx2")))
  inline void divr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2, sh;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
      sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
      sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
      sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
      ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
      ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
      sc1 = _mm256_mul_ps(ld3, ld3);  // in2[i+4]^2
      ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
      sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
      sh = _mm256_shuffle_ps(sc1, sc1, 0xb1);  // in2[i+4]^2 swap(re,im)
      sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
      sh = _mm256_add_ps(sc1, sh);    // |in2[i+4]|^2
      sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
      sh = _mm256_rcp_ps(sh);       // 1/|in2[i+4]|^2
      ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
      ld3 = _mm256_mul_ps(ld3, sh);   // (in1*conj(in2[i+4]))/|in2[i+4]|^2
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        sc1 = _mm256_mul_ps(ld3, ld3);  // in2[i+4]^2
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
        sh = _mm256_shuffle_ps(sc1, sc1, 0xb1);  // in2[i+4]^2 swap(re,im)
        sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
        sh = _mm256_add_ps(sc1, sh);    // |in2[i+4]|^2
        sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
        sh = _mm256_rcp_ps(sh);       // 1/|in2[i+4]|^2
        ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, sh);   // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
        sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
        sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
        ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // AVX2 real * recip(complex)
  #if !defined(DISABLE_AVX) // AVX real * recip(complex)
  __attribute__((__target__("avx")))
  inline void divr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2, sh;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
      sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
      ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
      ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
      sc1 = _mm256_mul_ps(ld3, ld3);  // in2[i+4]^2
      ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
      sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
      sh = _mm256_shuffle_ps(sc1, sc1, 0xb1);  // in2[i+4]^2 swap(re,im)
      sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
      sh = _mm256_add_ps(sc1, sh);    // |in2[i+4]|^2
      sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
      sh = _mm256_rcp_ps(sh);       // 1/|in2[i+4]|^2
      ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
      ld3 = _mm256_mul_ps(ld3, sh);   // (in1*conj(in2[i+4]))/|in2[i+4]|^2
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);  // negate every other element
        ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        sc1 = _mm256_mul_ps(ld3, ld3);  // in2[i+4]^2
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
        sh = _mm256_shuffle_ps(sc1, sc1, 0xb1);  // in2[i+4]^2 swap(re,im)
        sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
        sh = _mm256_add_ps(sc1, sh);    // |in2[i+4]|^2
        sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
        sh = _mm256_rcp_ps(sh);       // 1/|in2[i+4]|^2
        ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, sh);   // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc1 = _mm256_xor_ps(sc1, neg);  // negate every other element
        ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
        sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
        sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
        ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX real * recip(complex)
  __attribute__((__target__("default"))) // default real * recip(complex)
  inline void divr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    // complex(in1[i], 0) without the complex<float>(re, im) constructor, see
    // the icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> a = zero; // imaginary part stays 0
    for(int i=0; i<len; ++i)
    {
      a.real(in1[i]);
      out[i] = a/in2[i];
    }
  }
  
  // complex / real using reciprocal approximation
  #if !defined(DISABLE_AVX512) // AVX512 complex * recip(real)
  __attribute__((__target__("avx512f")))
  inline void divr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      ld1 = _mm512_rcp14_ps(ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      ld1 = _mm512_maskz_rcp14_ps(mask16(rem), ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      ld1 = _mm512_maskz_rcp14_ps(mask16(rem), ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 complex * recip(real)
  #if !defined(DISABLE_AVX2) // AVX2 complex * recip(real)
  __attribute__((__target__("avx2")))
  inline void divr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld1 = _mm256_rcp_ps(ld1);
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld1 = _mm256_rcp_ps(ld1);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        ld1 = _mm256_rcp_ps(ld1);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // AVX2 complex * recip(real)
  #if !defined(DISABLE_AVX) // AVX complex * recip(real)
  __attribute__((__target__("avx")))
  inline void divr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld1 = _mm256_rcp_ps(ld1);
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld1 = _mm256_rcp_ps(ld1);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        ld1 = _mm256_rcp_ps(ld1);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX complex * recip(real)
  __attribute__((__target__("default"))) // default complex * recip(real)
  inline void divr(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }

  // complex conjugate / complex
  #if !defined(DISABLE_AVX512) // AVX512 conj(complex) / complex
  __attribute__((__target__("avx512f")))
  inline void divc1(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, ld4, sh, re, im, rc;
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 8)
    {
      const __m512 one = _mm512_set1_ps(1.0f);
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld4 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B), note the casting is done
        // to avoid requiring the AVX512DQ instruction _mm512_xor_ps (casting
        // should be free)
        ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm512_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in1[i+8]));
        ld4 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i]+8));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B), note the casting is done
        // to avoid requiring the AVX512DQ instruction _mm512_xor_ps (casting
        // should be free)
        ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load A and B and compute conj(A)/B
    ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i])); // A
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i])); // B
    // conj(A), note the casting is done to avoid requiring the AVX512DQ
    // instruction _mm512_xor_ps (casting should be free)
    ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
    sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
    im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
    re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
    ld2 = _mm512_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br7)^2,(Bi7)^2]
    im = _mm512_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai7*Bi7,Ar7*Bi7]
    sh = _mm512_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
    ld2 = _mm512_add_ps(ld2, sh); // |B|^2
    // (conj(A)*conj(B))/|B|^2 == conj(A)/B
    ld1 = _mm512_maskz_div_ps(mk, ld1, ld2);
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    return;
  }
  #endif // AVX512 conj(complex) / complex
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void divc1(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B)
        ld1 = _mm256_xor_ps(ld1, neg); // conj(A)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm256_xor_ps(ld3, neg); // conj(C)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm256_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B)
        ld1 = _mm256_xor_ps(ld1, neg); // conj(A)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm256_xor_ps(ld3, neg); // conj(C)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm256_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute conj(A)/B
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    ld1 = _mm256_xor_ps(ld1, neg); // conj(A)
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai3*Bi3,Ar3*Bi3]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    ld1 = _mm256_div_ps(ld1, ld2);// (conj(A)*conj(B))/|B|^2 == conj(A)/B
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX2 conj(complex) / complex
  #if !defined(DISABLE_AVX) // AVX conj(complex) / complex
  __attribute__((__target__("avx")))
  inline void divc1(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A*B) == conj(A)*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        ld1 = _mm256_xor_ps(ld1, neg); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C*D) == conj(C)*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        ld3 = _mm256_xor_ps(ld3, neg); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm256_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A*B) == conj(A)*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        ld1 = _mm256_xor_ps(ld1, neg); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C*D) == conj(C)*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        ld3 = _mm256_xor_ps(ld3, neg); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm256_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute conj(A)/B
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
    ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
    ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
    ld1 = _mm256_xor_ps(ld1, neg); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    ld1 = _mm256_div_ps(ld1, ld2);// (conj(A)*conj(B))/|B|^2 == conj(A)/B
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX conj(complex) / complex
  __attribute__((__target__("default"))) // default conj(complex) / complex
  inline void divc1(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    // conj(in1[i]) without the complex<float>(re, im) constructor, see the
    // icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> a = zero;
    for(int i=0; i<len; ++i)
    {
      a.real(in1[i].real());
      a.imag(-in1[i].imag()); // conj(in1[i])
      out[i] = a/in2[i];
    }
  }
  
  // complex conjugate / real
  #if !defined(DISABLE_AVX512) // AVX512 conj(complex) / real
  __attribute__((__target__("avx512f")))
  inline void divc1(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    const __m512 one = _mm512_set1_ps(1.0f);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      // negate every other element, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_xor_ps(casting should be free)
      ld2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld2),neg));
      ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
      // compute 1/in2 with a single divide, then spread it to each complex
      ld1 = _mm512_div_ps(one, ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      ld2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld2),neg));
      ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
      // compute 1/in2 with a single divide, then spread it to each complex
      ld1 = _mm512_maskz_div_ps(mask16(rem), one, ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_maskz_mul_ps(mk, ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      ld2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld2),neg));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_maskz_div_ps(mk, ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 conj(complex) / real
  #if !defined(DISABLE_AVX2) // AVX2 conj(complex) / real
  __attribute__((__target__("avx2")))
  inline void divc1(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    const __m256 one = _mm256_set1_ps(1.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
      ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
      // compute 1/in2 with a single divide, then spread it to each complex
      ld1 = _mm256_div_ps(one, ld1);
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
        // compute 1/in2 with a single divide, then spread it to each complex
        ld1 = _mm256_div_ps(one, ld1);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_div_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // AVX2 conj(complex) / real
  #if !defined(DISABLE_AVX) // AVX conj(complex) / real
  __attribute__((__target__("avx")))
  inline void divc1(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    const __m256 one = _mm256_set1_ps(1.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
      ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
      // compute 1/in2 with a single divide, then spread it to each complex
      ld1 = _mm256_div_ps(one, ld1);
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
        // compute 1/in2 with a single divide, then spread it to each complex
        ld1 = _mm256_div_ps(one, ld1);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_div_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX conj(complex) / real
  __attribute__((__target__("default"))) // default conj(complex) / real
  inline void divc1(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    // written out rather than conj(in1[i])/in2[i], see the icpx note above
    for(int i=0; i<len; ++i)
    {
      const float re = in1[i].real()/in2[i];
      const float im = -in1[i].imag()/in2[i];
      out[i].real(re);
      out[i].imag(im);
    }
  }
  
  // real / complex conjugate
  #if !defined(DISABLE_AVX512) // AVX512 real / conj(complex)
  __attribute__((__target__("avx512f")))
  inline void divc2(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    // permutation constants to copy each real element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 8)
    {
      const __m512 one = _mm512_set1_ps(1.0f);
      for(; i<len-15; i+=16) // process 16 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_loadu_ps(&in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * in2[i+8]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm512_mul_ps(ld3, rc); // (in1*in2[i+8])/|in2[i+8]|^2
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * in2[i+8]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc); // (in1*in2[i+8])/|in2[i+8]|^2
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load in1 and in2 and compute in1/conj(in2)
    ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i]));
    sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
    ld1 = _mm512_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * in2[i]
    sc2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm512_add_ps(sc2, ld1);  // |in2[i]|^2
    ld2 = _mm512_maskz_div_ps(mk, ld2, sc2);  // (in1*in2[i])/|in2[i]|^2
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    return;
  }
  #endif // end AVX512 real / conj(complex)
  #if !defined(DISABLE_AVX2) // AVX2 real / conj(complex)
  __attribute__((__target__("avx2,fma")))
  inline void divc2(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 8 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * in2[i+4]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*in2[i+4])/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        // note msk2 accounts for 2 reals per element for the complex buffer
        const __m256i msk1 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem]));
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_maskload_ps(&in1[i], msk1);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * in2[i+4]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*in2[i+4])/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // note msk2 accounts for 2 reals per element for the complex buffer
    const __m256i msk1 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem]));
    const __m256i msk2 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load in1 and in2 and compute in1/conj(in2)
    ld1 = _mm256_maskload_ps(&in1[i], msk1);
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
        msk2);
    sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
    ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
    sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
    ld2 = _mm256_div_ps(ld2, sc2);  // (in1*in2[i])/|in2[i]|^2
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
    return;
  }
  #endif // AVX2 real / conj(complex)
  #if !defined(DISABLE_AVX) // AVX real / conj(complex)
  __attribute__((__target__("avx")))
  inline void divc2(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 8 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * in2[i+4]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*in2[i+4])/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        // note msk2 accounts for 2 reals per element for the complex buffer
        const __m256i msk1 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem]));
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_maskload_ps(&in1[i], msk1);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * in2[i+4]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*in2[i+4])/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // note msk2 accounts for 2 reals per element for the complex buffer
    const __m256i msk1 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem]));
    const __m256i msk2 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load in1 and in2 and compute in1/conj(in2)
    ld1 = _mm256_maskload_ps(&in1[i], msk1);
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
        msk2);
    sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
    ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
    sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
    ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
    sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
    ld2 = _mm256_div_ps(ld2, sc2);  // (in1*in2[i])/|in2[i]|^2
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
    return;
  }
  #endif // end AVX real / conj(complex)
  __attribute__((__target__("default"))) // default real / conj(complex)
  inline void divc2(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    // complex(in1[i], 0) and conj(in2[i]) without the complex<float>(re, im)
    // constructor, see the icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> a = zero; // imaginary part stays 0
    complex<float> b = zero;
    for(int i=0; i<len; ++i)
    {
      a.real(in1[i]);
      b.real(in2[i].real());
      b.imag(-in2[i].imag()); // conj(in2[i])
      out[i] = a/b;
    }
  }

  // complex / complex conjugate
  #if !defined(DISABLE_AVX512) // AVX512 complex / conj(complex)
  __attribute__((__target__("avx512f")))
  inline void divc2(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, ld4, sh, re, im, rc;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 8)
    {
      const __m512 one = _mm512_set1_ps(1.0f);
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld4 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm512_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in1[i+8]));
        ld4 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i]+8));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load A and B and compute A/conj(B)
    ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i])); // A
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i])); // B
    sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
    im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
    re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
    ld2 = _mm512_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br7)^2,(Bi7)^2]
    im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
    sh = _mm512_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm512_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
    ld2 = _mm512_add_ps(ld2, sh); // |B|^2
    ld1 = _mm512_maskz_div_ps(mk, ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    return;
  }
  #endif // AVX512 complex / conj(complex)
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void divc2(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm256_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm256_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute A/conj(B)
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm256_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    ld1 = _mm256_div_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX2 complex / conj(complex)
  #if !defined(DISABLE_AVX) // AVX complex / conj(complex)
  __attribute__((__target__("avx")))
  inline void divc2(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // divide is the heaviest operation, and operating on 2 registers at once
    // can be done with a single divide instead of 2
    if(len > 4)
    {
      const __m256 one = _mm256_set1_ps(1.0f);
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm256_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm256_div_ps(one,rc);// [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm256_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute A/conj(B)
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
    ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
    ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    ld1 = _mm256_div_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX complex / conj(complex)
  __attribute__((__target__("default"))) // default complex / conj(complex)
  inline void divc2(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    // conj(in2[i]) without the complex<float>(re, im) constructor, see the
    // icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> b = zero;
    for(int i=0; i<len; ++i)
    {
      b.real(in2[i].real());
      b.imag(-in2[i].imag()); // conj(in2[i])
      out[i] = in1[i]/b;
    }
  }
  
  // complex conjugate / complex using reciprocal approximation
  #if !defined(DISABLE_AVX512) // AVX512 conj(complex) * recip(complex)
  __attribute__((__target__("avx512f")))
  inline void divc1nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, ld4, sh, re, im, rc;
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m512 two = _mm512_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld4 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B), note the casting is done
        // to avoid requiring the AVX512DQ instruction _mm512_xor_ps (casting
        // should be free)
        ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm512_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in1[i+8]));
        ld4 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i]+8));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B), note the casting is done
        // to avoid requiring the AVX512DQ instruction _mm512_xor_ps (casting
        // should be free)
        ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load A and B and compute conj(A)/B
    ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i])); // A
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i])); // B
    // conj(A), note the casting is done to avoid requiring the AVX512DQ
    // instruction _mm512_xor_ps (casting should be free)
    ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
    sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
    im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
    re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
    ld2 = _mm512_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br7)^2,(Bi7)^2]
    im = _mm512_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai7*Bi7,Ar7*Bi7]
    sh = _mm512_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
    ld2 = _mm512_add_ps(ld2, sh); // |B|^2
    // compute 1/|B|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    re = _mm512_maskz_rcp14_ps(mk, ld2);    // inv
    ld2 = _mm512_fnmadd_ps(ld2, re, two); // 2-orig*inv
    ld2 = _mm512_mul_ps(re, ld2); // inv*(2-orig*inv) == 1/|B|^2
    // (conj(A)*conj(B))/|B|^2 == conj(A)/B
    ld1 = _mm512_mul_ps(ld1, ld2);
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    return;
  }
  #endif // AVX512 conj(complex) * recip(complex) with NR
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void divc1nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B)
        ld1 = _mm256_xor_ps(ld1, neg); // conj(A)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm256_xor_ps(ld3, neg); // conj(C)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm256_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B)
        ld1 = _mm256_xor_ps(ld1, neg); // conj(A)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm256_xor_ps(ld3, neg); // conj(C)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm256_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute conj(A)/B
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    ld1 = _mm256_xor_ps(ld1, neg); // conj(A)
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai3*Bi3,Ar3*Bi3]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm256_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    // compute 1/|B|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    re = _mm256_rcp_ps(ld2);    // inv
    ld2 = _mm256_fnmadd_ps(ld2, re, two); // 2-orig*inv
    ld2 = _mm256_mul_ps(re, ld2); // inv*(2-orig*inv) == 1/|B|^2
    ld1 = _mm256_mul_ps(ld1, ld2);// (conj(A)*conj(B))/|B|^2 == conj(A)/B
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX2 conj(complex) * recip(complex) with NR
  #if !defined(DISABLE_AVX) // AVX conj(complex) * recip(complex) with NR
  __attribute__((__target__("avx")))
  inline void divc1nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A*B) == conj(A)*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        ld1 = _mm256_xor_ps(ld1, neg); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C*D) == conj(C)*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        ld3 = _mm256_xor_ps(ld3, neg); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm256_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A*B) == conj(A)*conj(B)
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        ld1 = _mm256_xor_ps(ld1, neg); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C*D) == conj(C)*conj(D)
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        ld3 = _mm256_xor_ps(ld3, neg); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm256_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute conj(A)/B
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
    ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
    ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
    ld1 = _mm256_xor_ps(ld1, neg); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    // compute 1/|B|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    re = _mm256_rcp_ps(ld2);    // inv
    ld2 = _mm256_mul_ps(ld2, re); // orig*inv
    ld2 = _mm256_sub_ps(two, ld2);// 2-orig*inv
    ld2 = _mm256_mul_ps(re, ld2); // inv*(2-orig*inv) == 1/|B|^2
    ld1 = _mm256_mul_ps(ld1, ld2);// (conj(A)*conj(B))/|B|^2 == conj(A)/B
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX conj(complex) * recip(complex) with NR
  __attribute__((__target__("default"))) // default conj(complex)*recip(complex)
  inline void divc1nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    // conj(in1[i]) without the complex<float>(re, im) constructor, see the
    // icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> a = zero;
    for(int i=0; i<len; ++i)
    {
      a.real(in1[i].real());
      a.imag(-in1[i].imag()); // conj(in1[i])
      out[i] = a/in2[i];
    }
  }

  // complex conjugate / real using reciprocal approximation and 1 stage
  // Newton-Raphson
  #if !defined(DISABLE_AVX512) // AVX512 conj(complex) * recip(real) with NR
  __attribute__((__target__("avx512f")))
  inline void divc1nr(complex<float> const * const in1,
      float const * const in2, const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m512 two = _mm512_set1_ps(2.0f);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc = _mm512_rcp14_ps(ld1);
      ld1 = _mm512_fnmadd_ps(ld1, sc, two); // 2-orig*inv
      ld1 = _mm512_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
      // negate every other element, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_xor_ps (casting should be free)
      ld2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld2),neg));
      ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc = _mm512_maskz_rcp14_ps(mask16(rem), ld1);
      ld1 = _mm512_fnmadd_ps(ld1, sc, two); // 2-orig*inv
      ld1 = _mm512_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
      ld2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld2),neg));
      ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps (mask16(rem), &in2[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc = _mm512_maskz_rcp14_ps(mask16(rem), ld1);
      ld1 = _mm512_fnmadd_ps(ld1, sc, two); // 2-orig*inv
      ld1 = _mm512_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
      ld2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld2),neg));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 conj(complex) * recip(real) with NR
  #if !defined(DISABLE_AVX2) // AVX2 conj(complex) * recip(real) with NR
  __attribute__((__target__("avx2,fma")))
  inline void divc1nr(complex<float> const * const in1,
      float const * const in2, const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld1 = _mm256_loadu_ps(&in2[i]);
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc = _mm256_rcp_ps(ld1);
      ld1 = _mm256_fnmadd_ps(ld1, sc, two); // 2-orig*inv
      ld1 = _mm256_mul_ps(sc, ld1);   // inv*(2-orig*inv) == 1/in2
      ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
      ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sc = _mm256_rcp_ps(ld1);
        ld1 = _mm256_fnmadd_ps(ld1, sc, two); // 2-orig*inv
        ld1 = _mm256_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sc = _mm256_rcp_ps(ld1);
        ld1 = _mm256_fnmadd_ps(ld1, sc, two); // 2-orig*inv
        ld1 = _mm256_mul_ps(sc, ld1); // inv*(2-orig*inv) == 1/in2
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX2 conj(complex) * recip(real) with NR
  #if !defined(DISABLE_AVX) // AVX conj(complex) * recip(real) with NR
  __attribute__((__target__("avx")))
  inline void divc1nr(complex<float> const * const in1,
      float const * const in2, const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
      sc1 = _mm256_rcp_ps(ld1);    // inv
      ld1 = _mm256_mul_ps(ld1, sc1); // orig*inv
      ld1 = _mm256_sub_ps(two, ld1); // 2-orig*inv
      ld1 = _mm256_mul_ps(sc1, ld1); // inv*(2-orig*inv) == 1/in2
      ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
      ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
         // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sc1 = _mm256_rcp_ps(ld1);    // inv
        ld1 = _mm256_mul_ps(ld1, sc1); // orig*inv
        ld1 = _mm256_sub_ps(two, ld1); // 2-orig*inv
        ld1 = _mm256_mul_ps(sc1, ld1); // inv*(2-orig*inv) == 1/in2
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        // compute 1/in2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sc1 = _mm256_rcp_ps(ld1);    // inv
        ld1 = _mm256_mul_ps(ld1, sc1); // orig*inv
        ld1 = _mm256_sub_ps(two, ld1); // 2-orig*inv
        ld1 = _mm256_mul_ps(sc1, ld1); // inv*(2-orig*inv) == 1/in2
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX conj(complex) * recip(real) with NR
  __attribute__((__target__("default"))) // default conj(complex)*recip(real)
  inline void divc1nr(complex<float> const * const in1,
      float const * const in2, const int len, complex<float> * const out)
  {
    // written out rather than conj(in1[i])/in2[i], see the icpx note above
    for(int i=0; i<len; ++i)
    {
      const float re = in1[i].real()/in2[i];
      const float im = -in1[i].imag()/in2[i];
      out[i].real(re);
      out[i].imag(im);
    }
  }
  
  // real / complex conjugate using reciprocal approximation and 1 stage of
  // Newton-Raphson refinement
  #if !defined(DISABLE_AVX512) // AVX512 real * recip(conj(complex)) with NR
  __attribute__((__target__("avx512f")))
  inline void divc2nr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    // permutation constants to copy each real element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m512 two = _mm512_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 16 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_loadu_ps(&in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * in2[i+8]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm512_mul_ps(ld3, rc); // (in1*in2[i+8])/|in2[i+8]|^2
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * in2[i+8]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc); // (in1*in2[i+8])/|in2[i+8]|^2
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load in1 and in2 and compute in1/conj(in2)
    ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i]));
    sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
    ld1 = _mm512_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * in2[i]
    sc2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm512_add_ps(sc2, ld1);  // |in2[i]|^2
    // compute 1/|in2[i]|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    sh = _mm512_maskz_rcp14_ps(mk, sc2);    // inv
    sc2 = _mm512_fnmadd_ps(sc2, sh, two); // 2-orig*inv
    sc2 = _mm512_mul_ps(sh, sc2); // inv*(2-orig*inv) == 1/|in2[i]|^2
    ld2 = _mm512_mul_ps(ld2, sc2);  // (in1*in2[i])/|in2[i]|^2
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    return;
  }
  #endif // end AVX512 real * recip(conj(complex)) with NR
  #if !defined(DISABLE_AVX2) // AVX2 real * recip(conj(complex)) with NR
  __attribute__((__target__("avx2,fma")))
  inline void divc2nr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 8 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * in2[i+4]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*in2[i+4])/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        // note msk2 accounts for 2 reals per element for the complex buffer
        const __m256i msk1 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem]));
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_maskload_ps(&in1[i], msk1);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * in2[i+4]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*in2[i+4])/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // note msk2 accounts for 2 reals per element for the complex buffer
    const __m256i msk1 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem]));
    const __m256i msk2 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load in1 and in2 and compute in1/conj(in2)
    ld1 = _mm256_maskload_ps(&in1[i], msk1);
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
        msk2);
    sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
    ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
    sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
    // compute 1/|in2[i]|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    sh = _mm256_rcp_ps(sc2);    // inv
    sc2 = _mm256_fnmadd_ps(sc2, sh, two); // 2-orig*inv
    sc2 = _mm256_mul_ps(sh, sc2); // inv*(2-orig*inv) == 1/|in2[i]|^2
    ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*in2[i])/|in2[i]|^2
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
    return;
  }
  #endif // AVX2 real * recip(conj(complex)) with NR
  #if !defined(DISABLE_AVX) // AVX real * recip(conj(complex)) with NR
  __attribute__((__target__("avx")))
  inline void divc2nr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 8 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * in2[i+4]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*in2[i+4])/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        // note msk2 accounts for 2 reals per element for the complex buffer
        const __m256i msk1 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem]));
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load in1, B = in2[i], and D = in2[i+4]
        ld1 = _mm256_maskload_ps(&in1[i], msk1);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * in2[i+4]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm256_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, rc); // (in1*in2[i+4])/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // note msk2 accounts for 2 reals per element for the complex buffer
    const __m256i msk1 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem]));
    const __m256i msk2 = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load in1 and in2 and compute in1/conj(in2)
    ld1 = _mm256_maskload_ps(&in1[i], msk1);
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
        msk2);
    sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
    ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
    sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
    ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
    sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
    // compute 1/|in2[i]|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    sh = _mm256_rcp_ps(sc2);    // inv
    sc2 = _mm256_mul_ps(sc2, sh); // orig*inv
    sc2 = _mm256_sub_ps(two, sc2);// 2-orig*inv
    sc2 = _mm256_mul_ps(sh, sc2); // inv*(2-orig*inv) == 1/|in2[i]|^2
    ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*in2[i])/|in2[i]|^2
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
    return;
  }
  #endif // end AVX real * recip(conj(complex)) with NR
  __attribute__((__target__("default"))) // default real * recip(conj(complex))
  inline void divc2nr(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    // complex(in1[i], 0) and conj(in2[i]) without the complex<float>(re, im)
    // constructor, see the icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> a = zero; // imaginary part stays 0
    complex<float> b = zero;
    for(int i=0; i<len; ++i)
    {
      a.real(in1[i]);
      b.real(in2[i].real());
      b.imag(-in2[i].imag()); // conj(in2[i])
      out[i] = a/b;
    }
  }

  // complex / complex conjugate using reciprocal approximation with 1 stage of
  // Newton-Raphson refinement
  #if !defined(DISABLE_AVX512) // AVX512 complex*recip(conj(complex)) with NR
  __attribute__((__target__("avx512f")))
  inline void divc2nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, ld4, sh, re, im, rc;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m512 two = _mm512_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld4 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm512_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in1[i+8]));
        ld4 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i]+8));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm512_rcp14_ps(rc);   // inv
        rc = _mm512_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm512_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load A and B and compute A/conj(B)
    ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i])); // A
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i])); // B
    sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
    im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
    re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
    ld2 = _mm512_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br7)^2,(Bi7)^2]
    im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
    sh = _mm512_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm512_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
    ld2 = _mm512_add_ps(ld2, sh); // |B|^2
    // compute 1/|B|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    re = _mm512_maskz_rcp14_ps(mk, ld2);    // inv
    ld2 = _mm512_fnmadd_ps(ld2, re, two); // 2-orig*inv
    ld2 = _mm512_mul_ps(re, ld2); // inv*(2-orig*inv) == 1/|B|^2
    ld1 = _mm512_mul_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    return;
  }
  #endif // AVX512 complex*recip(conj(complex)) with NR
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void divc2nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm256_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm256_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_fnmadd_ps(rc, sh, two); // 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm256_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute A/conj(B)
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm256_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    // compute 1/|B|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    re = _mm256_rcp_ps(ld2);    // inv
    ld2 = _mm256_fnmadd_ps(ld2, re, two); // 2-orig*inv
    ld2 = _mm256_mul_ps(re, ld2); // inv*(2-orig*inv) == 1/|B|^2
    ld1 = _mm256_mul_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX2 complex*recip(conj(complex)) with NR
  #if !defined(DISABLE_AVX) // AVX complex*recip(conj(complex)) with NR
  __attribute__((__target__("avx")))
  inline void divc2nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, ld3, ld4, sh, re, im, rc;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld4 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm256_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
      }
      if(rem > 4) // one full and one partial register
      {
        const __m256i msk = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk);
        ld4 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk);
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm256_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        re = _mm256_mul_ps(re, re);// [Br0^2,Br1^2,Dr0^2,Dr1^2,...]
        im = _mm256_mul_ps(im, im);// [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm256_add_ps(re, im);// [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        // 1/|B|^2 and 1/|D|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
        sh = _mm256_rcp_ps(rc);   // inv
        rc = _mm256_mul_ps(rc, sh); // orig*inv
        rc = _mm256_sub_ps(two, rc);// 2-orig*inv
        rc = _mm256_mul_ps(sh, rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
        im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
        re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
        im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
        ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
        ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm256_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci3,Cr3]
        im = _mm256_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di3,Di3]
        re = _mm256_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr3,Dr3]
        im = _mm256_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci3*Di3,Cr3*Di3]
        ld3 = _mm256_mul_ps(ld3, re);   // [Cr0*Dr0,Ci0*Dr0,...,Cr3*Dr3,Ci3*Dr3]
        ld3 = _mm256_addsub_ps(ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm256_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm256_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm256_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm256_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __m256i msk = _mm256_load_si256(
        reinterpret_cast<__m256i const *>(masks[rem<<1]));
    // load A and B and compute A/conj(B)
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in1[i]), msk); // A
    ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
        &in2[i]), msk); // B
    sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
    im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
    re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
    ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
    im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
    ld1 = _mm256_mul_ps(ld1, re);   // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
    ld1 = _mm256_addsub_ps(ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
    sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld2 = _mm256_add_ps(ld2, sh); // |B|^2
    // compute 1/|B|^2 with 1 stage Newton-Raphson (inv*(2-orig*inv))
    re = _mm256_rcp_ps(ld2);    // inv
    ld2 = _mm256_mul_ps(ld2, re); // orig*inv
    ld2 = _mm256_sub_ps(two, ld2);// 2-orig*inv
    ld2 = _mm256_mul_ps(re, ld2); // inv*(2-orig*inv) == 1/|B|^2
    ld1 = _mm256_mul_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    return;
  }
  #endif // end AVX complex*recip(conj(complex)) with NR
  __attribute__((__target__("default"))) // default complex*recip(conj(complex))
  inline void divc2nr(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    // conj(in2[i]) without the complex<float>(re, im) constructor, see the
    // icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> b = zero;
    for(int i=0; i<len; ++i)
    {
      b.real(in2[i].real());
      b.imag(-in2[i].imag()); // conj(in2[i])
      out[i] = in1[i]/b;
    }
  }

  // complex conjugate / complex using reciprocal approximation
  #if !defined(DISABLE_AVX512) // AVX512 conj(complex) * recip(complex)
  __attribute__((__target__("avx512f")))
  inline void divc1r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, ld4, sh, re, im, rc;
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld4 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B), note the casting is done
        // to avoid requiring the AVX512DQ instruction _mm512_xor_ps (casting
        // should be free)
        ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        ld3 = _mm512_mul_ps(ld3, rc); // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load A, B, C, and D, and compute conj(A)/B and conj(C)/D
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in1[i+8]));
        ld4 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i]+8));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator conj(A)*conj(B), note the casting is done
        // to avoid requiring the AVX512DQ instruction _mm512_xor_ps (casting
        // should be free)
        ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
        // compute the complex numerator conj(C)*conj(D)
        ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [-Ci0,Cr0,-Ci1,Cr1,...,-Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [-Ci0*Di0,Cr0*Di0,...,-Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmsubadd_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,-Dr0*Ci0-Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (conj(A)*conj(B))/|B|^2 == conj(A)/B
        // (conj(C)*conj(D))/|D|^2 == conj(C)/D
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load A and B and compute conj(A)/B
    ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i])); // A
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i])); // B
    // conj(A), note the casting is done to avoid requiring the AVX512DQ
    // instruction _mm512_xor_ps (casting should be free)
    ld1 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld1),neg));
    sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
    im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
    re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
    ld2 = _mm512_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br7)^2,(Bi7)^2]
    im = _mm512_mul_ps(sh, im);   // [-Ai0*Bi0,Ar0*Bi0,...,-Ai7*Bi7,Ar7*Bi7]
    sh = _mm512_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm512_fmsubadd_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
    ld2 = _mm512_add_ps(ld2, sh); // |B|^2
    ld2 = _mm512_maskz_rcp14_ps(mk, ld2);  // 1/|B|^2
    // (conj(A)*conj(B))/|B|^2 == conj(A)/B
    ld1 = _mm512_mul_ps(ld1, ld2);
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    return;
  }
  #endif // AVX512 conj(complex) * recip(complex)
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void divc1r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld1 = _mm256_xor_ps(ld1, neg);// conj(A)
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1);// [-Ai0,Ar0,-Ai1,Ar1,...,-Ai7,Ar7]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_fmsubadd_ps(re, ld1, im);// [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (conj(A)*conj(B))/|B|^2 == conj(A)/B
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      ld1 = _mm256_xor_ps(ld1, neg);// conj(A)
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1);// [-Ai0,Ar0,-Ai1,Ar1,...,-Ai3,Ar3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_fmsubadd_ps(re, ld1, im);// [Br0*Ar0-Bi0Ai0,-Br0*Ai0-Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (conj(A)*conj(B))/|B|^2 == conj(A)/B
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 conj(complex) * recip(complex)
  #if !defined(DISABLE_AVX) // AVX conj(complex) * recip(complex)
  __attribute__((__target__("avx")))
  inline void divc1r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_addsub_ps(ld1, im);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      ld1 = _mm256_xor_ps(ld1, neg);// conj(A*B) == conj(A)*conj(B)
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (conj(A)*conj(B))/|B|^2 == conj(A)/B
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_addsub_ps(ld1, im);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      ld1 = _mm256_xor_ps(ld1, neg);// conj(A*B) == conj(A)*conj(B)
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (conj(A)*conj(B))/|B|^2 == conj(A)/B
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX conj(complex) * recip(complex)
  __attribute__((__target__("default"))) // default conj(complex)*recip(complex)
  inline void divc1r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    // conj(in1[i]) without the complex<float>(re, im) constructor, see the
    // icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> a = zero;
    for(int i=0; i<len; ++i)
    {
      a.real(in1[i].real());
      a.imag(-in1[i].imag()); // conj(in1[i])
      out[i] = a/in2[i];
    }
  }

  // complex conjugate / real using reciprocal approximation
  #if !defined(DISABLE_AVX512) // AVX512 conj(complex) * recip(real)
  __attribute__((__target__("avx512f")))
  inline void divc1r(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // constant to negate every other value for complex conjugation
    const __m512i neg = _mm512_setr_epi32(0, 0x80000000, 0, 0x80000000, 0,
        0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000,
        0, 0x80000000);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      // negate every other element, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_xor_ps (casting should be free)
      ld2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld2),neg));
      ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
      ld1 = _mm512_rcp14_ps(ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      ld2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld2),neg));
      ld3 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld3),neg));
      ld1 = _mm512_maskz_rcp14_ps(mask16(rem), ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      ld2 = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(ld2),neg));
      ld1 = _mm512_maskz_rcp14_ps(mask16(rem), ld1);
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 conj(complex) * recip(real)
  #if !defined(DISABLE_AVX2) // AVX2 conj(complex) * recip(real)
  __attribute__((__target__("avx2")))
  inline void divc1r(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
      ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
      ld1 = _mm256_rcp_ps(ld1);
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
        ld1 = _mm256_rcp_ps(ld1);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        ld1 = _mm256_rcp_ps(ld1);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // AVX2 conj(complex) * recip(real)
  #if !defined(DISABLE_AVX) // AVX conj(complex) * recip(real)
  __attribute__((__target__("avx")))
  inline void divc1r(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld1 = _mm256_loadu_ps(&in2[i]);
      ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
      ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
      ld1 = _mm256_rcp_ps(ld1);
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        ld3 = _mm256_xor_ps(ld3, neg);  // negate every other element
        ld1 = _mm256_rcp_ps(ld1);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        ld2 = _mm256_xor_ps(ld2, neg);  // negate every other element
        ld1 = _mm256_rcp_ps(ld1);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX conj(complex) * recip(real)
  __attribute__((__target__("default"))) // default conj(complex) * recip(real)
  inline void divc1r(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    // written out rather than conj(in1[i])/in2[i], see the icpx note above
    for(int i=0; i<len; ++i)
    {
      const float re = in1[i].real()/in2[i];
      const float im = -in1[i].imag()/in2[i];
      out[i].real(re);
      out[i].imag(im);
    }
  }
  

  // real / complex conjugate using reciprocal approximation
  #if !defined(DISABLE_AVX512) // AVX512 real * recip(conj(complex))
  __attribute__((__target__("avx512f")))
  inline void divc2r(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, sc1, sc2, sh, re, im, rc;
    // permutation constants to copy each real element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 16 real elements per register
      {
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_loadu_ps(&in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * in2[i+8]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm512_mul_ps(ld3, rc); // (in1*in2[i+8])/|in2[i+8]|^2
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load in1, B = in2[i], and D = in2[i+8]
        ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld3,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld3,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute the numerators in1*B and in1*D
        sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
        sc2 = _mm512_permutexvar_ps(p2, ld1); // repeat each real
        ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * in2[i]
        ld3 = _mm512_mul_ps(ld3, sc2);  // in1 * in2[i+8]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld2 = _mm512_mul_ps(ld2, sh); // (in1*in2[i])/|in2[i]|^2
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc); // (in1*in2[i+8])/|in2[i+8]|^2
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load in1 and in2 and compute in1/conj(in2)
    ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i]));
    sc1 = _mm512_permutexvar_ps(p1, ld1); // repeat each real
    ld1 = _mm512_mul_ps(ld2, ld2);  // in2[i]^2
    ld2 = _mm512_mul_ps(ld2, sc1);  // in1 * in2[i]
    sc2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
    sc2 = _mm512_add_ps(sc2, ld1);  // |in2[i]|^2
    sc2 = _mm512_maskz_rcp14_ps(mk, sc2);  // 1/|in2[i]|^2
    ld2 = _mm512_mul_ps(ld2, sc2);  // (in1*in2[i])/|in2[i]|^2
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    return;
  }
  #endif // end AVX512 real * recip(conj(complex))
  #if !defined(DISABLE_AVX2) // AVX2 real * recip(conj(complex))
  __attribute__((__target__("avx2")))
  inline void divc2r(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2, sh;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
      sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
      ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
      ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * in2[i]
      sc1 = _mm256_mul_ps(ld3, ld3);  // in2[i+4]^2
      ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * in2[i+4]
      sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
      sh = _mm256_shuffle_ps(sc1, sc1, 0xb1);  // in2[i+4]^2 swap(re,im)
      sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
      sh = _mm256_add_ps(sc1, sh);    // |in2[i+4]|^2
      sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
      sh = _mm256_rcp_ps(sh);       // 1/|in2[i+4]|^2
      ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*in2[i])/|in2[i]|^2
      ld3 = _mm256_mul_ps(ld3, sh);   // (in1*in2[i+4])/|in2[i+4]|^2
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);  // repeat first 4 reals
        sc2 = _mm256_permutevar8x32_ps(ld1, p2);  // repeat last 4 reals
        ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        sc1 = _mm256_mul_ps(ld3, ld3);  // in2[i+4]^2
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
        sh = _mm256_shuffle_ps(sc1, sc1, 0xb1);  // in2[i+4]^2 swap(re,im)
        sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
        sh = _mm256_add_ps(sc1, sh);    // |in2[i+4]|^2
        sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
        sh = _mm256_rcp_ps(sh);       // 1/|in2[i+4]|^2
        ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, sh);   // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc1 = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
        sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
        sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
        ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // AVX2 real * recip(conj(complex))
  #if !defined(DISABLE_AVX) // AVX real * recip(conj(complex))
  __attribute__((__target__("avx")))
  inline void divc2r(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc1, sc2, sh;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
      ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
      sc1 = _mm256_mul_ps(ld3, ld3);  // in2[i+4]^2
      ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
      sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
      sh = _mm256_shuffle_ps(sc1, sc1, 0xb1);  // in2[i+4]^2 swap(re,im)
      sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
      sh = _mm256_add_ps(sc1, sh);    // |in2[i+4]|^2
      sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
      sh = _mm256_rcp_ps(sh);       // 1/|in2[i+4]|^2
      ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
      ld3 = _mm256_mul_ps(ld3, sh);   // (in1*conj(in2[i+4]))/|in2[i+4]|^2
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        sc1 = _mm256_mul_ps(ld3, ld3);  // in2[i+4]^2
        ld3 = _mm256_mul_ps(ld3, sc2);  // in1 * conj(in2[i+4])
        sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
        sh = _mm256_shuffle_ps(sc1, sc1, 0xb1);  // in2[i+4]^2 swap(re,im)
        sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
        sh = _mm256_add_ps(sc1, sh);    // |in2[i+4]|^2
        sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
        sh = _mm256_rcp_ps(sh);       // 1/|in2[i+4]|^2
        ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
        ld3 = _mm256_mul_ps(ld3, sh);   // (in1*conj(in2[i+4]))/|in2[i+4]|^2
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld1 = _mm256_mul_ps(ld2, ld2);  // in2[i]^2
        ld2 = _mm256_mul_ps(ld2, sc1);  // in1 * conj(in2[i])
        sc2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // in2[i]^2 swap(re,im)
        sc2 = _mm256_add_ps(sc2, ld1);  // |in2[i]|^2
        sc2 = _mm256_rcp_ps(sc2);     // 1/|in2[i]|^2
        ld2 = _mm256_mul_ps(ld2, sc2);  // (in1*conj(in2[i]))/|in2[i]|^2
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX real * recip(conj(complex))
  __attribute__((__target__("default"))) // default real * recip(conj(complex))
  inline void divc2r(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    // complex(in1[i], 0) and conj(in2[i]) without the complex<float>(re, im)
    // constructor, see the icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> a = zero; // imaginary part stays 0
    complex<float> b = zero;
    for(int i=0; i<len; ++i)
    {
      a.real(in1[i]);
      b.real(in2[i].real());
      b.imag(-in2[i].imag()); // conj(in2[i])
      out[i] = a/b;
    }
  }

  // complex / complex conjugate using reciprocal approximation
  #if !defined(DISABLE_AVX512) // AVX512 complex*recip(conj(complex))
  __attribute__((__target__("avx512f")))
  inline void divc2r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, ld3, ld4, sh, re, im, rc;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once as the
    // 1/|B|^2 values for both registers fit in a single register, so only 1
    // reciprocal approximation is needed instead of 2
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld4 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm512_mul_ps(ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
      }
      if(rem > 8) // one full and one partial register
      {
        const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
        // load A, B, C, and D, and compute A/conj(B) and C/conj(D)
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in1[i+8]));
        ld4 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
            &in2[i]+8));
        // compute 1/|B|^2 and 1/|D|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld2,ld4,0x88); // [Br0,Br1,Dr0,Dr1,...]
        im = _mm512_shuffle_ps(ld2,ld4,0xdd); // [Bi0,Bi1,Di0,Di1,...]
        im = _mm512_mul_ps(im,im); // [Bi0^2,Bi1^2,Di0^2,Di1^2,...]
        rc = _mm512_fmadd_ps(re,re,im); // [|B0|^2,|B1|^2,|D0|^2,|D1|^2,...]
        rc = _mm512_rcp14_ps(rc); // [1/|B0|^2,1/|B1|^2,1/|D0|^2,1/|D1|^2,...]
        // compute complex numerator A*B
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
        im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
        re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
        im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
        ld1 = _mm512_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
        // compute the complex numerator C*D
        sh = _mm512_shuffle_ps(ld3, ld3, 0xb1); // [Ci0,Cr0,Ci1,Cr1,...,Ci7,Cr7]
        im = _mm512_movehdup_ps(ld4); // [Di0,Di0,Di1,Di1,...,Di7,Di7]
        re = _mm512_moveldup_ps(ld4); // [Dr0,Dr0,Dr1,Dr1,...,Dr7,Dr7]
        im = _mm512_mul_ps(sh, im);   // [Ci0*Di0,Cr0*Di0,...,Ci7*Di7,Cr7*Di7]
        ld3 = _mm512_fmaddsub_ps(re, ld3, im); // [Dr0*Cr0-Di0Ci0,Dr0*Ci0+Cr0*Di0]
        // spread reciprocals into 2 registers and scale the numerators
        sh = _mm512_unpacklo_ps(rc, rc); // [1/|B0|^2,1/|B0|^2,1/|B1|^2,...]
        rc = _mm512_unpackhi_ps(rc, rc); // [1/|D0|^2,1/|D0|^2,1/|D1|^2,...]
        ld1 = _mm512_mul_ps(ld1, sh); // (A*B)/|B|^2 == A/conj(B)
        ld3 = _mm512_maskz_mul_ps(mk, ld3, rc); // (C*D)/|D|^2 == C/conj(D)
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
            ld3);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    const __mmask16 mk = mask16(rem<<1);
    // load A and B and compute A/conj(B)
    ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in1[i])); // A
    ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
        &in2[i])); // B
    sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
    im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
    re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
    ld2 = _mm512_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br7)^2,(Bi7)^2]
    im = _mm512_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
    sh = _mm512_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
    ld1 = _mm512_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
    ld2 = _mm512_add_ps(ld2, sh); // |B|^2
    ld2 = _mm512_maskz_rcp14_ps(mk, ld2);  // 1/|B|^2
    ld1 = _mm512_mul_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    return;
  }
  #endif // AVX512 complex*recip(conj(complex))
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void divc2r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);   // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_fmaddsub_ps(re, ld1, im); // [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex*recip(conj(complex))
  #if !defined(DISABLE_AVX) // AVX complex*recip(conj(complex))
  __attribute__((__target__("avx")))
  inline void divc2r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_addsub_ps(ld1, im);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
          msk);
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(ld2, ld2);// [(Br0)^2,(Bi0)^2,...,(Br3)^2,(Bi3)^2]
      im = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      sh = _mm256_shuffle_ps(ld2, ld2, 0xb1);// [(Bi0)^2,(Br0)^2,...]
      ld1 = _mm256_addsub_ps(ld1, im);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      ld2 = _mm256_add_ps(ld2, sh); // |B|^2
      ld2 = _mm256_rcp_ps(ld2);   // 1/|B|^2
      ld1 = _mm256_mul_ps(ld1, ld2);// (A*B)/|B|^2 == A/conj(B)
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex*recip(conj(complex))
  __attribute__((__target__("default"))) // default complex*recip(conj(complex))
  inline void divc2r(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    // conj(in2[i]) without the complex<float>(re, im) constructor, see the
    // icpx note above
    static constexpr complex<float> zero(0.0f, 0.0f);
    complex<float> b = zero;
    for(int i=0; i<len; ++i)
    {
      b.real(in2[i].real());
      b.imag(-in2[i].imag()); // conj(in2[i])
      out[i] = in1[i]/b;
    }
  }

  //=====================================//
  // Dot-Product Function Implementation //
  //=====================================//
  // real*real dot-product
  #if !defined(DISABLE_AVX512) // AVX512 real*real dotp
  __attribute__((__target__("avx512f")))
  inline float dotp(float const * const in1, float const * const in2,
      const int len)
  {
    if(len<=0) return 0;
    __m512 ld1, ld2;
    __m512 sum = _mm512_set1_ps(0.0f);
    int i = 0;
    // note: with a single accumulator every FMA waits on the previous one,
    // instead use 4 independent accumulators to keep the pipeline full
    if(len > 63)
    {
      __m512 sum1 = sum;
      __m512 sum2 = sum;
      __m512 sum3 = sum;
      for(; i<len-63; i+=64) // process 4 registers (64 real elements)
      {
        ld1 = _mm512_loadu_ps(&in1[i]);
        ld2 = _mm512_loadu_ps(&in2[i]);
        sum = _mm512_fmadd_ps(ld1, ld2, sum);
        ld1 = _mm512_loadu_ps(&in1[i+16]);
        ld2 = _mm512_loadu_ps(&in2[i+16]);
        sum1 = _mm512_fmadd_ps(ld1, ld2, sum1);
        ld1 = _mm512_loadu_ps(&in1[i+32]);
        ld2 = _mm512_loadu_ps(&in2[i+32]);
        sum2 = _mm512_fmadd_ps(ld1, ld2, sum2);
        ld1 = _mm512_loadu_ps(&in1[i+48]);
        ld2 = _mm512_loadu_ps(&in2[i+48]);
        sum3 = _mm512_fmadd_ps(ld1, ld2, sum3);
      }
      sum = _mm512_add_ps(sum, sum1);
      sum2 = _mm512_add_ps(sum2, sum3);
      sum = _mm512_add_ps(sum, sum2);
    }
    // 0-3 full registers remain and up to 1 partial register, compute these
    // individually as using a loop will cause intel compiler to save/restore
    // additional registers
    // bit 5 indicates there are at least 32 elements that weren't taken care
    // of by the loop above
    if(len&32) // process 2 registers (32 real elements)
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      sum = _mm512_fmadd_ps(ld1, ld2, sum);
      ld1 = _mm512_loadu_ps(&in1[i+16]);
      ld2 = _mm512_loadu_ps(&in2[i+16]);
      sum = _mm512_fmadd_ps(ld1, ld2, sum);
      i += 32;
    }
    // bit 4 indicates there were an additional 16 elements not taken care of
    // by the loop above
    if(len&16) // 1 register (16 real elements)
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      sum = _mm512_fmadd_ps(ld1, ld2, sum);
      i += 16;
    }
    // handle remaining elements in partial register (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      sum = _mm512_fmadd_ps(ld1, ld2, sum);
    }
    return _mm512_reduce_add_ps(sum);
  }
  #endif // end AVX512 real*real dotp
  #if !defined(DISABLE_AVX2) // AVX2 real*real dotp
  __attribute__((__target__("avx2,fma")))
  inline float dotp(float const * const in1, float const * const in2,
      const int len)
  {
    if(len<=0) return 0;
    __m256 ld1, ld2;
    __m256 sum = _mm256_set1_ps(0.0f);
    __m128 hs;  // horizontal sum to finish dot-product
    int i = 0;
    // note: with a single accumulator every FMA waits on the previous one,
    // instead use 4 independent accumulators to keep the pipeline full
    if(len > 31)
    {
      __m256 sum1 = sum;
      __m256 sum2 = sum;
      __m256 sum3 = sum;
      for(; i<len-31; i+=32) // process 4 registers (32 real elements)
      {
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(&in2[i]);
        sum = _mm256_fmadd_ps(ld1, ld2, sum);
        ld1 = _mm256_loadu_ps(&in1[i+8]);
        ld2 = _mm256_loadu_ps(&in2[i+8]);
        sum1 = _mm256_fmadd_ps(ld1, ld2, sum1);
        ld1 = _mm256_loadu_ps(&in1[i+16]);
        ld2 = _mm256_loadu_ps(&in2[i+16]);
        sum2 = _mm256_fmadd_ps(ld1, ld2, sum2);
        ld1 = _mm256_loadu_ps(&in1[i+24]);
        ld2 = _mm256_loadu_ps(&in2[i+24]);
        sum3 = _mm256_fmadd_ps(ld1, ld2, sum3);
      }
      sum = _mm256_add_ps(sum, sum1);
      sum2 = _mm256_add_ps(sum2, sum3);
      sum = _mm256_add_ps(sum, sum2);
    }
    // 0-3 full registers remain and up to 1 partial register, compute these
    // individually as using a loop will cause intel compiler to save/restore
    // additional registers
    // bit 4 indicates there are at least 16 elements that weren't taken care
    // of by the loop above
    if(len&16) // 2 registers (16 real elements)
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      sum = _mm256_fmadd_ps(ld1, ld2, sum);
      ld1 = _mm256_loadu_ps(&in1[i+8]);
      ld2 = _mm256_loadu_ps(&in2[i+8]);
      sum = _mm256_fmadd_ps(ld1, ld2, sum);
      i += 16;
    }
    // bit 3 indicates there were an additional 8 elements not taken care of
    // by the loop above
    if(len&8) // 1 register (8 real elements)
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      sum = _mm256_fmadd_ps(ld1, ld2, sum);
      i += 8;
    }
    // handle remaining elements in partial register (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      sum = _mm256_fmadd_ps(ld1, ld2, sum);
    }
    // horizontal sum
    hs = _mm256_extractf128_ps(sum, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(sum), hs); // [0+4, 1+5, 2+6, 3+7]
    hs = _mm_add_ps(hs, _mm_movehl_ps(hs, hs)); // [0+4+2+6, 1+5+3+7, ...]
    hs = _mm_add_ss(hs, _mm_movehdup_ps(hs));   // [0+1+2+3+4+5+6+7, ...]
    return _mm_cvtss_f32(hs); // extract the first element in the register
  }
  #endif // end AVX2 real*real dotp
  #if !defined(DISABLE_AVX) // AVX real*real dotp
  __attribute__((__target__("avx")))
  inline float dotp(float const * const in1, float const * const in2,
      const int len)
  {
    if(len<=0) return 0;
    __m256 ld1, ld2;
    __m256 sum = _mm256_set1_ps(0.0f);
    __m128 hs;  // horizontal sum to finish dot-product
    int i = 0;
    // note: with a single accumulator every mul/add waits on the previous one,
    // instead use 4 independent accumulators to keep the pipeline full
    if(len > 31)
    {
      __m256 sum1 = sum;
      __m256 sum2 = sum;
      __m256 sum3 = sum;
      for(; i<len-31; i+=32) // process 4 registers (32 real elements)
      {
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(&in2[i]);
        ld1 = _mm256_mul_ps(ld1, ld2);
        sum = _mm256_add_ps(ld1, sum);
        ld1 = _mm256_loadu_ps(&in1[i+8]);
        ld2 = _mm256_loadu_ps(&in2[i+8]);
        ld1 = _mm256_mul_ps(ld1, ld2);
        sum1 = _mm256_add_ps(ld1, sum1);
        ld1 = _mm256_loadu_ps(&in1[i+16]);
        ld2 = _mm256_loadu_ps(&in2[i+16]);
        ld1 = _mm256_mul_ps(ld1, ld2);
        sum2 = _mm256_add_ps(ld1, sum2);
        ld1 = _mm256_loadu_ps(&in1[i+24]);
        ld2 = _mm256_loadu_ps(&in2[i+24]);
        ld1 = _mm256_mul_ps(ld1, ld2);
        sum3 = _mm256_add_ps(ld1, sum3);
      }
      sum = _mm256_add_ps(sum, sum1);
      sum2 = _mm256_add_ps(sum2, sum3);
      sum = _mm256_add_ps(sum, sum2);
    }
    // 0-3 full registers remain and up to 1 partial register, compute these
    // individually as using a loop will cause intel compiler to save/restore
    // additional registers
    // bit 4 indicates there are at least 16 elements that weren't taken care
    // of by the loop above
    if(len&16) // 2 registers (16 real elements)
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_mul_ps(ld1, ld2);
      sum = _mm256_add_ps(ld1, sum);
      ld1 = _mm256_loadu_ps(&in1[i+8]);
      ld2 = _mm256_loadu_ps(&in2[i+8]);
      ld1 = _mm256_mul_ps(ld1, ld2);
      sum = _mm256_add_ps(ld1, sum);
      i += 16;
    }
    // bit 3 indicates there were an additional 8 elements not taken care of
    // by the loop above
    if(len&8) // 1 register (8 real elements)
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_mul_ps(ld1, ld2);
      sum = _mm256_add_ps(ld1, sum);
      i += 8;
    }
    // handle remaining elements in partial register (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld2);
      sum = _mm256_add_ps(ld1, sum);
    }
    // horizontal sum
    hs = _mm256_extractf128_ps(sum, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(sum), hs); // [0+4, 1+5, 2+6, 3+7]
    hs = _mm_add_ps(hs, _mm_movehl_ps(hs, hs)); // [0+4+2+6, 1+5+3+7, ...]
    hs = _mm_add_ss(hs, _mm_movehdup_ps(hs));   // [0+1+2+3+4+5+6+7, ...]
    return _mm_cvtss_f32(hs); // extract the first element in the register
  }
  #endif // end AVX real*real dotp
  __attribute__((__target__("default"))) // default real*real dotp
  inline float dotp(float const * const in1, float const * const in2,
      const int len)
  {
    float sum = 0.0f;
    for(int i=0; i<len; ++i) sum += (in1[i]*in2[i]);
    return sum;
  }





/*
  #if !defined(DISABLE_AVX512) // AVX512 complex*complex dotp
  __attribute__((__target__("avx512f")))
  
  #endif // end AVX512 complex*complex dotp
  #if !defined(DISABLE_AVX2) // AVX2 complex*complex dotp
  __attribute__((__target__("avx2,fma")))
  
  #endif // end AVX2 complex*complex dotp
  #if !defined(DISABLE_AVX) // AVX complex*complex dotp
  __attribute__((__target__("avx")))

  #endif // end AVX complex*complex dotp
  __attribute__((__target__("default"))) // default complex*complex dotp
  inline complex<float> dotp(complex<float> const * const in1,
      complex<float> const * const in2, const int len)
  {

  }
  



  #if !defined(DISABLE_AVX512) // AVX512 real*complex dotp
  __attribute__((__target__("avx512f")))
  
  #endif // end AVX512 real*complex dotp
  #if !defined(DISABLE_AVX2) // AVX2 real*complex dotp
  __attribute__((__target__("avx2,fma")))
  
  #endif // end AVX2 real*complex dotp
  #if !defined(DISABLE_AVX) // AVX real*complex dotp
  __attribute__((__target__("avx")))

  #endif // end AVX real*complex dotp
  __attribute__((__target__("default"))) // default real*complex dotp
  inline complex<float> dotp(float const * const in1,
      complex<float> const * const in2, const int len)
  {

  }
  


  #if !defined(DISABLE_AVX512) // AVX512 complex*conj(complex) dotp
  __attribute__((__target__("avx512f")))
  
  #endif // end AVX512 complex*conj(complex) dotp
  #if !defined(DISABLE_AVX2) // AVX2 complex*conj(complex) dotp
  __attribute__((__target__("avx2,fma")))
  
  #endif // end AVX2 complex*conj(complex) dotp
  #if !defined(DISABLE_AVX) // AVX complex*conj(complex) dotp
  __attribute__((__target__("avx")))

  #endif // end AVX complex*conj(complex) dotp
  __attribute__((__target__("default"))) // default complex*conj(complex) dotp
  inline complex<float> dotpc(complex<float> const * const in1,
      complex<float> const * const in2, const int len)
  {

  }
  



  #if !defined(DISABLE_AVX512) // AVX512 real*conj(complex) dotp
  __attribute__((__target__("avx512f")))
  
  #endif // end AVX512 real*conj(complex) dotp
  #if !defined(DISABLE_AVX2) // AVX2 real*conj(complex) dotp
  __attribute__((__target__("avx2,fma")))
  
  #endif // end AVX2 real*conj(complex) dotp
  #if !defined(DISABLE_AVX) // AVX real*conj(complex) dotp
  __attribute__((__target__("avx")))

  #endif // end AVX real*conj(complex) dotp
  __attribute__((__target__("default"))) // default real*conj(complex) dotp
  inline complex<float> dotpc(float const * const in1,
      complex<float> const * const in2, const int len)
  {

  }
*/


  //========================================//
  // Raise to Power Function Implementation //
  //========================================//

  //=====================================================//
  // Magnitude/Magnitude Squared Function Implementation //
  //=====================================================//

  //===============================================//
  // Specialized Functions Function Implementation //
  //===============================================//

} // namespace EVM
#endif // EFFICIENT_VECTOR_MATH_H