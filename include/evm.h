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
//  Independently, the polar/cartesian conversions use the Intel Short Vector
//  Math Library (SVML) for their sine/cosine of radians and atan2 when it is
//  available, otherwise their own SIMD versions in evm_trig.h, which are
//  slightly less accurate (up to 1.6-3.3 ulp vs. 0.6-0.9 ulp).  polToCart
//  (phase in cycles) always uses evm_trig.h, which is faster than SVML and
//  much more accurate for large phases.  SVML will automatically be enabled
//  when the Intel compiler is detected but may be manually disabled with:
//    -DDISABLE_SVML (does not affect SIMD)
//
//  Note that currently there are only implementations for AVX, AVX2, and AVX512
//  branches along with a default (non-SIMD) branch, so disabling the SSE
//  branches will not have any practical effect.
//
//  The AVX512 versions of the 8 and 16-bit integer conversions also require
//  AVX512BW and AVX512VL, CPUs without them (Xeon Phi) use the AVX2 versions.
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
//  Comparison:
//    min = point-by-point comparison of two buffers, placing the minimum in
//        the output buffer
//    max = point-by-point comparison of two buffers, placing the maximum in
//        the output buffer
//    minMax = compute the minimum and maximum values of an input buffer or
//        minimum and maximum or each element of two buffers
//    zpos = set the positive values of a buffer to zero
//    zneg = set the negative values of a buffer to zero
//    signedAbsMin = point-by-point minimum of the absolute value of each
//        element in a series of buffers
//    signedAbsMax = point-by-point maximum of the absolute value of each
//        element in a series of buffers
//    threshCmp = comparison of each element of a buffer with a threshold,
//        choosing to output the element if it is above the threshold or the
//        value of a second buffer if it is below the threshold
//    blend = comparison of each value of a buffer to a threshold, outputting
//        the value of either of two other buffers depending on if the value
//        is greater than or less than the threshold
//
//  Multiplication:
//    scale = multiply an input buffer by a constant scale factor
//    mult = point-by-point vector multiply
//    multc = point-by-point complex conjugate multiply, second input conjugated
//    multcs = point-by-point complex conjugate multiply, second input
//        conjugated, with a constant scale factor applied
//    multe = multiply-expand, a point-by-point multiply of two real-valued
//        buffers with the result placed in the real part of each element of
//        a complex output buffer, imaginary part left unset
//    remult = point-by-point multiply of a complex-valued input buffer with
//        the real part of a second complex input buffer
//    immult = point-by-point multiply of a complex-valued input buffer with
//        the imaginary part of a second complex input buffer
//    multSum = point-by-point multiply of the values in two input buffers and
//        sum with the values of a third input buffer
//    multcSum = point-by-point complex-conjugate multiply of the values in
//        two input buffers (second input conjugated) and sum with the values
//        of a third input buffer
//    multMult = point-by-point multiply of two complex-valued buffers with
//        the result multiplied by the values of a third real-valued buffer
//    multDiv = point-by-point multiply of a complex-valued buffer with a
//        real-valued buffer, with the result divided by another real-valued
//        buffer
//    multDivr = same as multDiv but uses a reciprocal approximation insted of
//        true division
//    multDivnr = same as multDiv but uses a reciprocal approximation with 1
//        Newton-Raphson stage for refinement instead of true division, within
//        2 ulp of multDiv with AVX512 (equal ~99.8% of the time), 3 ulp with
//        AVX2 and 4 ulp with AVX (within 1 ulp ~97% and ~92% of the time)
//
//  Division:
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
//    dotpc = conjugate dot-product of two input buffers, second input
//        conjugated
//
//  Power/Magnitude:
//    pow2 = squares each element of an input buffer
//    pow4 = raises each element of an input buffer to the 4th power
//    pow8 = raises each element of an input buffer to the 8th power
//    mag = magnitude of input buffer
//    magr = magnitude of a complex input buffer using a reciprocal square
//        root approximation (trades accuracy for speed)
//    magSq = magnitude squared of a complex input buffer, faster than mag as
//        it doesn't have to compute a square root
//    normSq = compute the Euclidean norm squared of an input buffer or the
//        product of two input buffers
//
//  Conversion:
//    cnvrtInt8ToFloat = convert signed 8-bit integers to 32-bit floats
//    cnvrtInt16ToFloat = convert signed 16-bit integers to 32-bit floats
//    cnvrtInt32ToFloat = convert signed 32-bit integers to 32-bit floats
//    polToCart = convert a buffer of complex values from polar to cartesian
//        coordinate system, phase in cycles
//    polToCartRad = same as polToCart, phase in radians
//    cartToPol = convert a buffer of complex values from cartesian to polar
//        coordinate system, phase in cycles
//    cartToPolRad = same as cartToPol, phase in radians
//
//  Specialized:
//    tune = tunes complex data with a numerically controlled oscillator (NCO)
//    fs2Shift = Fs/2 frequency shift for complex data
//    agc = run an automatic gain control over an input buffer, adaptively
//        adjusting its magnitude squared to try and maintain a constant value
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
//    4. The magnitude of complex functions that use a fast reciprocal square
//       root approximation have a relative error up to about 6e-5 with AVX512
//       and 6e-4 with AVX/AVX2.  For both the regular magnitude (mag) and
//       reciprocal square root approximation (magr) complex methods the
//       results are only reliable while |in| is between about 1.1e-19 and
//       1.8e19: smaller values lose precision or give 0 (|in| = 0 always
//       gives exactly 0) and larger ones give INF with AVX/AVX2 or NaN with
//       AVX512 due to computing |in|^2 first then computing the magnitude from
//       that.
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

#define _USE_MATH_DEFINES // Must be defined BEFORE including cmath for M_PI
#include <cmath>     // for M_PI
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
#include "evm_trig.h"

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

// check for SVML and only enable if desired
#if ((defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)) && !defined(DISABLE_SVML))
  // icc (classic) / icx, dpcpp (LLVM-based)
  #define EVM_SVML
#endif

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
  #define M_TWO_PI (2*M_PI)
  #define RAD_TO_CYC (1.0/M_TWO_PI)
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

  // reverse mask table for the AVX/AVX2 kernels that handle the final partial
  // register by loading the last full register (overlapping the one before),
  // rows can be cast as __m256i* directly and row n selects the last n floats
  alignas(32) inline constexpr int32_t rmasks[8][8] = {
    { 0, 0, 0, 0, 0, 0, 0, 0},   // remainder 0 (unused)
    { 0, 0, 0, 0, 0, 0, 0,-1},   // remainder 1
    { 0, 0, 0, 0, 0, 0,-1,-1},   // remainder 2
    { 0, 0, 0, 0, 0,-1,-1,-1},   // remainder 3
    { 0, 0, 0, 0,-1,-1,-1,-1},   // remainder 4
    { 0, 0, 0,-1,-1,-1,-1,-1},   // remainder 5
    { 0, 0,-1,-1,-1,-1,-1,-1},   // remainder 6
    { 0,-1,-1,-1,-1,-1,-1,-1}    // remainder 7
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

  //==================//
  // Value Comparison //
  //==================//
  // Compare the individual values of two buffers, placing the minimum in the
  // output buffer
  //    in1 = first input buffer to compare
  //    in2 = second input buffer to compare
  //    len = number of elements to compare
  //    out = minimum of in1 and in2
  inline void min(float const * const in1, float const * const in2,
      const int len, float * const out);
  
  // Compare the individual values of two buffers, placing the maximum in the
  // output buffer
  //    in1 = first input buffer to compare
  //    in2 = second input buffer to compare
  //    len = number of elements to compare
  //    out = maximum of in1 and in2
  inline void max(float const * const in1, float const * const in2,
      const int len, float * const out);
  
  // Compute both the minimum and maximum values of an input buffer
  //    in = buffer to compute max and min over
  //    len = number of elements to compute max and min over
  //    minOut = minimum value
  //    maxOut = maximum value
  //    noInit = if true, minOut and maxOut will not be initialized to the
  //        first value in the input buffer and will instead retain their input
  //        values at the start of the comparison
  inline void minMax(float const * const in, const int len, float &minOut,
      float &maxOut, const bool noInit);

  // Compare the individual elements of two input buffers, placing the min and
  // max of each pair in the respective output buffers
  //    in1 = first input buffer to compare
  //    in2 = second input buffer to compare
  //    len = number of elements to compare
  //    minOut = minimum of in1 and in2
  //    maxOut = maximum of in1 and in2
  inline void minMax(float const * const in1, float const * const in2,
      const int len, float * const minOut, float * const maxOut);

  // Set negative values of a buffer to zero
  //    in = buffer to set negative values to zero
  //    len = number of elements to process
  inline void zneg(float * const in, const int len);

  // Set positive values of a buffer to zero
  //    in = buffer to set positive values to zero
  //    len = number of elements to process
  inline void zpos(float * const in, const int len);
  
  // Compute the minimum of the absolute value of each element across a series
  // of buffers, with the original (signed) value placed in the output buffer
  //    in = buffer of pointers to input buffers to compare, this buffer should
  //        contain at least nbufs pointers which each point to a buffer that
  //        contains at least len samples
  //    nbufs = number of pointers (buffers) in input pointer
  //    len = number of elements to process in each of the nbufs input buffers
  //    out = output buffer to hold the element that had the min absolute value
  //        at each sample location across the nbufs buffers, must have room
  //        for len output samples
  inline void signedAbsMin(float const * const * const in, const int nbufs,
      const int len, float * const out);

  // Compute the maximum of the absolute value of each element across a series
  // of buffers, with the original (signed) value placed in the output buffer
  //    in = buffer of pointers to input buffers to compare, this buffer should
  //        contain at least nbufs pointers which each point to a buffer that
  //        contains at least len samples
  //    nbufs = number of pointers (buffers) in input pointer
  //    len = number of elements to process in each of the nbufs input buffers
  //    out = output buffer to hold the element that had the max absolute value
  //        at each sample location across the nbufs buffers, must have room
  //        for len output samples
  inline void signedAbsMax(float const * const * const in, const int nbufs,
      const int len, float * const out);

  // Given a threshold and two input buffers, compares elements of the first
  // input buffer with the threshold, replacing any values less than the
  // threshold with the value from the second input buffer and putting the
  // result in the output buffer (out[i] = in1[i]>=thresh? in1[i]:in2[i])
  //    in1 = input buffer to compare with threshold
  //    in2 = input buffer to place in output if the corresponding value of in1
  //        is less than the threshold
  //    len = number of elements to compare
  //    thresh = threshold to compare against
  //    out = output buffer to hold values of in1 (if greater than or equal to
  //        the threshold) or values of in2 (if in1 less than the threshold)
  inline void threshCmp(float const * const in1, float const * const in2,
      const int len, const float thresh, float * const out);

  // Places the values of either in1 or in2 into an output buffer depending on
  // whether the values in a third input buffer are greater than or less than
  // a threshold
  //    in1 = buffer of values to use for output if in3 >= threshold
  //    in2 = buffer of values to use for output if in3 < threshold
  //    in3 = buffer of values to compare with threshold
  //    len = number of values to process
  //    thresh = threshold to compare against
  //    out = output buffer to hold values of in1 or in2 depending on if in3
  //        is greater than or less than threshold
  inline void blend(float const * const in1, float const * const in2,
      float const * const in3, const int len, const float thresh,
      float * const out);

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
  inline void scale(float const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out);
  inline void scale(complex<float> const * const in1,
      const complex<float> scaleFactor, const int len,
      complex<float> * const out);

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
  
  //==================================//
  // Vector Multiply and Scale/Expand //
  //==================================//
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
  
  // Point by point multiply of two real-valued buffers with the result placed
  // in the real part of each element of a complex output buffer, imaginary
  // part left unset.  The output buffer may be the same as one of the input
  // buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply
  //    len = number of elements to multiply
  //    out = output buffer to recieve result of multiply in the real part of
  //        each element, imaginary part unmodified
  inline void multe(float const * const in1, float const * const in2,
      const int len, complex<float> * const out);
  
  // Point by point multiply of a complex-valued input buffer with the real
  // part of a second complex input buffer (imaginary part of the second buffer
  // is ignored).  The output buffer may be the same as one of the input
  // buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply, only real part used
  //    len = number of elements to multiply
  //    out = output buffer to recieve result of multiply
  inline void remult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);
  
  // Point by point multiply of a complex-valued input buffer with the
  // imaginary part of a second complex input buffer (real part of the second
  // buffer is ignored).  The output buffer may be the same as one of the input
  // buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply, only imaginary part used
  //    len = number of elements to multiply
  //    out = output buffer to recieve result of multiply
  inline void immult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out);


  //=============================//
  // Three-Input Vector Multiply //
  //=============================//
  // Point by point multiply of two input buffers with and sum with a third
  // input buffer.  The output buffer may be the same as one of the input
  // buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply
  //    in3 = third input buffer to sum with result of multiply
  //    len = number of elements to process
  //    out = output buffer to recieve result of in1*in2+in3
  inline void multSum(float const * const in1, float const * const in2,
      float const * const in3, const int len, float * const out);
  inline void multSum(complex<float> const * const in1, float const * const in2,
      complex<float> const * const in3, const int len,
      complex<float> * const out);
  inline void multSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out);
  
  // Point by point complex-conjugate multiply of two input buffers (second
  // input conjugated) with and sum with a third input buffer.  The output
  // buffer may be the same as one of the input buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply, conjugated
  //    in3 = third input buffer to sum with result of multiply
  //    len = number of elements to process
  //    out = output buffer to recieve result of in1*conj(in2)+in3
  inline void multcSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out);
  
  // Point by point multiply of two complex-valued buffers with the result
  // multiplied by the values of a third real-valued buffer.  The output buffer
  // may be the same as one of the input buffers
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply
  //    in3 = third input buffer for multiply
  //    len = number of elements to process
  //    out = output buffer to recieve result of in1*in2*in3
  inline void multMult(complex<float> const * const in1,
      complex<float> const * const in2, float const * const in3, const int len,
      complex<float> * const out);
  
  // Point by point multiply of a complex-valued buffer and a real-valued
  // buffer with the result divided by the values of a third real-valued
  // buffer.  The output buffer may be the same as one of the input buffers
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply
  //    in3 = third input buffer, divisor
  //    len = number of elements to process
  //    out = output buffer to recieve result of in1*in2/in3
  inline void multDiv(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out);
  // Alternate implementation using a fast reciprocal approximation and
  // multiply rather than a true divide
  inline void multDivr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out);
  // Alternate implementation using a fast reciprocal approximation with 1
  // Newton-Raphson stage for refinement to achieve an answer comparable to
  // true division but with increased efficiency.  Note that due to the
  // reciprocal approximation and Newton-Raphson stage, this method may
  // produce different results than normal division when the divisor is a
  // special floating point value (0 yields NaN rather than INF, INF yields
  // NaN rather than 0, and subnormals may yield NaN or -+INF with the wrong
  // sign).
  inline void multDivnr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
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
  // Dot-product of two input buffers, if len <= 0, out will be set to 0
  //    in1 = first input buffer for dot-product
  //    in2 = second input buffer for dot-product
  //    len = number of elements to compute dot-product over
  //    out = dot-product of in1 and in2
  inline void dotp(float const * const in1, float const * const in2,
      const int len, float &out);
  inline void dotp(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out);
  inline void dotp(float const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out);
  // Wrappers for dotp functions that return by value, these may be slightly
  // slower than the reference versions due to the extra copy
  inline float dotp(float const * const in1, float const * const in2,
      const int len);
  inline complex<float> dotp(complex<float> const * const in1,
      complex<float> const * const in2, const int len);
  inline complex<float> dotp(float const * const in1,
      complex<float> const * const in2, const int len);
  
  // Point by point complex-conjugate dot-product of two input buffers, second
  // input conjugated, if len <= 0, out will be set to 0
  //    in1 = first input buffer for dot-product
  //    in2 = second input buffer for dot-product, conjugated
  //    len = number of elements to compute dot-product over
  //    out = dot-product of in1 and conj(in2)
  inline void dotpc(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out);
  inline void dotpc(float const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out);
  // Wrappers for dotpc functions that return by value, these may be slightly
  // slower than the reference versions due to the extra copy
  inline complex<float> dotpc(complex<float> const * const in1,
      complex<float> const * const in2, const int len);
  inline complex<float> dotpc(float const * const in1,
      complex<float> const * const in2, const int len);

  //================//
  // Raise to Power //
  //================//
  // Square each input sample.  The output buffer may be the same as the input
  // buffer.
  //    in = input buffer to square
  //    len = number of elements to process
  //    out = output buffer to recieve result of in^2
  inline void pow2(float const * const in, const int len, float * const out);
  inline void pow2(complex<float> const * const in, const int len,
      complex<float> * const out);
  
  // Raise each input sample to the 4th power.  The output buffer may be the
  // same as the input buffer.
  //    in = input buffer to raise to 4th power
  //    len = number of elements to process
  //    out = output buffer to recieve result of in^4
  inline void pow4(float const * const in, const int len, float * const out);
  inline void pow4(complex<float> const * const in, const int len,
      complex<float> * const out);

  // Raise each input sample to the 8th power.  The output buffer may be the
  // same as the input buffer.
  //    in = input buffer to raise to 8th power
  //    len = number of elements to process
  //    out = output buffer to recieve result of in^8
  inline void pow8(float const * const in, const int len, float * const out);
  inline void pow8(complex<float> const * const in, const int len,
      complex<float> * const out);

  //=============================//
  // Magnitude/Magnitude Squared //
  //=============================//
  // Compute the magnitude squared of each input sample.  The output buffer
  // may be the same as the input buffer.
  //    in = input buffer to compute magnitude squared of
  //    len = number of elements to process
  //    out = output buffer to recieve result of |in|^2
  inline void magSq(complex<float> const * const in, const int len,
      float * const out);
  // convenience wrapper for magnitude squared of real valued data
  inline void magSq(float const * const in, const int len, float * const out);

  // Compute the magnitude of each input sample.  The output buffer may be the
  // same as the input buffer.  For complex input |in|^2 is computed first, so
  // results are only reliable while |in| is between about 1.1e-19 and 1.8e19:
  // smaller values lose precision or give 0 and larger ones give INF.
  //    in = input buffer to compute magnitude of
  //    len = number of elements to process
  //    out = output buffer to recieve result of |in|
  inline void mag(float const * const in, const int len, float * const out);
  inline void mag(complex<float> const * const in, const int len,
      float * const out);
  // Convenience wrapper for magnitude of real
  inline void abs(float const * const in, const int len, float * const out);
  
  // Compute the magnitude of each input sample using a fast reciprocal square
  // root approximation then multiplying magnitude squared by this
  // approximation.  The output buffer may be the same as the input buffer.
  // The relative error is up to about 6e-5 with AVX512 and 6e-4 with AVX/AVX2.
  // Results are only reliable while |in| is between about 1.1e-19 and 1.8e19:
  // smaller values lose precision or give 0 (|in| = 0 always gives exactly 0)
  // and larger ones give INF with AVX/AVX2 or NaN with AVX512.
  //    in = input buffer to compute magnitude of
  //    len = number of elements to process
  //    out = output buffer to recieve result of |in|
  inline void magr(complex<float> const * const in, const int len,
      float * const out);

  //========================//
  // Euclidean Norm Squared //
  //========================//
  // Euclidean norm squared of an input buffer, if len <= 0, out will be set
  // to 0
  //    in = buffer to compute norm squared of
  //    len = number of elements to compute over
  //    out = Euclidean norm squared of input
  inline void normSq(float const * const in, const int len, float &out);
  inline void normSq(complex<float> const * const in, const int len,
      float &out);
  // Wrappers for normSq functions that return by value, these may be slightly
  // slower than the reference versions due to the extra copy
  inline float normSq(float const * const in, const int len);
  inline float normSq(complex<float> const * const in, const int len);

  //===================================//
  // Unit/Coordinate System Conversion //
  //===================================//
  // Convert 8-bit integer elements to 32-bit floating point
  //    in = input buffer of 8-bit integers to convert
  //    len = number of elements to convert
  //    out = output buffer of floating point elements
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out);
  // Wrapper that uses a common buffer for input and output
  inline void cnvrtInt8ToFloat(void * const inout, const int len);
  
  // Overloaded function that returns min and max elements of the buffer in
  // addition to converting to float
  //    in = input buffer of 8-bit integers to convert
  //    len = number of elements to convert
  //    out = output buffer of floating point elements
  //    minOut = minimum of all output values
  //    maxOut = maximum of all output values
  //    noInit = if true, minOut and maxOut will not be initialized to the
  //        first value in the input buffer and will instead retain their input
  //        values at the start of the comparison
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit);

  // Convert 16-bit integer elements to 32-bit floating point
  //    in = input buffer of 16-bit integers to convert
  //    len = number of elements to convert
  //    out = output buffer of floating point elements
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out);
  // Wrapper that uses a common buffer for input and output
  inline void cnvrtInt16ToFloat(void * const inout, const int len);
  
  // Overloaded function that returns min and max elements of the buffer in
  // addition to converting to float
  //    in = input buffer of 16-bit integers to convert
  //    len = number of elements to convert
  //    out = output buffer of floating point elements
  //    minOut = minimum of all output values
  //    maxOut = maximum of all output values
  //    noInit = if true, minOut and maxOut will not be initialized to the
  //        first value in the input buffer and will instead retain their input
  //        values at the start of the comparison
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit);

  // Convert 32-bit integer elements to 32-bit floating point
  //    in = input buffer of 32-bit integers to convert
  //    len = number of elements to convert
  //    out = output buffer of floating point elements
  inline void cnvrtInt32ToFloat(int32_t const * const in, const int len,
      float * const out);
  // Wrapper that uses a common buffer for input and output
  inline void cnvrtInt32ToFloat(void * const inout, const int len);

  // Overloaded function that returns min and max elements of the buffer in
  // addition to converting to float
  //    in = input buffer of 32-bit integers to convert
  //    len = number of elements to convert
  //    out = output buffer of floating point elements
  //    minOut = minimum of all output values
  //    maxOut = maximum of all output values
  //    noInit = if true, minOut and maxOut will not be initialized to the
  //        first value in the input buffer and will instead retain their input
  //        values at the start of the comparison
  inline void cnvrtInt32ToFloat(int32_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit);

  // Convert a buffer of complex values from polar (mag, phase) to cartesian
  // (x, y) coordinate system.  Phase is expected to be in cycles, output
  // buffer may be same as input buffer.  Whole cycles are removed from the
  // phase, so large phases keep their accuracy.
  //    in = input buffer of polar coordinates to convert
  //    len = number of elements to process
  //    out = output buffer of cartesian coordinates
  inline void polToCart(complex<float> const * const in, const int len,
      complex<float> * const out);

  // Convert a buffer of complex values from polar (mag, phase) to cartesian
  // (x, y) coordinate system.  Phase is expected to be in radians, output
  // buffer may be same as input buffer.
  //    in = input buffer of polar coordinates to convert
  //    len = number of elements to process
  //    out = output buffer of cartesian coordinates
  inline void polToCartRad(complex<float> const * const in, const int len,
      complex<float> * const out);

  // Convert a buffer of complex values from cartesian (x, y) to polar
  // (mag, phase) coordinate system, returned phase is in cycles, output
  // buffer may be same as input buffer.
  //    in = input buffer of cartesian coordinates to convert
  //    len = number of elements to process
  //    out = output buffer of polar coordinates
  inline void cartToPol(complex<float> const * const in, const int len,
      complex<float> * const out);

  // Convert a buffer of complex values from cartesian (x, y) to polar
  // (mag, phase) coordinate system, returned phase is in radians, output
  // buffer may be same as input buffer.
  //    in = input buffer of cartesian coordinates to convert
  //    len = number of elements to process
  //    out = output buffer of polar coordinates
  inline void cartToPolRad(complex<float> const * const in, const int len,
      complex<float> * const out);

  //=======================//
  // Specialized Functions //
  //=======================//
  // Multiply complex input data with a numerically controlled oscillator
  // (NCO) to induce a frequency shift.  The NCO is implemented via creating
  // a complex exponential of the starting phase and multiplying it by a
  // complex exponential with the per-sample phase difference on each sample.
  // Since initial phase exponential is only computed once at the start of this
  // function, this can produce a sawtooth effect on the phase will will slowly
  // deviate from the true phase of the complex exponential due to minute
  // rounding differences in the complex multiply vs. a true phase accumulator.
  // The function will return the starting phase that should be used on the
  // next call to the function for continuity in the NCO.  The returned phase
  // is computed independently from the NCO using a double-precision multiply
  // and modulo division, so it will be accurate even if the NCO slowly drifts
  // away from the desired phase.
  //    in = input buffer of data to tune
  //    len = number of elements to process
  //    phInit = initial phase of the tuner for this buffer in cycles
  //    phIncr = phase increment per sample in cycles
  //    out = output buffer of tuned data
  //    EVM::tune = initial phase to use for next call this method in cycles
  inline double tune(complex<float> const * const in, const int len,
      const float phInit, const float phIncr, complex<float> * const out);
  inline double tune(float const * const in, const int len,
      const float phInit, const float phIncr, complex<float> * const out);
  // Overloaded method that also returns the complex phasor values used to tune
  //    in = input buffer of data to tune
  //    len = number of elements to process
  //    phInit = initial phase of the tuner for this buffer in cycles
  //    phIncr = phase increment per sample in cycles
  //    out = output buffer of tuned data
  //    phOut = complex phasor values used to tune the input data
  //    EVM::tune = initial phase to use for next call this method in cycles
  inline double tune(complex<float> const * const in, const int len,
      const float phInit, const float phIncr, complex<float> * const out,
      complex<float> * const phOut);
  inline double tune(float const * const in, const int len,
      const float phInit, const float phIncr, complex<float> * const out,
      complex<float> * const phOut);

  // Implements a Fs/2 frequency shift (where Fs is the sampling rate) for
  // complex data
  //    in = input buffer to apply frequency shift to
  //    len = number of elements to process
  //    out = output buffer of shifted data
  //    init = if set to 0, start the phase shift on 0 cycles, if set to 1 (or
  //        any other value) set the first phase shift to .5 cycles
  //    EVM::fs2shift = initialization value for next call to this method
  inline int fs2Shift(complex<float> const * const in, const int len,
      complex<float> * const out, const int init);

  // Run an automatic gain control over an input buffer, adaptively adjusting
  // its magnitude squared to keep it as close to the specified constant value
  // as possible.
  //    in = input buffer of values to apply AGC to
  //    len = number of values to process
  //    init = initial gain to apply (returned value from previous call or 1
  //        if this is the first call)
  //    desired = desired average magnitude squared value to achieve
  //    out = output buffer of values after AGC applied
  //    EVM::agc = initial gain to apply to next input buffer
  inline float agc(float const * const in, const int len, const float init,
      const float desired, float * const out);
  inline float agc(complex<float> const * const in, const int len,
      const float init, const float desired, complex<float> * const out);
  // Overloaded methods that also output the gain values applied
  //    in = input buffer of values to apply AGC to
  //    len = number of values to process
  //    init = initial gain to apply (returned value from previous call or 1
  //        if this is the first call)
  //    desired = desired average magnitude squared value to achieve
  //    out = output buffer of values after AGC applied
  //    gains = output buffer of gain values applied
  //    EVM::agc = initial gain to apply to next input buffer
  inline float agc(float const * const in, const int len, const float init,
      const float desired, float * const out, float * const gains);
  inline float agc(complex<float> const * const in, const int len,
      const float init, const float desired, complex<float> * const out,
      float * const gains);

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
  
  //==================//
  // Value Comparison //
  //==================//
  // compare each element of two input buffers, output min
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void min(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    int i = 0;
    // note: min_ps returns the second operand on a tie or NaN, so in2 goes
    // first to match std::min(in1, in2) in the default version
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld1 = _mm512_min_ps(ld2, ld1);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld1 = _mm512_min_ps(ld2, ld1);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 min
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void min(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2;
    int i = 0;
    // note: min_ps returns the second operand on a tie or NaN, so in2 goes
    // first to match std::min(in1, in2) in the default version
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_min_ps(ld2, ld1);
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
      ld1 = _mm256_min_ps(ld2, ld1);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX min
  __attribute__((__target__("default")))
  inline void min(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = std::min(in1[i], in2[i]);
    return;
  }

  // compare each element of two input buffers, output max
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void max(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    int i = 0;
    // note: max_ps returns the second operand on a tie or NaN, so in2 goes
    // first to match std::max(in1, in2) in the default version
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld1 = _mm512_max_ps(ld2, ld1);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld1 = _mm512_max_ps(ld2, ld1);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 max
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void max(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2;
    int i = 0;
    // note: max_ps returns the second operand on a tie or NaN, so in2 goes
    // first to match std::max(in1, in2) in the default version
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_max_ps(ld2, ld1);
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
      ld1 = _mm256_max_ps(ld2, ld1);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX max
  __attribute__((__target__("default")))
  inline void max(float const * const in1, float const * const in2,
      const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = std::max(in1[i], in2[i]);
    return;
  }
  
  // determine minimum and maximum value of an input buffer
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void minMax(float const * const in, const int len, float &minOut,
      float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    __m512 ld;
    // start from the current min/max (noInit) or the first element
    __m512 mn = _mm512_set1_ps(noInit? minOut : in[0]);
    __m512 mx = _mm512_set1_ps(noInit? maxOut : in[0]);
    int i = 0;
    // note: the loaded value is the first operand so that a NaN element is
    // skipped (min_ps/max_ps return the second operand on a NaN)
    // note: with a single min and max accumulator every min/max waits on the
    // previous one, instead use 4 independent pairs to keep the pipeline full
    if(len > 63)
    {
      __m512 mn1 = mn;
      __m512 mn2 = mn;
      __m512 mn3 = mn;
      __m512 mx1 = mx;
      __m512 mx2 = mx;
      __m512 mx3 = mx;
      for(; i<len-63; i+=64) // process 4 registers (64 real elements)
      {
        ld = _mm512_loadu_ps(&in[i]);
        mn = _mm512_min_ps(ld, mn);
        mx = _mm512_max_ps(ld, mx);
        ld = _mm512_loadu_ps(&in[i+16]);
        mn1 = _mm512_min_ps(ld, mn1);
        mx1 = _mm512_max_ps(ld, mx1);
        ld = _mm512_loadu_ps(&in[i+32]);
        mn2 = _mm512_min_ps(ld, mn2);
        mx2 = _mm512_max_ps(ld, mx2);
        ld = _mm512_loadu_ps(&in[i+48]);
        mn3 = _mm512_min_ps(ld, mn3);
        mx3 = _mm512_max_ps(ld, mx3);
      }
      mn = _mm512_min_ps(mn, mn1);
      mn2 = _mm512_min_ps(mn2, mn3);
      mn = _mm512_min_ps(mn, mn2);
      mx = _mm512_max_ps(mx, mx1);
      mx2 = _mm512_max_ps(mx2, mx3);
      mx = _mm512_max_ps(mx, mx2);
    }
    // 0-3 full registers remain and up to 1 partial register, compute these
    // individually as using a loop will cause intel compiler to save/restore
    // additional registers
    // bit 5 indicates there are at least 32 elements that weren't taken care
    // of by the loop above
    if(len&32) // process 2 registers (32 real elements)
    {
      ld = _mm512_loadu_ps(&in[i]);
      mn = _mm512_min_ps(ld, mn);
      mx = _mm512_max_ps(ld, mx);
      ld = _mm512_loadu_ps(&in[i+16]);
      mn = _mm512_min_ps(ld, mn);
      mx = _mm512_max_ps(ld, mx);
      i += 32;
    }
    // bit 4 indicates there were an additional 16 elements not taken care of
    // by the loop above
    if(len&16) // 1 register (16 real elements)
    {
      ld = _mm512_loadu_ps(&in[i]);
      mn = _mm512_min_ps(ld, mn);
      mx = _mm512_max_ps(ld, mx);
      i += 16;
    }
    // handle remaining elements in partial register (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld = _mm512_maskz_loadu_ps(mk, &in[i]);
      // only compare the loaded elements, the others keep their min/max
      mn = _mm512_mask_min_ps(mn, mk, ld, mn);
      mx = _mm512_mask_max_ps(mx, mk, ld, mx);
    }
    minOut = _mm512_reduce_min_ps(mn);
    maxOut = _mm512_reduce_max_ps(mx);
    return;
  }
  #endif // AVX512 minMax element
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void minMax(float const * const in, const int len, float &minOut,
      float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    __m256 ld;
    // start from the current min/max (noInit) or the first element
    __m256 mn = _mm256_set1_ps(noInit? minOut : in[0]);
    __m256 mx = _mm256_set1_ps(noInit? maxOut : in[0]);
    __m128 hs;  // horizontal min/max to finish
    int i = 0;
    // note: the loaded value is the first operand so that a NaN element is
    // skipped (min_ps/max_ps return the second operand on a NaN)
    // note: with a single min and max accumulator every min/max waits on the
    // previous one, instead use 4 independent pairs to keep the pipeline full
    if(len > 31)
    {
      __m256 mn1 = mn;
      __m256 mn2 = mn;
      __m256 mn3 = mn;
      __m256 mx1 = mx;
      __m256 mx2 = mx;
      __m256 mx3 = mx;
      for(; i<len-31; i+=32) // process 4 registers (32 real elements)
      {
        ld = _mm256_loadu_ps(&in[i]);
        mn = _mm256_min_ps(ld, mn);
        mx = _mm256_max_ps(ld, mx);
        ld = _mm256_loadu_ps(&in[i+8]);
        mn1 = _mm256_min_ps(ld, mn1);
        mx1 = _mm256_max_ps(ld, mx1);
        ld = _mm256_loadu_ps(&in[i+16]);
        mn2 = _mm256_min_ps(ld, mn2);
        mx2 = _mm256_max_ps(ld, mx2);
        ld = _mm256_loadu_ps(&in[i+24]);
        mn3 = _mm256_min_ps(ld, mn3);
        mx3 = _mm256_max_ps(ld, mx3);
      }
      mn = _mm256_min_ps(mn, mn1);
      mn2 = _mm256_min_ps(mn2, mn3);
      mn = _mm256_min_ps(mn, mn2);
      mx = _mm256_max_ps(mx, mx1);
      mx2 = _mm256_max_ps(mx2, mx3);
      mx = _mm256_max_ps(mx, mx2);
    }
    // 0-3 full registers remain and up to 1 partial register, compute these
    // individually as using a loop will cause intel compiler to save/restore
    // additional registers
    // bit 4 indicates there are at least 16 elements that weren't taken care
    // of by the loop above
    if(len&16) // 2 registers (16 real elements)
    {
      ld = _mm256_loadu_ps(&in[i]);
      mn = _mm256_min_ps(ld, mn);
      mx = _mm256_max_ps(ld, mx);
      ld = _mm256_loadu_ps(&in[i+8]);
      mn = _mm256_min_ps(ld, mn);
      mx = _mm256_max_ps(ld, mx);
      i += 16;
    }
    // bit 3 indicates there were an additional 8 elements not taken care of
    // by the loop above
    if(len&8) // 1 register (8 real elements)
    {
      ld = _mm256_loadu_ps(&in[i]);
      mn = _mm256_min_ps(ld, mn);
      mx = _mm256_max_ps(ld, mx);
      i += 8;
    }
    // handle remaining elements in partial register (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld = _mm256_maskload_ps(&in[i], msk);
      // maskload sets the unused elements to 0, which could become the min or
      // max, so set them to in[i] instead (always loaded)
      ld = _mm256_or_ps(ld, _mm256_andnot_ps(_mm256_castsi256_ps(msk),
          _mm256_broadcast_ss(&in[i])));
      mn = _mm256_min_ps(ld, mn);
      mx = _mm256_max_ps(ld, mx);
    }
    // horizontal min and max
    hs = _mm256_extractf128_ps(mn, 1); // extract upper half of register
    hs = _mm_min_ps(_mm256_castps256_ps128(mn), hs); // [min(0,4), min(1,5),..]
    hs = _mm_min_ps(hs, _mm_movehl_ps(hs, hs)); // [min(0,4,2,6), ...]
    hs = _mm_min_ss(hs, _mm_movehdup_ps(hs));   // [min(0,...,7), ...]
    minOut = _mm_cvtss_f32(hs); // extract the first element in the register
    hs = _mm256_extractf128_ps(mx, 1); // extract upper half of register
    hs = _mm_max_ps(_mm256_castps256_ps128(mx), hs); // [max(0,4), max(1,5),..]
    hs = _mm_max_ps(hs, _mm_movehl_ps(hs, hs)); // [max(0,4,2,6), ...]
    hs = _mm_max_ss(hs, _mm_movehdup_ps(hs));   // [max(0,...,7), ...]
    maxOut = _mm_cvtss_f32(hs); // extract the first element in the register
    return;
  }
  #endif // AVX minMax element
  __attribute__((__target__("default")))
  inline void minMax(float const * const in, const int len, float &minOut,
      float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    if(!noInit)
    {
      minOut = in[0];
      maxOut = in[0];
    }
    for(int i=0; i<len; ++i)
    {
      minOut = std::min(minOut, in[i]);
      maxOut = std::max(maxOut, in[i]);
    }
    return;
  }
  
  // compare each element of two input buffers, output min and max in separate
  // buffers
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void minMax(float const * const in1, float const * const in2,
      const int len, float * const minOut, float * const maxOut)
  {
    if(len<=0) return;
    // process 4 registers at once for each input buffer to pipeline stores
    __m512 ld1a, ld1b, ld1c, ld1d;
    __m512 ld2a, ld2b, ld2c, ld2d;
    __m512 mn;  // min value
    int i = 0;
    // note: min_ps/max_ps return the second operand on a tie or NaN, so in2
    // goes first to match std::min/max(in1, in2) in the default version
    // note: storing one register at a time to minOut then maxOut is about 2.5x
    // slower when the outputs are not 64-byte aligned and don't fit in L1
    // cache, instead store 4 registers to minOut and then 4 to maxOut
    for(; i<len-63; i+=64) // process 4 registers (64 real elements)
    {
      // load all 4 registers first in case minOut or maxOut is the same
      // buffer as in1 or in2
      ld1a = _mm512_loadu_ps(&in1[i]);
      ld1b = _mm512_loadu_ps(&in1[i+16]);
      ld1c = _mm512_loadu_ps(&in1[i+32]);
      ld1d = _mm512_loadu_ps(&in1[i+48]);
      ld2a = _mm512_loadu_ps(&in2[i]);
      ld2b = _mm512_loadu_ps(&in2[i+16]);
      ld2c = _mm512_loadu_ps(&in2[i+32]);
      ld2d = _mm512_loadu_ps(&in2[i+48]);
      _mm512_storeu_ps(&minOut[i], _mm512_min_ps(ld2a, ld1a));
      _mm512_storeu_ps(&minOut[i+16], _mm512_min_ps(ld2b, ld1b));
      _mm512_storeu_ps(&minOut[i+32], _mm512_min_ps(ld2c, ld1c));
      _mm512_storeu_ps(&minOut[i+48], _mm512_min_ps(ld2d, ld1d));
      _mm512_storeu_ps(&maxOut[i], _mm512_max_ps(ld2a, ld1a));
      _mm512_storeu_ps(&maxOut[i+16], _mm512_max_ps(ld2b, ld1b));
      _mm512_storeu_ps(&maxOut[i+32], _mm512_max_ps(ld2c, ld1c));
      _mm512_storeu_ps(&maxOut[i+48], _mm512_max_ps(ld2d, ld1d));
    }
    // 0-3 full registers remain and up to 1 partial register
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1a = _mm512_loadu_ps(&in1[i]);
      ld2a = _mm512_loadu_ps(&in2[i]);
      mn = _mm512_min_ps(ld2a, ld1a);
      ld1a = _mm512_max_ps(ld2a, ld1a);
      _mm512_storeu_ps(&minOut[i], mn);
      _mm512_storeu_ps(&maxOut[i], ld1a);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1a = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2a = _mm512_maskz_loadu_ps(mk, &in2[i]);
      mn = _mm512_min_ps(ld2a, ld1a);
      ld1a = _mm512_max_ps(ld2a, ld1a);
      _mm512_mask_storeu_ps(&minOut[i], mk, mn);
      _mm512_mask_storeu_ps(&maxOut[i], mk, ld1a);
    }
    return;
  }
  #endif // AVX512 minMax buffer
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void minMax(float const * const in1, float const * const in2,
      const int len, float * const minOut, float * const maxOut)
  {
    if(len<=0) return;
    // process 4 registers at once for each input buffer to pipeline stores
    __m256 ld1a, ld1b, ld1c, ld1d;
    __m256 ld2a, ld2b, ld2c, ld2d;
    __m256 mn;  // min value
    int i = 0;
    // note: min_ps/max_ps return the second operand on a tie or NaN, so in2
    // goes first to match std::min/max(in1, in2) in the default version
    // note: storing one register at a time to minOut then maxOut is about 3x
    // slower when the outputs don't fit in L1 cache, instead store 4
    // registers to minOut and then 4 to maxOut
    for(; i<len-31; i+=32) // process 4 registers (32 real elements)
    {
      // load all 4 registers first in case minOut or maxOut is the same
      // buffer as in1 or in2
      ld1a = _mm256_loadu_ps(&in1[i]);
      ld1b = _mm256_loadu_ps(&in1[i+8]);
      ld1c = _mm256_loadu_ps(&in1[i+16]);
      ld1d = _mm256_loadu_ps(&in1[i+24]);
      ld2a = _mm256_loadu_ps(&in2[i]);
      ld2b = _mm256_loadu_ps(&in2[i+8]);
      ld2c = _mm256_loadu_ps(&in2[i+16]);
      ld2d = _mm256_loadu_ps(&in2[i+24]);
      _mm256_storeu_ps(&minOut[i], _mm256_min_ps(ld2a, ld1a));
      _mm256_storeu_ps(&minOut[i+8], _mm256_min_ps(ld2b, ld1b));
      _mm256_storeu_ps(&minOut[i+16], _mm256_min_ps(ld2c, ld1c));
      _mm256_storeu_ps(&minOut[i+24], _mm256_min_ps(ld2d, ld1d));
      _mm256_storeu_ps(&maxOut[i], _mm256_max_ps(ld2a, ld1a));
      _mm256_storeu_ps(&maxOut[i+8], _mm256_max_ps(ld2b, ld1b));
      _mm256_storeu_ps(&maxOut[i+16], _mm256_max_ps(ld2c, ld1c));
      _mm256_storeu_ps(&maxOut[i+24], _mm256_max_ps(ld2d, ld1d));
    }
    // 0-3 full registers remain and up to 1 partial register
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1a = _mm256_loadu_ps(&in1[i]);
      ld2a = _mm256_loadu_ps(&in2[i]);
      mn = _mm256_min_ps(ld2a, ld1a);
      ld1a = _mm256_max_ps(ld2a, ld1a);
      _mm256_storeu_ps(&minOut[i], mn);
      _mm256_storeu_ps(&maxOut[i], ld1a);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1a = _mm256_maskload_ps(&in1[i], msk);
      ld2a = _mm256_maskload_ps(&in2[i], msk);
      mn = _mm256_min_ps(ld2a, ld1a);
      ld1a = _mm256_max_ps(ld2a, ld1a);
      _mm256_maskstore_ps(&minOut[i], msk, mn);
      _mm256_maskstore_ps(&maxOut[i], msk, ld1a);
    }
    return;
  }
  #endif // AVX minMax buffer
  __attribute__((__target__("default")))
  inline void minMax(float const * const in1, float const * const in2,
      const int len, float * const minOut, float * const maxOut)
  {
    float val1;
    float val2;
    for(int i=0; i<len; ++i)
    {
      // copy the inputs first in case minOut is the same buffer as in1 or in2
      val1 = in1[i];
      val2 = in2[i];
      minOut[i] = std::min(val1, val2);
      maxOut[i] = std::max(val1, val2);
    }
    return;
  }
  
  // set negative values to zero
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void zneg(float * const in, const int len)
  {
    if(len<=0) return;
    __m512 ld;
    const __m512 zero = _mm512_set1_ps(0.0f);
    int i = 0;
    // note: 0 goes first so that NaN and -0.0 are kept (max_ps returns the
    // second operand on a tie or NaN), as with std::max(in[i], 0.0f)
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld = _mm512_loadu_ps(&in[i]);
      ld = _mm512_max_ps(zero, ld);
      _mm512_storeu_ps(&in[i], ld);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld = _mm512_maskz_loadu_ps(mk, &in[i]);
      ld = _mm512_max_ps(zero, ld);
      _mm512_mask_storeu_ps(&in[i], mk, ld);
    }
    return;
  }
  #endif // AVX512 zneg
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void zneg(float * const in, const int len)
  {
    if(len<=0) return;
    __m256 ld;
    const __m256 zero = _mm256_set1_ps(0.0f);
    int i = 0;
    // note: 0 goes first so that NaN and -0.0 are kept (max_ps returns the
    // second operand on a tie or NaN), as with std::max(in[i], 0.0f)
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld = _mm256_loadu_ps(&in[i]);
      ld = _mm256_max_ps(zero, ld);
      _mm256_storeu_ps(&in[i], ld);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld = _mm256_maskload_ps(&in[i], msk);
      ld = _mm256_max_ps(zero, ld);
      _mm256_maskstore_ps(&in[i], msk, ld);
    }
    return;
  }
  #endif // AVX zneg
  __attribute__((__target__("default")))
  inline void zneg(float * const in, const int len)
  {
    for(int i=0; i<len; ++i) in[i] = std::max(in[i], 0.0f);
    return;
  }
  
  // set positive elements to zero
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void zpos(float * const in, const int len)
  {
    if(len<=0) return;
    __m512 ld;
    const __m512 zero = _mm512_set1_ps(0.0f);
    int i = 0;
    // note: 0 goes first so that NaN is kept (min_ps returns the second
    // operand on a tie or NaN), as with std::min(in[i], 0.0f)
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld = _mm512_loadu_ps(&in[i]);
      ld = _mm512_min_ps(zero, ld);
      _mm512_storeu_ps(&in[i], ld);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld = _mm512_maskz_loadu_ps(mk, &in[i]);
      ld = _mm512_min_ps(zero, ld);
      _mm512_mask_storeu_ps(&in[i], mk, ld);
    }
    return;
  }
  #endif // AVX512 zpos
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void zpos(float * const in, const int len)
  {
    if(len<=0) return;
    __m256 ld;
    const __m256 zero = _mm256_set1_ps(0.0f);
    int i = 0;
    // note: 0 goes first so that NaN is kept (min_ps returns the second
    // operand on a tie or NaN), as with std::min(in[i], 0.0f)
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld = _mm256_loadu_ps(&in[i]);
      ld = _mm256_min_ps(zero, ld);
      _mm256_storeu_ps(&in[i], ld);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld = _mm256_maskload_ps(&in[i], msk);
      ld = _mm256_min_ps(zero, ld);
      _mm256_maskstore_ps(&in[i], msk, ld);
    }
    return;
  }
  #endif // AVX zpos
  __attribute__((__target__("default")))
  inline void zpos(float * const in, const int len)
  {
    for(int i=0; i<len; ++i) in[i] = std::min(in[i], 0.0f);
    return;
  }
  
  // output elements with min absolute value across N buffers
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void signedAbsMin(float const * const * const in, const int nbufs,
      const int len, float * const out)
  {
    if(nbufs<=0 || len<=0) return;
    // use a 2 register loop for independent compare loops
    __m512 ld1, absLd1, val1, absVal1;
    __m512 ld2, absLd2, val2, absVal2;
    __mmask16 lt1, lt2;
    // constant to make the sign bit positive
    const __m512i pos = _mm512_set1_epi32(0x7fffffff);
    int i = 0;
    // note: a strict less-than so a tie keeps the value from the earlier
    // buffer, as in the default version
    // note: 2 registers per pass over the buffers halves the loop overhead
    // and gives 2 independent compare chains
    for(; i<len-31; i+=32) // process 2 registers (32 real elements)
    {
      val1 = _mm512_loadu_ps(&in[0][i]);
      val2 = _mm512_loadu_ps(&in[0][i+16]);
      // force sign bit positive, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_and_ps (casting should be free)
      absVal1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(val1),
          pos));
      absVal2 = _mm512_castsi512_ps(_mm512_and_si512(
          _mm512_castps_si512(val2), pos));
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm512_loadu_ps(&in[j][i]);
        // load both registers before using either one, otherwise g++ reads
        // ld from memory a second time
        ld2 = _mm512_loadu_ps(&in[j][i+16]);
        absLd1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(ld1),
            pos));
        lt1 = _mm512_cmp_ps_mask(absLd1, absVal1, _CMP_LT_OQ); // |ld| < |val|
        absVal1 = _mm512_mask_mov_ps(absVal1, lt1, absLd1);
        val1 = _mm512_mask_mov_ps(val1, lt1, ld1);
        absLd2 = _mm512_castsi512_ps(_mm512_and_si512(
            _mm512_castps_si512(ld2), pos));
        lt2 = _mm512_cmp_ps_mask(absLd2, absVal2, _CMP_LT_OQ);
        absVal2 = _mm512_mask_mov_ps(absVal2, lt2, absLd2);
        val2 = _mm512_mask_mov_ps(val2, lt2, ld2);
      }
      _mm512_storeu_ps(&out[i], val1);
      _mm512_storeu_ps(&out[i+16], val2);
    }
    // 0-1 full registers remain and up to 1 partial register
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      val1 = _mm512_loadu_ps(&in[0][i]);
      // force sign bit positive, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_and_ps (casting should be free)
      absVal1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(val1),
          pos));
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm512_loadu_ps(&in[j][i]);
        absLd1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(ld1),
            pos));
        lt1 = _mm512_cmp_ps_mask(absLd1, absVal1, _CMP_LT_OQ); // |ld| < |val|
        absVal1 = _mm512_mask_mov_ps(absVal1, lt1, absLd1);
        val1 = _mm512_mask_mov_ps(val1, lt1, ld1);
      }
      _mm512_storeu_ps(&out[i], val1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      val1 = _mm512_maskz_loadu_ps(mk, &in[0][i]);
      // force sign bit positive, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_and_ps (casting should be free)
      absVal1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(val1),
          pos));
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm512_maskz_loadu_ps(mk, &in[j][i]);
        absLd1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(ld1),
            pos));
        lt1 = _mm512_cmp_ps_mask(absLd1, absVal1, _CMP_LT_OQ); // |ld| < |val|
        absVal1 = _mm512_mask_mov_ps(absVal1, lt1, absLd1);
        val1 = _mm512_mask_mov_ps(val1, lt1, ld1);
      }
      _mm512_mask_storeu_ps(&out[i], mk, val1);
    }
    return;
  }
  #endif // AVX512 signedAbsMin
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void signedAbsMin(float const * const * const in, const int nbufs,
      const int len, float * const out)
  {
    if(nbufs<=0 || len<=0) return;
    // use a 2 register loop for independent compare loops
    __m256 ld1, absLd1, val1, absVal1, lt1;
    __m256 ld2, absLd2, val2, absVal2, lt2;
    // constant to make the sign bit positive
    const __m256 pos = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    int i = 0;
    // note: a strict less-than so a tie keeps the value from the earlier
    // buffer, as in the default version
    // note: 2 registers per pass over the buffers halves the loop overhead
    // and gives 2 independent compare chains
    for(; i<len-15; i+=16) // process 2 registers (16 real elements)
    {
      val1 = _mm256_loadu_ps(&in[0][i]);
      val2 = _mm256_loadu_ps(&in[0][i+8]);
      absVal1 = _mm256_and_ps(val1, pos);
      absVal2 = _mm256_and_ps(val2, pos);
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm256_loadu_ps(&in[j][i]);
        // load both registers before using either one, otherwise g++ reads
        // ld from memory a second time
        ld2 = _mm256_loadu_ps(&in[j][i+8]);
        absLd1 = _mm256_and_ps(ld1, pos);
        lt1 = _mm256_cmp_ps(absLd1, absVal1, _CMP_LT_OQ); // |ld| < |val|
        absVal1 = _mm256_min_ps(absLd1, absVal1); // same choice as lt
        // lt? ld : val (avoids _mm256_blendv_ps, see setUnusedToOne)
        val1 = _mm256_or_ps(_mm256_and_ps(lt1, ld1),
            _mm256_andnot_ps(lt1, val1));
        absLd2 = _mm256_and_ps(ld2, pos);
        lt2 = _mm256_cmp_ps(absLd2, absVal2, _CMP_LT_OQ);
        absVal2 = _mm256_min_ps(absLd2, absVal2);
        val2 = _mm256_or_ps(_mm256_and_ps(lt2, ld2),
            _mm256_andnot_ps(lt2, val2));
      }
      _mm256_storeu_ps(&out[i], val1);
      _mm256_storeu_ps(&out[i+8], val2);
    }
    // 0-1 full registers remain and up to 1 partial register
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      val1 = _mm256_loadu_ps(&in[0][i]);
      absVal1 = _mm256_and_ps(val1, pos);
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm256_loadu_ps(&in[j][i]);
        absLd1 = _mm256_and_ps(ld1, pos);
        lt1 = _mm256_cmp_ps(absLd1, absVal1, _CMP_LT_OQ); // |ld| < |val|
        absVal1 = _mm256_min_ps(absLd1, absVal1); // same choice as lt
        // lt1? ld1 : val1 (avoids _mm256_blendv_ps, see setUnusedToOne)
        val1 = _mm256_or_ps(_mm256_and_ps(lt1, ld1),
            _mm256_andnot_ps(lt1, val1));
      }
      _mm256_storeu_ps(&out[i], val1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      val1 = _mm256_maskload_ps(&in[0][i], msk);
      absVal1 = _mm256_and_ps(val1, pos);
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm256_maskload_ps(&in[j][i], msk);
        absLd1 = _mm256_and_ps(ld1, pos);
        lt1 = _mm256_cmp_ps(absLd1, absVal1, _CMP_LT_OQ); // |ld| < |val|
        absVal1 = _mm256_min_ps(absLd1, absVal1); // same choice as lt
        // lt? ld : val (avoids _mm256_blendv_ps, see setUnusedToOne)
        val1 = _mm256_or_ps(_mm256_and_ps(lt1, ld1),
            _mm256_andnot_ps(lt1, val1));
      }
      _mm256_maskstore_ps(&out[i], msk, val1);
    }
    return;
  }
  #endif // AVX signedAbsMin
  __attribute__((__target__("default")))
  inline void signedAbsMin(float const * const * const in, const int nbufs,
      const int len, float * const out)
  {
    if(nbufs<=0 || len<=0) return;
    float val;
    float absVal;
    float tmpVal;
    for(int i=0; i<len; ++i)
    {
      // hold the result in val until every buffer is read in case out is the
      // same buffer as one of the inputs
      val = in[0][i];
      absVal = fabsf(val);
      for(int j=1; j<nbufs; ++j)
      {
        tmpVal = fabsf(in[j][i]);
        if(tmpVal < absVal)
        {
          absVal = tmpVal;
          val = in[j][i];
        }
      }
      out[i] = val;
    }
    return;
  }

  // output elements with max absolute value across N buffers
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void signedAbsMax(float const * const * const in, const int nbufs,
      const int len, float * const out)
  {
    if(nbufs<=0 || len<=0) return;
    // use a 2 register loop for independent compare loops
    __m512 ld1, absLd1, val1, absVal1;
    __m512 ld2, absLd2, val2, absVal2;
    __mmask16 gt1, gt2;
    // constant to make the sign bit positive
    const __m512i pos = _mm512_set1_epi32(0x7fffffff);
    int i = 0;
    // note: a strict greater-than so a tie keeps the value from the earlier
    // buffer, as in the default version
    // note: 2 registers per pass over the buffers halves the loop overhead
    // and gives 2 independent compare chains
    for(; i<len-31; i+=32) // process 2 registers (32 real elements)
    {
      val1 = _mm512_loadu_ps(&in[0][i]);
      val2 = _mm512_loadu_ps(&in[0][i+16]);
      // force sign bit positive, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_and_ps (casting should be free)
      absVal1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(val1),
          pos));
      absVal2 = _mm512_castsi512_ps(_mm512_and_si512(
          _mm512_castps_si512(val2), pos));
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm512_loadu_ps(&in[j][i]);
        // load both registers before using either one, otherwise g++ reads
        // ld from memory a second time
        ld2 = _mm512_loadu_ps(&in[j][i+16]);
        absLd1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(ld1),
            pos));
        gt1 = _mm512_cmp_ps_mask(absLd1, absVal1, _CMP_GT_OQ); // |ld| > |val|
        absVal1 = _mm512_mask_mov_ps(absVal1, gt1, absLd1);
        val1 = _mm512_mask_mov_ps(val1, gt1, ld1);
        absLd2 = _mm512_castsi512_ps(_mm512_and_si512(
            _mm512_castps_si512(ld2), pos));
        gt2 = _mm512_cmp_ps_mask(absLd2, absVal2, _CMP_GT_OQ);
        absVal2 = _mm512_mask_mov_ps(absVal2, gt2, absLd2);
        val2 = _mm512_mask_mov_ps(val2, gt2, ld2);
      }
      _mm512_storeu_ps(&out[i], val1);
      _mm512_storeu_ps(&out[i+16], val2);
    }
    // 0-1 full registers remain and up to 1 partial register
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      val1 = _mm512_loadu_ps(&in[0][i]);
      // force sign bit positive, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_and_ps (casting should be free)
      absVal1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(val1),
          pos));
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm512_loadu_ps(&in[j][i]);
        absLd1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(ld1),
            pos));
        gt1 = _mm512_cmp_ps_mask(absLd1, absVal1, _CMP_GT_OQ); // |ld| > |val|
        absVal1 = _mm512_mask_mov_ps(absVal1, gt1, absLd1);
        val1 = _mm512_mask_mov_ps(val1, gt1, ld1);
      }
      _mm512_storeu_ps(&out[i], val1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      val1 = _mm512_maskz_loadu_ps(mk, &in[0][i]);
      // force sign bit positive, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_and_ps (casting should be free)
      absVal1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(val1),
          pos));
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm512_maskz_loadu_ps(mk, &in[j][i]);
        absLd1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(ld1),
            pos));
        gt1 = _mm512_cmp_ps_mask(absLd1, absVal1, _CMP_GT_OQ); // |ld| > |val|
        absVal1 = _mm512_mask_mov_ps(absVal1, gt1, absLd1);
        val1 = _mm512_mask_mov_ps(val1, gt1, ld1);
      }
      _mm512_mask_storeu_ps(&out[i], mk, val1);
    }
    return;
  }
  #endif // AVX512 signedAbsMax
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void signedAbsMax(float const * const * const in, const int nbufs,
      const int len, float * const out)
  {
    if(nbufs<=0 || len<=0) return;
    // use a 2 register loop for independent compare loops
    __m256 ld1, absLd1, val1, absVal1, gt1;
    __m256 ld2, absLd2, val2, absVal2, gt2;
    // constant to make the sign bit positive
    const __m256 pos = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    int i = 0;
    // note: a strict greater-than so a tie keeps the value from the earlier
    // buffer, as in the default version
    // note: 2 registers per pass over the buffers halves the loop overhead
    // and gives 2 independent compare chains
    for(; i<len-15; i+=16) // process 2 registers (16 real elements)
    {
      val1 = _mm256_loadu_ps(&in[0][i]);
      val2 = _mm256_loadu_ps(&in[0][i+8]);
      absVal1 = _mm256_and_ps(val1, pos);
      absVal2 = _mm256_and_ps(val2, pos);
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm256_loadu_ps(&in[j][i]);
        // load both registers before using either one, otherwise g++ reads
        // ld from memory a second time
        ld2 = _mm256_loadu_ps(&in[j][i+8]);
        absLd1 = _mm256_and_ps(ld1, pos);
        gt1 = _mm256_cmp_ps(absLd1, absVal1, _CMP_GT_OQ); // |ld| > |val|
        absVal1 = _mm256_max_ps(absLd1, absVal1); // same choice as gt
        // gt? ld : val (avoids _mm256_blendv_ps, see setUnusedToOne)
        val1 = _mm256_or_ps(_mm256_and_ps(gt1, ld1),
            _mm256_andnot_ps(gt1, val1));
        absLd2 = _mm256_and_ps(ld2, pos);
        gt2 = _mm256_cmp_ps(absLd2, absVal2, _CMP_GT_OQ);
        absVal2 = _mm256_max_ps(absLd2, absVal2);
        val2 = _mm256_or_ps(_mm256_and_ps(gt2, ld2),
            _mm256_andnot_ps(gt2, val2));
      }
      _mm256_storeu_ps(&out[i], val1);
      _mm256_storeu_ps(&out[i+8], val2);
    }
    // 0-1 full registers remain and up to 1 partial register
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      val1 = _mm256_loadu_ps(&in[0][i]);
      absVal1 = _mm256_and_ps(val1, pos);
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm256_loadu_ps(&in[j][i]);
        absLd1 = _mm256_and_ps(ld1, pos);
        gt1 = _mm256_cmp_ps(absLd1, absVal1, _CMP_GT_OQ); // |ld| > |val|
        absVal1 = _mm256_max_ps(absLd1, absVal1); // same choice as gt
        // gt1? ld1 : val1 (avoids _mm256_blendv_ps, see setUnusedToOne)
        val1 = _mm256_or_ps(_mm256_and_ps(gt1, ld1),
            _mm256_andnot_ps(gt1, val1));
      }
      _mm256_storeu_ps(&out[i], val1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      val1 = _mm256_maskload_ps(&in[0][i], msk);
      absVal1 = _mm256_and_ps(val1, pos);
      for(int j=1; j<nbufs; ++j)
      {
        ld1 = _mm256_maskload_ps(&in[j][i], msk);
        absLd1 = _mm256_and_ps(ld1, pos);
        gt1 = _mm256_cmp_ps(absLd1, absVal1, _CMP_GT_OQ); // |ld| > |val|
        absVal1 = _mm256_max_ps(absLd1, absVal1); // same choice as gt
        // gt? ld : val (avoids _mm256_blendv_ps, see setUnusedToOne)
        val1 = _mm256_or_ps(_mm256_and_ps(gt1, ld1),
            _mm256_andnot_ps(gt1, val1));
      }
      _mm256_maskstore_ps(&out[i], msk, val1);
    }
    return;
  }
  #endif // AVX signedAbsMax
  __attribute__((__target__("default")))
  inline void signedAbsMax(float const * const * const in, const int nbufs,
      const int len, float * const out)
  {
    if(nbufs<=0 || len<=0) return;
    float val;
    float absVal;
    float tmpVal;
    for(int i=0; i<len; ++i)
    {
      // hold the result in val until every buffer is read in case out is the
      // same buffer as one of the inputs
      val = in[0][i];
      absVal = fabsf(val);
      for(int j=1; j<nbufs; ++j)
      {
        tmpVal = fabsf(in[j][i]);
        if(tmpVal > absVal)
        {
          absVal = tmpVal;
          val = in[j][i];
        }
      }
      out[i] = val;
    }
    return;
  }

  // Replace values of in1 with in2 if in1 < thresh
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void threshCmp(float const * const in1, float const * const in2,
      const int len, const float thresh, float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    __mmask16 ge;
    const __m512 th = _mm512_set1_ps(thresh);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      ge = _mm512_cmp_ps_mask(ld1, th, _CMP_GE_OQ); // in1 >= thresh
      ld1 = _mm512_mask_blend_ps(ge, ld2, ld1); // ge? in1 : in2
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ge = _mm512_cmp_ps_mask(ld1, th, _CMP_GE_OQ); // in1 >= thresh
      ld1 = _mm512_mask_blend_ps(ge, ld2, ld1); // ge? in1 : in2
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 threshCmp
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void threshCmp(float const * const in1, float const * const in2,
      const int len, const float thresh, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ge;
    const __m256 th = _mm256_set1_ps(thresh);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ge = _mm256_cmp_ps(ld1, th, _CMP_GE_OQ); // in1 >= thresh
      // ge? in1 : in2 (avoids _mm256_blendv_ps, see setUnusedToOne)
      ld1 = _mm256_or_ps(_mm256_and_ps(ge, ld1), _mm256_andnot_ps(ge, ld2));
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
      ge = _mm256_cmp_ps(ld1, th, _CMP_GE_OQ); // in1 >= thresh
      // ge? in1 : in2 (avoids _mm256_blendv_ps, see setUnusedToOne)
      ld1 = _mm256_or_ps(_mm256_and_ps(ge, ld1), _mm256_andnot_ps(ge, ld2));
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX threshCmp
  __attribute__((__target__("default")))
  inline void threshCmp(float const * const in1, float const * const in2,
      const int len, const float thresh, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = (in1[i]>=thresh? in1[i]:in2[i]);
    return;
  }

  // Output values of in1 if in3 >= thresh, or in2 if in3 < thresh
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void blend(float const * const in1, float const * const in2,
      float const * const in3, const int len, const float thresh,
      float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3;
    __mmask16 ge;
    const __m512 th = _mm512_set1_ps(thresh);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld3 = _mm512_loadu_ps(&in3[i]);
      ge = _mm512_cmp_ps_mask(ld3, th, _CMP_GE_OQ); // in3 >= thresh
      ld1 = _mm512_mask_blend_ps(ge, ld2, ld1); // ge? in1 : in2
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld3 = _mm512_maskz_loadu_ps(mk, &in3[i]);
      ge = _mm512_cmp_ps_mask(ld3, th, _CMP_GE_OQ); // in3 >= thresh
      ld1 = _mm512_mask_blend_ps(ge, ld2, ld1); // ge? in1 : in2
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 blend
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void blend(float const * const in1, float const * const in2,
      float const * const in3, const int len, const float thresh,
      float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, ge;
    const __m256 th = _mm256_set1_ps(thresh);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld3 = _mm256_loadu_ps(&in3[i]);
      ge = _mm256_cmp_ps(ld3, th, _CMP_GE_OQ); // in3 >= thresh
      // ge? in1 : in2 (avoids _mm256_blendv_ps, see setUnusedToOne)
      ld1 = _mm256_or_ps(_mm256_and_ps(ge, ld1), _mm256_andnot_ps(ge, ld2));
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
      ld3 = _mm256_maskload_ps(&in3[i], msk);
      ge = _mm256_cmp_ps(ld3, th, _CMP_GE_OQ); // in3 >= thresh
      // ge? in1 : in2 (avoids _mm256_blendv_ps, see setUnusedToOne)
      ld1 = _mm256_or_ps(_mm256_and_ps(ge, ld1), _mm256_andnot_ps(ge, ld2));
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX blend
  __attribute__((__target__("default")))
  inline void blend(float const * const in1, float const * const in2,
      float const * const in3, const int len, const float thresh,
      float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = (in3[i]>=thresh? in1[i]:in2[i]);
    return;
  }

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
  
  //===================================================//
  // Vector Multiply and Scale Function Implementation //
  //===================================================//
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

  // real x real into real part of complex
  #if !defined(DISABLE_AVX512) // AVX512 multiply-expand
  __attribute__((__target__("avx512f")))
  inline void multe(float const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    // note: the output is 2x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded, and only the real
    // parts are stored (masked) so the imaginary parts are left as they were
    __m512 ld1, ld2;
    // indices to copy each product to both parts of a complex element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    const __mmask16 re = 0x5555; // real parts of 8 complex elements
    int i = len&~15; // start of the final partial register
    // handle remaining elements first (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld1 = _mm512_mul_ps(ld1, ld2); // [P0,P1,...,P15]
      if(rem>8) // if remainder is > 8, need 2 registers worth
      {
        const __mmask16 mk2 = re & mask16((rem-8)<<1); // 2 floats per complex
        ld2 = _mm512_permutexvar_ps(p2, ld1); // [P8,P8,P9,P9,...,P15,P15]
        ld1 = _mm512_permutexvar_ps(p1, ld1); // [P0,P0,P1,P1,...,P7,P7]
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk2, ld2);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), re, ld1);
      }
      else
      {
        const __mmask16 mk2 = re & mask16(rem<<1); // 2 floats per complex
        ld1 = _mm512_permutexvar_ps(p1, ld1); // [P0,P0,P1,P1,...,P7,P7]
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk2, ld1);
      }
    }
    for(i-=16; i>=0; i-=16) // process 16 real elements (2 output registers)
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld1 = _mm512_mul_ps(ld1, ld2); // [P0,P1,...,P15]
      ld2 = _mm512_permutexvar_ps(p2, ld1); // [P8,P8,P9,P9,...,P15,P15]
      ld1 = _mm512_permutexvar_ps(p1, ld1); // [P0,P0,P1,P1,...,P7,P7]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), re, ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), re, ld1);
    }
    return;
  }
  #endif // end AVX512 multiply-expand
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2")))
  inline void multe(float const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    // note: the output is 2x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded, and only the real
    // parts are stored (masked) so the imaginary parts are left as they were
    __m256 ld1, ld2;
    // indices to copy each product to both parts of a complex element
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256i re = _mm256_setr_epi32(-1,0,-1,0,-1,0,-1,0); // real parts
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld2); // [P0,P1,P2,P3,P4,P5,P6,P7]
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_and_si256(re, _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1])));
        ld2 = _mm256_permutevar8x32_ps(ld1, p2); // [P4,P4,P5,P5,P6,P6,P7,P7]
        ld1 = _mm256_permutevar8x32_ps(ld1, p1); // [P0,P0,P1,P1,P2,P2,P3,P3]
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), re, ld1);
      }
      else
      {
        const __m256i msk2 = _mm256_and_si256(re, _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1])));
        ld1 = _mm256_permutevar8x32_ps(ld1, p1); // [P0,P0,P1,P1,P2,P2,P3,P3]
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld1);
      }
    }
    for(i-=8; i>=0; i-=8) // process 8 real elements (2 output registers)
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_mul_ps(ld1, ld2); // [P0,P1,P2,P3,P4,P5,P6,P7]
      ld2 = _mm256_permutevar8x32_ps(ld1, p2); // [P4,P4,P5,P5,P6,P6,P7,P7]
      ld1 = _mm256_permutevar8x32_ps(ld1, p1); // [P0,P0,P1,P1,P2,P2,P3,P3]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), re, ld2);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), re, ld1);
    }
    return;
  }
  #endif // end AVX2 multiply-expand
  #if !defined(DISABLE_AVX) // AVX multiply-expand
  __attribute__((__target__("avx")))
  inline void multe(float const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    // note: the output is 2x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded, and only the real
    // parts are stored (masked) so the imaginary parts are left as they were
    __m256 ld1, ld2, sh;
    const __m256i re = _mm256_setr_epi32(-1,0,-1,0,-1,0,-1,0); // real parts
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer,
      // integer and needs AVX2, so the masks are and'ed as floats
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld2); // [P0,P1,P2,P3,P4,P5,P6,P7]
      sh = _mm256_unpacklo_ps(ld1, ld1); // [P0,P0,P1,P1,P4,P4,P5,P5]
      ld2 = _mm256_unpackhi_ps(ld1, ld1); // [P2,P2,P3,P3,P6,P6,P7,P7]
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_castps_si256(_mm256_and_ps(
            _mm256_castsi256_ps(re), _mm256_castsi256_ps(_mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1])))));
        ld1 = _mm256_permute2f128_ps(sh, ld2, 0x20); // [P0,P0,...,P3,P3]
        ld2 = _mm256_permute2f128_ps(sh, ld2, 0x31); // [P4,P4,...,P7,P7]
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), re, ld1);
      }
      else
      {
        const __m256i msk2 = _mm256_castps_si256(_mm256_and_ps(
            _mm256_castsi256_ps(re), _mm256_castsi256_ps(_mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1])))));
        ld1 = _mm256_permute2f128_ps(sh, ld2, 0x20); // [P0,P0,...,P3,P3]
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld1);
      }
    }
    for(i-=8; i>=0; i-=8) // process 8 real elements (2 output registers)
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld1 = _mm256_mul_ps(ld1, ld2); // [P0,P1,P2,P3,P4,P5,P6,P7]
      sh = _mm256_unpacklo_ps(ld1, ld1); // [P0,P0,P1,P1,P4,P4,P5,P5]
      ld2 = _mm256_unpackhi_ps(ld1, ld1); // [P2,P2,P3,P3,P6,P6,P7,P7]
      ld1 = _mm256_permute2f128_ps(sh, ld2, 0x20); // [P0,P0,...,P3,P3]
      ld2 = _mm256_permute2f128_ps(sh, ld2, 0x31); // [P4,P4,...,P7,P7]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), re, ld2);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), re, ld1);
    }
    return;
  }
  #endif // end AVX multiply-expand
  // default multiply-expand
  __attribute__((__target__("default")))
  inline void multe(float const * const in1, float const * const in2,
      const int len, complex<float> * const out)
  {
    for(int i=len-1; i>=0; --i) out[i].real(in1[i]*in2[i]);
    return;
  }
  
  
  // complex x real part of complex
  #if !defined(DISABLE_AVX512) // AVX512 remult
  __attribute__((__target__("avx512f")))
  inline void remult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld2 = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld1 = _mm512_mul_ps(ld1, ld2); // [Ar0*Br0,Ai0*Br0,...,Ar7*Br7,Ai7*Br7]
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
      ld2 = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld1 = _mm512_mul_ps(ld1, ld2); // [Ar0*Br0,Ai0*Br0,...,Ar7*Br7,Ai7*Br7]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // end AVX512 remult
  #if !defined(DISABLE_AVX) // AVX
  __attribute__((__target__("avx")))
  inline void remult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld2 = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld1 = _mm256_mul_ps(ld1, ld2); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
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
      ld2 = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld1 = _mm256_mul_ps(ld1, ld2); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX remult
  // default remult
  __attribute__((__target__("default")))
  inline void remult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i].real();
    return;
  }
  
  // complex x imaginary
  #if !defined(DISABLE_AVX512) // AVX512 immult
  __attribute__((__target__("avx512f")))
  inline void immult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld2 = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      ld1 = _mm512_mul_ps(ld1, ld2); // [Ar0*Bi0,Ai0*Bi0,...,Ar7*Bi7,Ai7*Bi7]
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
      ld2 = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      ld1 = _mm512_mul_ps(ld1, ld2); // [Ar0*Bi0,Ai0*Bi0,...,Ar7*Bi7,Ai7*Bi7]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // end AVX512 immult
  #if !defined(DISABLE_AVX) // AVX immult
  __attribute__((__target__("avx")))
  inline void immult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld2 = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      ld1 = _mm256_mul_ps(ld1, ld2); // [Ar0*Bi0,Ai0*Bi0,...,Ar3*Bi3,Ai3*Bi3]
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
      ld2 = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      ld1 = _mm256_mul_ps(ld1, ld2); // [Ar0*Bi0,Ai0*Bi0,...,Ar3*Bi3,Ai3*Bi3]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX immult
  // default immult
  __attribute__((__target__("default")))
  inline void immult(complex<float> const * const in1,
      complex<float> const * const in2, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i].imag();
    return;
  }

  //=============================//
  // Three-Input Vector Multiply //
  //=============================//
  // Multiply and sum, real x real + real
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void multSum(float const * const in1, float const * const in2,
      float const * const in3, const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld3 = _mm512_loadu_ps(&in3[i]);
      ld1 = _mm512_fmadd_ps(ld1, ld2, ld3); // in1*in2 + in3
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld3 = _mm512_maskz_loadu_ps(mk, &in3[i]);
      ld1 = _mm512_fmadd_ps(ld1, ld2, ld3); // in1*in2 + in3
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 real x real + real
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void multSum(float const * const in1, float const * const in2,
      float const * const in3, const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld3 = _mm256_loadu_ps(&in3[i]);
      ld1 = _mm256_fmadd_ps(ld1, ld2, ld3); // in1*in2 + in3
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
      ld3 = _mm256_maskload_ps(&in3[i], msk);
      ld1 = _mm256_fmadd_ps(ld1, ld2, ld3); // in1*in2 + in3
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX2 real x real + real
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void multSum(float const * const in1, float const * const in2,
      float const * const in3, const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(&in2[i]);
      ld3 = _mm256_loadu_ps(&in3[i]);
      ld1 = _mm256_mul_ps(ld1, ld2);
      ld1 = _mm256_add_ps(ld1, ld3); // in1*in2 + in3
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
      ld3 = _mm256_maskload_ps(&in3[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld2);
      ld1 = _mm256_add_ps(ld1, ld3); // in1*in2 + in3
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX real x real + real
  __attribute__((__target__("default")))
  inline void multSum(float const * const in1, float const * const in2,
      float const * const in3, const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i]+in3[i];
    return;
  }
  
  // Multiply and sum, complex x real + complex
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void multSum(complex<float> const * const in1, float const * const in2,
      complex<float> const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1a, ld1b, ld2a, ld2b, ld3a, ld3b;
    // indices to copy each real of in2 to both parts of a complex element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 complex elements (2 registers)
    {
      ld2a = _mm512_loadu_ps(&in2[i]);
      ld1a = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld1b = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      ld3a = _mm512_loadu_ps(reinterpret_cast<float const *>(&in3[i]));
      ld3b = _mm512_loadu_ps(reinterpret_cast<float const *>(&in3[i+8]));
      ld2b = _mm512_permutexvar_ps(p2, ld2a); // [B8,B8,B9,B9,...,B15,B15]
      ld2a = _mm512_permutexvar_ps(p1, ld2a); // [B0,B0,B1,B1,...,B7,B7]
      ld1a = _mm512_fmadd_ps(ld1a, ld2a, ld3a); // A*B + C
      ld1b = _mm512_fmadd_ps(ld1b, ld2b, ld3b);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld1b);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld2a = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld1a = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld1b = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      ld3a = _mm512_loadu_ps(reinterpret_cast<float const *>(&in3[i]));
      ld3b = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in3[i+8]));
      ld2b = _mm512_permutexvar_ps(p2, ld2a); // [B8,B8,B9,B9,...,B15,B15]
      ld2a = _mm512_permutexvar_ps(p1, ld2a); // [B0,B0,B1,B1,...,B7,B7]
      ld1a = _mm512_fmadd_ps(ld1a, ld2a, ld3a); // A*B + C
      ld1b = _mm512_fmadd_ps(ld1b, ld2b, ld3b);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk,
          ld1b);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex
      ld2a = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      ld1a = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      ld3a = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in3[i]));
      ld2a = _mm512_permutexvar_ps(p1, ld2a); // [B0,B0,B1,B1,...,B7,B7]
      ld1a = _mm512_fmadd_ps(ld1a, ld2a, ld3a); // A*B + C
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1a);
    }
    return;
  }
  #endif // AVX512 complex x real + complex
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void multSum(complex<float> const * const in1, float const * const in2,
      complex<float> const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1a, ld1b, ld2a, ld2b, ld3a, ld3b;
    // indices to copy each real of in2 to both parts of a complex element
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld2a = _mm256_loadu_ps(&in2[i]);
      ld1a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld1b = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld3a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i]));
      ld3b = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i+4]));
      ld2b = _mm256_permutevar8x32_ps(ld2a, p2); // [B4,B4,B5,B5,B6,B6,B7,B7]
      ld2a = _mm256_permutevar8x32_ps(ld2a, p1); // [B0,B0,B1,B1,B2,B2,B3,B3]
      ld1a = _mm256_fmadd_ps(ld1a, ld2a, ld3a); // A*B + C
      ld1b = _mm256_fmadd_ps(ld1b, ld2b, ld3b);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld1b);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffers
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld2a = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld1a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld1b = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld3a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i]));
        ld3b = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in3[i+4]), msk2);
        ld2b = _mm256_permutevar8x32_ps(ld2a, p2); // [B4,B4,B5,B5,B6,B6,B7,B7]
        ld2a = _mm256_permutevar8x32_ps(ld2a, p1); // [B0,B0,B1,B1,B2,B2,B3,B3]
        ld1a = _mm256_fmadd_ps(ld1a, ld2a, ld3a); // A*B + C
        ld1b = _mm256_fmadd_ps(ld1b, ld2b, ld3b);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld1b);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld1a = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        ld3a = _mm256_maskload_ps(reinterpret_cast<float const *>(&in3[i]),
            msk2);
        ld2a = _mm256_permutevar8x32_ps(ld2a, p1); // [B0,B0,B1,B1,B2,B2,B3,B3]
        ld1a = _mm256_fmadd_ps(ld1a, ld2a, ld3a); // A*B + C
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld1a);
      }
    }
    return;
  }
  #endif // AVX2 complex x real + complex
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void multSum(complex<float> const * const in1, float const * const in2,
      complex<float> const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1a, ld1b, ld2a, ld2b, ld3a, ld3b, sh;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld2a = _mm256_loadu_ps(&in2[i]);
      ld1a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld1b = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld3a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i]));
      ld3b = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i+4]));
      sh = _mm256_unpacklo_ps(ld2a, ld2a); // [0,0,1,1,4,4,5,5]
      ld2b = _mm256_unpackhi_ps(ld2a, ld2a); // [2,2,3,3,6,6,7,7]
      ld2a = _mm256_permute2f128_ps(sh, ld2b, 0x20); // [0,0,1,1,2,2,3,3]
      ld2b = _mm256_permute2f128_ps(sh, ld2b, 0x31); // [4,4,5,5,6,6,7,7]
      ld1a = _mm256_mul_ps(ld1a, ld2a);
      ld1b = _mm256_mul_ps(ld1b, ld2b);
      ld1a = _mm256_add_ps(ld1a, ld3a); // A*B + C
      ld1b = _mm256_add_ps(ld1b, ld3b);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld1b);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffers
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld2a = _mm256_maskload_ps(&in2[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld1a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld1b = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld3a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i]));
        ld3b = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in3[i+4]), msk2);
        sh = _mm256_unpacklo_ps(ld2a, ld2a); // [0,0,1,1,4,4,5,5]
        ld2b = _mm256_unpackhi_ps(ld2a, ld2a); // [2,2,3,3,6,6,7,7]
        ld2a = _mm256_permute2f128_ps(sh, ld2b, 0x20); // [0,0,1,1,2,2,3,3]
        ld2b = _mm256_permute2f128_ps(sh, ld2b, 0x31); // [4,4,5,5,6,6,7,7]
        ld1a = _mm256_mul_ps(ld1a, ld2a);
        ld1b = _mm256_mul_ps(ld1b, ld2b);
        ld1a = _mm256_add_ps(ld1a, ld3a); // A*B + C
        ld1b = _mm256_add_ps(ld1b, ld3b);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2,
            ld1b);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld1a = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]),
            msk2);
        ld3a = _mm256_maskload_ps(reinterpret_cast<float const *>(&in3[i]),
            msk2);
        sh = _mm256_unpacklo_ps(ld2a, ld2a); // [0,0,1,1,4,4,5,5]
        ld2b = _mm256_unpackhi_ps(ld2a, ld2a); // [2,2,3,3,6,6,7,7]
        ld2a = _mm256_permute2f128_ps(sh, ld2b, 0x20); // [0,0,1,1,2,2,3,3]
        ld1a = _mm256_mul_ps(ld1a, ld2a);
        ld1a = _mm256_add_ps(ld1a, ld3a); // A*B + C
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2,
            ld1a);
      }
    }
    return;
  }
  #endif // AVX complex x real + complex
  __attribute__((__target__("default")))
  inline void multSum(complex<float> const * const in1, float const * const in2,
      complex<float> const * const in3, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i]+in3[i];
    return;
  }
  
  // Multiply and sum, complex x complex + complex
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void multSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sh, re, im;
    int i = 0;
    // note: C is folded into the 2 fused multiplies, the first one gives
    // [Ai*Bi-Cr, Ar*Bi+Ci] and the second Ar*Br minus/plus that
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in3[i]));// C
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld3 = _mm512_fmaddsub_ps(sh, im, ld3); // [Ai0*Bi0-Cr0,Ar0*Bi0+Ci0,...]
      ld1 = _mm512_fmaddsub_ps(ld1, re, ld3);// [Ar0*Br0-Ai0*Bi0+Cr0,...]
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
      ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in3[i])); // C
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld3 = _mm512_fmaddsub_ps(sh, im, ld3); // [Ai0*Bi0-Cr0,Ar0*Bi0+Ci0,...]
      ld1 = _mm512_fmaddsub_ps(ld1, re, ld3);// [Ar0*Br0-Ai0*Bi0+Cr0,...]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x complex + complex
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void multSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sh, re, im;
    int i = 0;
    // note: C is folded into the 2 fused multiplies, the first one gives
    // [Ai*Bi-Cr, Ar*Bi+Ci] and the second Ar*Br minus/plus that
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i]));// C
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld3 = _mm256_fmaddsub_ps(sh, im, ld3); // [Ai0*Bi0-Cr0,Ar0*Bi0+Ci0,...]
      ld1 = _mm256_fmaddsub_ps(ld1, re, ld3);// [Ar0*Br0-Ai0*Bi0+Cr0,...]
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
      ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in3[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld3 = _mm256_fmaddsub_ps(sh, im, ld3); // [Ai0*Bi0-Cr0,Ar0*Bi0+Ci0,...]
      ld1 = _mm256_fmaddsub_ps(ld1, re, ld3);// [Ar0*Br0-Ai0*Bi0+Cr0,...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX2 complex x complex + complex
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void multSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sh, re, im;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i]));// C
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      ld1 = _mm256_add_ps(ld1, ld3); // A*B + C
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
      ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in3[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      ld1 = _mm256_add_ps(ld1, ld3); // A*B + C
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX complex x complex + complex
  __attribute__((__target__("default")))
  inline void multSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i)
    {
      const float re = in1[i].real()*in2[i].real() -
          in1[i].imag()*in2[i].imag() + in3[i].real();
      const float im = in1[i].imag()*in2[i].real() +
          in1[i].real()*in2[i].imag() + in3[i].imag();
      out[i].real(re);
      out[i].imag(im);
    }
    return;
  }

  // Multiply and sum, complex x conj(complex) + complex
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  inline void multcSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sh, re, im;
    int i = 0;
    // note: C is folded into the 2 fused multiplies, the first one gives
    // [Ai*Bi+Cr, Ar*Bi-Ci] and the second Ar*Br plus/minus that
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in3[i]));// C
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld3 = _mm512_fmsubadd_ps(sh, im, ld3); // [Ai0*Bi0+Cr0,Ar0*Bi0-Ci0,...]
      ld1 = _mm512_fmsubadd_ps(ld1, re, ld3);// [Ar0*Br0+Ai0*Bi0+Cr0,...]
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
      ld3 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(
          &in3[i])); // C
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld3 = _mm512_fmsubadd_ps(sh, im, ld3); // [Ai0*Bi0+Cr0,Ar0*Bi0-Ci0,...]
      ld1 = _mm512_fmsubadd_ps(ld1, re, ld3);// [Ar0*Br0+Ai0*Bi0+Cr0,...]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x conj(complex) + complex
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  inline void multcSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sh, re, im;
    int i = 0;
    // note: C is folded into the 2 fused multiplies, the first one gives
    // [Ai*Bi+Cr, Ar*Bi-Ci] and the second Ar*Br plus/minus that
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i]));// C
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld3 = _mm256_fmsubadd_ps(sh, im, ld3); // [Ai0*Bi0+Cr0,Ar0*Bi0-Ci0,...]
      ld1 = _mm256_fmsubadd_ps(ld1, re, ld3);// [Ar0*Br0+Ai0*Bi0+Cr0,...]
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
      ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in3[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld3 = _mm256_fmsubadd_ps(sh, im, ld3); // [Ai0*Bi0+Cr0,Ar0*Bi0-Ci0,...]
      ld1 = _mm256_fmsubadd_ps(ld1, re, ld3);// [Ar0*Br0+Ai0*Bi0+Cr0,...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX2 complex x conj(complex) + complex
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  inline void multcSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sh, re, im;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in3[i]));// C
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // -[Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // -[Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0+Ai0*Bi0,Ai0*Br0-Ar0*Bi0]
      ld1 = _mm256_add_ps(ld1, ld3); // A*conj(B) + C
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
      ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in3[i]),
          msk);
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // -[Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // -[Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0+Ai0*Bi0,Ai0*Br0-Ar0*Bi0]
      ld1 = _mm256_add_ps(ld1, ld3); // A*conj(B) + C
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX complex x conj(complex) + complex
  __attribute__((__target__("default")))
  inline void multcSum(complex<float> const * const in1,
      complex<float> const * const in2, complex<float> const * const in3,
      const int len, complex<float> * const out)
  {
    // written out rather than in1[i]*conj(in2[i]), see the icpx note above
    for(int i=0; i<len; ++i)
    {
      const float re = in1[i].real()*in2[i].real() +
          in1[i].imag()*in2[i].imag() + in3[i].real();
      const float im = in1[i].imag()*in2[i].real() -
          in1[i].real()*in2[i].imag() + in3[i].imag();
      out[i].real(re);
      out[i].imag(im);
    }
    return;
  }

  // complex x complex x real
  #if !defined(DISABLE_AVX512) // AVX512 multMult
  __attribute__((__target__("avx512f")))
  inline void multMult(complex<float> const * const in1,
      complex<float> const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1a, ld1b, ld2a, ld2b, ld3a, ld3b, sh, re, im;
    // indices to copy each real to both parts of a complex element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 complex elements (2 registers)
    {
      ld3a = _mm512_loadu_ps(&in3[i]);
      ld1a = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));   // A
      ld1b = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      ld2a = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));   // B
      ld2b = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
      ld3b = _mm512_permutexvar_ps(p2, ld3a); // [C8,C8,C9,C9,...,C15,C15]
      ld3a = _mm512_permutexvar_ps(p1, ld3a); // [C0,C0,C1,C1,...,C7,C7]
      sh = _mm512_shuffle_ps(ld1a, ld1a, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2a); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2a); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2a = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1a = _mm512_fmaddsub_ps(re, ld1a, ld2a); // A*B
      ld1a = _mm512_mul_ps(ld1a, ld3a); // A*B*C
      sh = _mm512_shuffle_ps(ld1b, ld1b, 0xb1);
      im = _mm512_movehdup_ps(ld2b);
      re = _mm512_moveldup_ps(ld2b);
      ld2b = _mm512_mul_ps(sh, im);
      ld1b = _mm512_fmaddsub_ps(re, ld1b, ld2b);
      ld1b = _mm512_mul_ps(ld1b, ld3b);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld1b);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld3a = _mm512_maskz_loadu_ps(mask16(rem), &in3[i]);
      ld1a = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld1b = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      ld2a = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld2b = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i+8]));
      ld3b = _mm512_permutexvar_ps(p2, ld3a); // [C8,C8,C9,C9,...,C15,C15]
      ld3a = _mm512_permutexvar_ps(p1, ld3a); // [C0,C0,C1,C1,...,C7,C7]
      sh = _mm512_shuffle_ps(ld1a, ld1a, 0xb1);
      im = _mm512_movehdup_ps(ld2a);
      re = _mm512_moveldup_ps(ld2a);
      ld2a = _mm512_mul_ps(sh, im);
      ld1a = _mm512_fmaddsub_ps(re, ld1a, ld2a);
      ld1a = _mm512_mul_ps(ld1a, ld3a);
      sh = _mm512_shuffle_ps(ld1b, ld1b, 0xb1);
      im = _mm512_movehdup_ps(ld2b);
      re = _mm512_moveldup_ps(ld2b);
      ld2b = _mm512_mul_ps(sh, im);
      ld1b = _mm512_fmaddsub_ps(re, ld1b, ld2b);
      ld1b = _mm512_mul_ps(ld1b, ld3b);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk, ld1b);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex
      ld3a = _mm512_maskz_loadu_ps(mask16(rem), &in3[i]);
      ld1a = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      ld2a = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i]));
      ld3a = _mm512_permutexvar_ps(p1, ld3a); // [C0,C0,C1,C1,...,C7,C7]
      sh = _mm512_shuffle_ps(ld1a, ld1a, 0xb1);
      im = _mm512_movehdup_ps(ld2a);
      re = _mm512_moveldup_ps(ld2a);
      ld2a = _mm512_mul_ps(sh, im);
      ld1a = _mm512_fmaddsub_ps(re, ld1a, ld2a);
      ld1a = _mm512_mul_ps(ld1a, ld3a);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1a);
    }
    return;
  }
  #endif // end AVX512 multMult
  #if !defined(DISABLE_AVX2) // AVX2 multMult
  __attribute__((__target__("avx2,fma")))
  inline void multMult(complex<float> const * const in1,
      complex<float> const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1a, ld1b, ld2a, ld2b, ld3a, ld3b, sh, re, im;
    // indices to copy each real to both parts of a complex element
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld3a = _mm256_loadu_ps(&in3[i]);
      ld1a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));   // A
      ld1b = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld2a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));   // B
      ld2b = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      ld3b = _mm256_permutevar8x32_ps(ld3a, p2); // [C4,C4,C5,C5,C6,C6,C7,C7]
      ld3a = _mm256_permutevar8x32_ps(ld3a, p1); // [C0,C0,C1,C1,C2,C2,C3,C3]
      sh = _mm256_shuffle_ps(ld1a, ld1a, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2a); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2a); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2a = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1a = _mm256_fmaddsub_ps(re, ld1a, ld2a); // A*B
      ld1a = _mm256_mul_ps(ld1a, ld3a); // A*B*C
      sh = _mm256_shuffle_ps(ld1b, ld1b, 0xb1);
      im = _mm256_movehdup_ps(ld2b);
      re = _mm256_moveldup_ps(ld2b);
      ld2b = _mm256_mul_ps(sh, im);
      ld1b = _mm256_fmaddsub_ps(re, ld1b, ld2b);
      ld1b = _mm256_mul_ps(ld1b, ld3b);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld1b);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffers
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld3a = _mm256_maskload_ps(&in3[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld1a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld1b = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld2a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld2b = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        ld3b = _mm256_permutevar8x32_ps(ld3a, p2); // [C4,C4,...,C7,C7]
        ld3a = _mm256_permutevar8x32_ps(ld3a, p1); // [C0,C0,...,C3,C3]
        sh = _mm256_shuffle_ps(ld1a, ld1a, 0xb1);
        im = _mm256_movehdup_ps(ld2a);
        re = _mm256_moveldup_ps(ld2a);
        ld2a = _mm256_mul_ps(sh, im);
        ld1a = _mm256_fmaddsub_ps(re, ld1a, ld2a);
        ld1a = _mm256_mul_ps(ld1a, ld3a);
        sh = _mm256_shuffle_ps(ld1b, ld1b, 0xb1);
        im = _mm256_movehdup_ps(ld2b);
        re = _mm256_moveldup_ps(ld2b);
        ld2b = _mm256_mul_ps(sh, im);
        ld1b = _mm256_fmaddsub_ps(re, ld1b, ld2b);
        ld1b = _mm256_mul_ps(ld1b, ld3b);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld1b);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld1a = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i]), msk2);
        ld2a = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i]), msk2);
        ld3a = _mm256_permutevar8x32_ps(ld3a, p1); // [C0,C0,...,C3,C3]
        sh = _mm256_shuffle_ps(ld1a, ld1a, 0xb1);
        im = _mm256_movehdup_ps(ld2a);
        re = _mm256_moveldup_ps(ld2a);
        ld2a = _mm256_mul_ps(sh, im);
        ld1a = _mm256_fmaddsub_ps(re, ld1a, ld2a);
        ld1a = _mm256_mul_ps(ld1a, ld3a);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld1a);
      }
    }
    return;
  }
  #endif // end AVX2 multMult
  #if !defined(DISABLE_AVX) // AVX multMult
  __attribute__((__target__("avx")))
  inline void multMult(complex<float> const * const in1,
      complex<float> const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1a, ld1b, ld2a, ld2b, ld3a, ld3b, sh, re, im;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld3a = _mm256_loadu_ps(&in3[i]);
      ld1a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));   // A
      ld1b = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      ld2a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));   // B
      ld2b = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sh = _mm256_unpacklo_ps(ld3a, ld3a);   // [C0,C0,C1,C1,C4,C4,C5,C5]
      ld3b = _mm256_unpackhi_ps(ld3a, ld3a); // [C2,C2,C3,C3,C6,C6,C7,C7]
      ld3a = _mm256_permute2f128_ps(sh, ld3b, 0x20); // [C0,C0,...,C3,C3]
      ld3b = _mm256_permute2f128_ps(sh, ld3b, 0x31); // [C4,C4,...,C7,C7]
      sh = _mm256_shuffle_ps(ld1a, ld1a, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2a); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2a); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2a = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1a = _mm256_mul_ps(ld1a, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1a = _mm256_addsub_ps(ld1a, ld2a); // A*B
      ld1a = _mm256_mul_ps(ld1a, ld3a); // A*B*C
      sh = _mm256_shuffle_ps(ld1b, ld1b, 0xb1);
      im = _mm256_movehdup_ps(ld2b);
      re = _mm256_moveldup_ps(ld2b);
      ld2b = _mm256_mul_ps(sh, im);
      ld1b = _mm256_mul_ps(ld1b, re);
      ld1b = _mm256_addsub_ps(ld1b, ld2b);
      ld1b = _mm256_mul_ps(ld1b, ld3b);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld1b);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffers
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld3a = _mm256_maskload_ps(&in3[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld1a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld1b = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        ld2a = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld2b = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i+4]), msk2);
        sh = _mm256_unpacklo_ps(ld3a, ld3a);   // [C0,C0,C1,C1,C4,C4,C5,C5]
        ld3b = _mm256_unpackhi_ps(ld3a, ld3a); // [C2,C2,C3,C3,C6,C6,C7,C7]
        ld3a = _mm256_permute2f128_ps(sh, ld3b, 0x20); // [C0,C0,...,C3,C3]
        ld3b = _mm256_permute2f128_ps(sh, ld3b, 0x31); // [C4,C4,...,C7,C7]
        sh = _mm256_shuffle_ps(ld1a, ld1a, 0xb1);
        im = _mm256_movehdup_ps(ld2a);
        re = _mm256_moveldup_ps(ld2a);
        ld2a = _mm256_mul_ps(sh, im);
        ld1a = _mm256_mul_ps(ld1a, re);
        ld1a = _mm256_addsub_ps(ld1a, ld2a);
        ld1a = _mm256_mul_ps(ld1a, ld3a);
        sh = _mm256_shuffle_ps(ld1b, ld1b, 0xb1);
        im = _mm256_movehdup_ps(ld2b);
        re = _mm256_moveldup_ps(ld2b);
        ld2b = _mm256_mul_ps(sh, im);
        ld1b = _mm256_mul_ps(ld1b, re);
        ld1b = _mm256_addsub_ps(ld1b, ld2b);
        ld1b = _mm256_mul_ps(ld1b, ld3b);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1a);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld1b);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld1a = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i]), msk2);
        ld2a = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in2[i]), msk2);
        sh = _mm256_unpacklo_ps(ld3a, ld3a);   // [C0,C0,C1,C1,C4,C4,C5,C5]
        ld3b = _mm256_unpackhi_ps(ld3a, ld3a); // [C2,C2,C3,C3,C6,C6,C7,C7]
        ld3a = _mm256_permute2f128_ps(sh, ld3b, 0x20); // [C0,C0,...,C3,C3]
        sh = _mm256_shuffle_ps(ld1a, ld1a, 0xb1);
        im = _mm256_movehdup_ps(ld2a);
        re = _mm256_moveldup_ps(ld2a);
        ld2a = _mm256_mul_ps(sh, im);
        ld1a = _mm256_mul_ps(ld1a, re);
        ld1a = _mm256_addsub_ps(ld1a, ld2a);
        ld1a = _mm256_mul_ps(ld1a, ld3a);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld1a);
      }
    }
    return;
  }
  #endif // end AVX multMult
  // default multMult
  __attribute__((__target__("default")))
  inline void multMult(complex<float> const * const in1,
      complex<float> const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i)
    {
      const float re = (in1[i].real()*in2[i].real() -
          in1[i].imag()*in2[i].imag())*in3[i];
      const float im = (in1[i].imag()*in2[i].real() +
          in1[i].real()*in2[i].imag())*in3[i];
      out[i].real(re);
      out[i].imag(im);
    }
    return;
  }
  
  // complex x real / real with true divide
  #if !defined(DISABLE_AVX512) // AVX512 multDiv
  __attribute__((__target__("avx512f")))
  inline void multDiv(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    // indices to copy each real to both parts of a complex element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 complex elements (2 registers)
    {
      ld1 = _mm512_loadu_ps(&in2[i]);
      sc = _mm512_loadu_ps(&in3[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      // compute in2/in3 with a single divide, then spread it to each complex
      ld1 = _mm512_div_ps(ld1, sc); // in2/in3
      sc = _mm512_permutexvar_ps(p1, ld1);  // [Q0,Q0,Q1,Q1,...,Q7,Q7]
      ld1 = _mm512_permutexvar_ps(p2, ld1); // [Q8,Q8,Q9,Q9,...,Q15,Q15]
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16), the unused lanes of
    // the divisor are masked out so they don't divide 0 by 0
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      sc = _mm512_maskz_loadu_ps(mask16(rem), &in3[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      ld1 = _mm512_maskz_div_ps(mask16(rem), ld1, sc);
      sc = _mm512_permutexvar_ps(p1, ld1);  // [Q0,Q0,Q1,Q1,...,Q7,Q7]
      ld1 = _mm512_permutexvar_ps(p2, ld1); // [Q8,Q8,Q9,Q9,...,Q15,Q15]
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk, ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      sc = _mm512_maskz_loadu_ps(mask16(rem), &in3[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      ld1 = _mm512_maskz_div_ps(mask16(rem), ld1, sc);
      sc = _mm512_permutexvar_ps(p1, ld1);  // [Q0,Q0,Q1,Q1,...,Q7,Q7]
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 multDiv
  #if !defined(DISABLE_AVX2) // AVX2 multDiv
  __attribute__((__target__("avx2")))
  inline void multDiv(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    // indices to copy each real to both parts of a complex element
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      sc = _mm256_loadu_ps(&in3[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // compute in2/in3 with a single divide, then spread it to each complex
      ld1 = _mm256_div_ps(ld1, sc); // in2/in3
      sc = _mm256_permutevar8x32_ps(ld1, p1);  // [Q0,Q0,Q1,Q1,Q2,Q2,Q3,Q3]
      ld1 = _mm256_permutevar8x32_ps(ld1, p2); // [Q4,Q4,Q5,Q5,Q6,Q6,Q7,Q7]
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
      sc = _mm256_maskload_ps(&in3[i], msk1);
      ld1 = _mm256_div_ps(ld1, sc); // in2/in3
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);  // [Q0,Q0,Q1,Q1,Q2,Q2,Q3,Q3]
        ld1 = _mm256_permutevar8x32_ps(ld1, p2); // [Q4,Q4,Q5,Q5,Q6,Q6,Q7,Q7]
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);  // [Q0,Q0,Q1,Q1,Q2,Q2,Q3,Q3]
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
      }
    }
    return;
  }
  #endif // end AVX2 
  #if !defined(DISABLE_AVX) // AVX multDiv
  __attribute__((__target__("avx")))
  inline void multDiv(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc, tmp;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      sc = _mm256_loadu_ps(&in3[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // compute in2/in3 with a single divide, then spread it to each complex
      ld1 = _mm256_div_ps(ld1, sc); // in2/in3
      sc = _mm256_unpacklo_ps(ld1, ld1);  // [Q0,Q0,Q1,Q1,Q4,Q4,Q5,Q5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [Q2,Q2,Q3,Q3,Q6,Q6,Q7,Q7]
      tmp = _mm256_permute2f128_ps(sc, ld1, 0x20); // [Q0,Q0,...,Q3,Q3]
      ld1 = _mm256_permute2f128_ps(sc, ld1, 0x31); // [Q4,Q4,...,Q7,Q7]
      sc = tmp;
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
      sc = _mm256_maskload_ps(&in3[i], msk1);
      ld1 = _mm256_div_ps(ld1, sc); // in2/in3
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        sc = _mm256_unpacklo_ps(ld1, ld1);  // [Q0,Q0,Q1,Q1,Q4,Q4,Q5,Q5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [Q2,Q2,Q3,Q3,Q6,Q6,Q7,Q7]
        tmp = _mm256_permute2f128_ps(sc, ld1, 0x20); // [Q0,Q0,...,Q3,Q3]
        ld1 = _mm256_permute2f128_ps(sc, ld1, 0x31); // [Q4,Q4,...,Q7,Q7]
        sc = tmp;
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i]), msk2);
        sc = _mm256_unpacklo_ps(ld1, ld1);  // [Q0,Q0,Q1,Q1,Q4,Q4,Q5,Q5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [Q2,Q2,Q3,Q3,Q6,Q6,Q7,Q7]
        sc = _mm256_permute2f128_ps(sc, ld1, 0x20); // [Q0,Q0,...,Q3,Q3]
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
      }
    }
    return;
  }
  #endif // end AVX multDiv
  // default multDiv
  __attribute__((__target__("default")))
  inline void multDiv(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i)
    {
      const float sc = in2[i]/in3[i];
      out[i].real(in1[i].real()*sc);
      out[i].imag(in1[i].imag()*sc);
    }
    return;
  }
  
  // complex x real / real with reciprocal approx.
  #if !defined(DISABLE_AVX512) // AVX512 multDivr
  __attribute__((__target__("avx512f")))
  inline void multDivr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc;
    // indices to copy each real to both parts of a complex element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 complex elements (2 registers)
    {
      ld1 = _mm512_loadu_ps(&in2[i]);
      sc = _mm512_loadu_ps(&in3[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      // in2/in3 from a reciprocal approximation, spread to each complex
      ld1 = _mm512_mul_ps(ld1, _mm512_rcp14_ps(sc)); // in2*(~1/in3)
      sc = _mm512_permutexvar_ps(p1, ld1);  // [Q0,Q0,Q1,Q1,...,Q7,Q7]
      ld1 = _mm512_permutexvar_ps(p2, ld1); // [Q8,Q8,Q9,Q9,...,Q15,Q15]
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16), the unused lanes of
    // the divisor are masked out so they don't divide 0 by 0
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      sc = _mm512_maskz_loadu_ps(mask16(rem), &in3[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      ld1 = _mm512_mul_ps(ld1, _mm512_maskz_rcp14_ps(mask16(rem), sc));
      sc = _mm512_permutexvar_ps(p1, ld1);  // [Q0,Q0,Q1,Q1,...,Q7,Q7]
      ld1 = _mm512_permutexvar_ps(p2, ld1); // [Q8,Q8,Q9,Q9,...,Q15,Q15]
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk, ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      sc = _mm512_maskz_loadu_ps(mask16(rem), &in3[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      ld1 = _mm512_mul_ps(ld1, _mm512_maskz_rcp14_ps(mask16(rem), sc));
      sc = _mm512_permutexvar_ps(p1, ld1);  // [Q0,Q0,Q1,Q1,...,Q7,Q7]
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 multDivr
  #if !defined(DISABLE_AVX2) // AVX2 multDivr
  __attribute__((__target__("avx2")))
  inline void multDivr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc;
    // indices to copy each real to both parts of a complex element
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      sc = _mm256_loadu_ps(&in3[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // in2/in3 from a reciprocal approximation, spread to each complex
      ld1 = _mm256_mul_ps(ld1, _mm256_rcp_ps(sc)); // in2*(~1/in3)
      sc = _mm256_permutevar8x32_ps(ld1, p1);  // [Q0,Q0,Q1,Q1,Q2,Q2,Q3,Q3]
      ld1 = _mm256_permutevar8x32_ps(ld1, p2); // [Q4,Q4,Q5,Q5,Q6,Q6,Q7,Q7]
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
      sc = _mm256_maskload_ps(&in3[i], msk1);
      ld1 = _mm256_mul_ps(ld1, _mm256_rcp_ps(sc)); // in2*(~1/in3)
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);  // [Q0,Q0,Q1,Q1,Q2,Q2,Q3,Q3]
        ld1 = _mm256_permutevar8x32_ps(ld1, p2); // [Q4,Q4,Q5,Q5,Q6,Q6,Q7,Q7]
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);  // [Q0,Q0,Q1,Q1,Q2,Q2,Q3,Q3]
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
      }
    }
    return;
  }
  #endif // end AVX2 multDivr
  #if !defined(DISABLE_AVX) // AVX multDivr
  __attribute__((__target__("avx")))
  inline void multDivr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc, tmp;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      sc = _mm256_loadu_ps(&in3[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // in2/in3 from a reciprocal approximation, spread to each complex
      ld1 = _mm256_mul_ps(ld1, _mm256_rcp_ps(sc)); // in2*(~1/in3)
      sc = _mm256_unpacklo_ps(ld1, ld1);  // [Q0,Q0,Q1,Q1,Q4,Q4,Q5,Q5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [Q2,Q2,Q3,Q3,Q6,Q6,Q7,Q7]
      tmp = _mm256_permute2f128_ps(sc, ld1, 0x20); // [Q0,Q0,...,Q3,Q3]
      ld1 = _mm256_permute2f128_ps(sc, ld1, 0x31); // [Q4,Q4,...,Q7,Q7]
      sc = tmp;
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
      sc = _mm256_maskload_ps(&in3[i], msk1);
      ld1 = _mm256_mul_ps(ld1, _mm256_rcp_ps(sc)); // in2*(~1/in3)
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        sc = _mm256_unpacklo_ps(ld1, ld1);  // [Q0,Q0,Q1,Q1,Q4,Q4,Q5,Q5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [Q2,Q2,Q3,Q3,Q6,Q6,Q7,Q7]
        tmp = _mm256_permute2f128_ps(sc, ld1, 0x20); // [Q0,Q0,...,Q3,Q3]
        ld1 = _mm256_permute2f128_ps(sc, ld1, 0x31); // [Q4,Q4,...,Q7,Q7]
        sc = tmp;
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i]), msk2);
        sc = _mm256_unpacklo_ps(ld1, ld1);  // [Q0,Q0,Q1,Q1,Q4,Q4,Q5,Q5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [Q2,Q2,Q3,Q3,Q6,Q6,Q7,Q7]
        sc = _mm256_permute2f128_ps(sc, ld1, 0x20); // [Q0,Q0,...,Q3,Q3]
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
      }
    }
    return;
  }
  #endif // end AVX multDivr
  // default multDivr
  __attribute__((__target__("default")))
  inline void multDivr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i)
    {
      const float sc = in2[i]/in3[i];
      out[i].real(in1[i].real()*sc);
      out[i].imag(in1[i].imag()*sc);
    }
    return;
  }

  // complex x real / real with reciprocal approx. and 1 stage Newton-Raphson
  // refinement
  #if !defined(DISABLE_AVX512) // AVX512 multDivnr
  __attribute__((__target__("avx512f")))
  inline void multDivnr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, ld3, sc, inv, q;
    // indices to copy each real to both parts of a complex element
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 complex elements (2 registers)
    {
      ld1 = _mm512_loadu_ps(&in2[i]);
      sc = _mm512_loadu_ps(&in3[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
      // in2/in3 from a reciprocal approximation r and one newton-raphson
      // stage applied to the quotient rather than to r (q = q0+r*(in2-in3*q0)
      // with q0 = in2*r), with FMA this costs the same and is more accurate,
      // then spread it to each complex
      inv = _mm512_rcp14_ps(sc);           // r
      q = _mm512_mul_ps(ld1, inv);         // q0 = in2*r
      ld1 = _mm512_fnmadd_ps(sc, q, ld1);  // in2-in3*q0
      ld1 = _mm512_fmadd_ps(ld1, inv, q);  // q0+r*(in2-in3*q0)
      sc = _mm512_permutexvar_ps(p1, ld1);  // [Q0,Q0,Q1,Q1,...,Q7,Q7]
      ld1 = _mm512_permutexvar_ps(p2, ld1); // [Q8,Q8,Q9,Q9,...,Q15,Q15]
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16), the unused lanes of
    // the divisor are masked out so they don't divide 0 by 0
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      sc = _mm512_maskz_loadu_ps(mask16(rem), &in3[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i+8]));
      inv = _mm512_maskz_rcp14_ps(mask16(rem), sc);
      q = _mm512_mul_ps(ld1, inv);
      ld1 = _mm512_fnmadd_ps(sc, q, ld1);
      ld1 = _mm512_fmadd_ps(ld1, inv, q);
      sc = _mm512_permutexvar_ps(p1, ld1);  // [Q0,Q0,Q1,Q1,...,Q7,Q7]
      ld1 = _mm512_permutexvar_ps(p2, ld1); // [Q8,Q8,Q9,Q9,...,Q15,Q15]
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk, ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in2[i]);
      sc = _mm512_maskz_loadu_ps(mask16(rem), &in3[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in1[i]));
      inv = _mm512_maskz_rcp14_ps(mask16(rem), sc);
      q = _mm512_mul_ps(ld1, inv);
      ld1 = _mm512_fnmadd_ps(sc, q, ld1);
      ld1 = _mm512_fmadd_ps(ld1, inv, q);
      sc = _mm512_permutexvar_ps(p1, ld1);  // [Q0,Q0,Q1,Q1,...,Q7,Q7]
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 multDivnr
  #if !defined(DISABLE_AVX2) // AVX2 multDivnr
  __attribute__((__target__("avx2,fma")))
  inline void multDivnr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc, inv, q;
    // indices to copy each real to both parts of a complex element
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      sc = _mm256_loadu_ps(&in3[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // in2/in3 from a reciprocal approximation r and one newton-raphson
      // stage applied to the quotient rather than to r (q = q0+r*(in2-in3*q0)
      // with q0 = in2*r), with FMA this costs the same and is more accurate,
      // then spread it to each complex
      inv = _mm256_rcp_ps(sc);             // r
      q = _mm256_mul_ps(ld1, inv);         // q0 = in2*r
      ld1 = _mm256_fnmadd_ps(sc, q, ld1);  // in2-in3*q0
      ld1 = _mm256_fmadd_ps(ld1, inv, q);  // q0+r*(in2-in3*q0)
      sc = _mm256_permutevar8x32_ps(ld1, p1);  // [Q0,Q0,Q1,Q1,Q2,Q2,Q3,Q3]
      ld1 = _mm256_permutevar8x32_ps(ld1, p2); // [Q4,Q4,Q5,Q5,Q6,Q6,Q7,Q7]
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
      sc = _mm256_maskload_ps(&in3[i], msk1);
      inv = _mm256_rcp_ps(sc);
      q = _mm256_mul_ps(ld1, inv);
      ld1 = _mm256_fnmadd_ps(sc, q, ld1);
      ld1 = _mm256_fmadd_ps(ld1, inv, q);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);  // [Q0,Q0,Q1,Q1,Q2,Q2,Q3,Q3]
        ld1 = _mm256_permutevar8x32_ps(ld1, p2); // [Q4,Q4,Q5,Q5,Q6,Q6,Q7,Q7]
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);  // [Q0,Q0,Q1,Q1,Q2,Q2,Q3,Q3]
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
      }
    }
    return;
  }
  #endif // end AVX2 multDivnr
  #if !defined(DISABLE_AVX) // AVX multDivnr
  __attribute__((__target__("avx")))
  inline void multDivnr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, ld3, sc, tmp, inv;
    // constant for newton-raphson refinement of reciprocal approximation
    const __m256 two = _mm256_set1_ps(2.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(&in2[i]);
      sc = _mm256_loadu_ps(&in3[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
      // in2/in3 from a reciprocal approximation and one newton-raphson
      // stage (new = old*(2-orig*old)), then spread it to each complex
      inv = _mm256_rcp_ps(sc);     // old
      sc = _mm256_mul_ps(sc, inv); // orig*old
      sc = _mm256_sub_ps(two, sc); // 2-orig*old
      sc = _mm256_mul_ps(sc, inv); // old*(2-orig*old)
      ld1 = _mm256_mul_ps(ld1, sc);
      sc = _mm256_unpacklo_ps(ld1, ld1);  // [Q0,Q0,Q1,Q1,Q4,Q4,Q5,Q5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [Q2,Q2,Q3,Q3,Q6,Q6,Q7,Q7]
      tmp = _mm256_permute2f128_ps(sc, ld1, 0x20); // [Q0,Q0,...,Q3,Q3]
      ld1 = _mm256_permute2f128_ps(sc, ld1, 0x31); // [Q4,Q4,...,Q7,Q7]
      sc = tmp;
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
      sc = _mm256_maskload_ps(&in3[i], msk1);
      inv = _mm256_rcp_ps(sc);
      sc = _mm256_mul_ps(sc, inv);
      sc = _mm256_sub_ps(two, sc);
      sc = _mm256_mul_ps(sc, inv);
      ld1 = _mm256_mul_ps(ld1, sc);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i+4]), msk2);
        sc = _mm256_unpacklo_ps(ld1, ld1);  // [Q0,Q0,Q1,Q1,Q4,Q4,Q5,Q5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [Q2,Q2,Q3,Q3,Q6,Q6,Q7,Q7]
        tmp = _mm256_permute2f128_ps(sc, ld1, 0x20); // [Q0,Q0,...,Q3,Q3]
        ld1 = _mm256_permute2f128_ps(sc, ld1, 0x31); // [Q4,Q4,...,Q7,Q7]
        sc = tmp;
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk2, ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(
            &in1[i]), msk2);
        sc = _mm256_unpacklo_ps(ld1, ld1);  // [Q0,Q0,Q1,Q1,Q4,Q4,Q5,Q5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [Q2,Q2,Q3,Q3,Q6,Q6,Q7,Q7]
        sc = _mm256_permute2f128_ps(sc, ld1, 0x20); // [Q0,Q0,...,Q3,Q3]
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk2, ld2);
      }
    }
    return;
  }
  #endif // end AVX multDivnr
  // default multDivnr
  __attribute__((__target__("default")))
  inline void multDivnr(complex<float> const * const in1,
      float const * const in2, float const * const in3, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i)
    {
      const float sc = in2[i]/in3[i];
      out[i].real(in1[i].real()*sc);
      out[i].imag(in1[i].imag()*sc);
    }
    return;
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
  inline void dotp(float const * const in1, float const * const in2,
      const int len, float &out)
  {
    if(len<=0)
    {
      out = 0.0f;
      return;
    }
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
    out = _mm512_reduce_add_ps(sum);
    return;
  }
  #endif // end AVX512 real*real dotp
  #if !defined(DISABLE_AVX2) // AVX2 real*real dotp
  __attribute__((__target__("avx2,fma")))
  inline void dotp(float const * const in1, float const * const in2,
      const int len, float &out)
  {
    if(len<=0)
    {
      out = 0.0f;
      return;
    }
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
    out = _mm_cvtss_f32(hs); // extract the first element in the register
    return;
  }
  #endif // end AVX2 real*real dotp
  #if !defined(DISABLE_AVX) // AVX real*real dotp
  __attribute__((__target__("avx")))
  inline void dotp(float const * const in1, float const * const in2,
      const int len, float &out)
  {
    if(len<=0)
    {
      out = 0.0f;
      return;
    }
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
    out = _mm_cvtss_f32(hs); // extract the first element in the register
    return;
  }
  #endif // end AVX real*real dotp
  __attribute__((__target__("default"))) // default real*real dotp
  inline void dotp(float const * const in1, float const * const in2,
      const int len, float &out)
  {
    // use a separate sum so the compiler knows that it won't alias with the
    // inputs (when summing directly into out the compiler can not assume this
    // and it costs load/store cycles)
    float sum = 0.0f;
    for(int i=0; i<len; ++i) sum += (in1[i]*in2[i]);
    out = sum;
    return;
  }
  // wrapper returning result by value, must be defined after every target
  // version so the calls go through the runtime dispatcher
  inline float dotp(float const * const in1, float const * const in2,
      const int len)
  {
    float out;
    dotp(in1, in2, len, out);
    return out;
  }

  // complex*complex dot-product
  #if !defined(DISABLE_AVX512) // AVX512 complex*complex dotp
  __attribute__((__target__("avx512f")))
  inline void dotp(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m512 ld1, ld2, sh;
    __m512 sumr = _mm512_set1_ps(0.0f); // real part of dot-product
    __m512 sumi = sumr;                 // imaginary part of dot-product
    __m256 hs2;  // horizontal sum, 8 elements (fold and sum halves of m512)
    __m128 hs1;  // horizontal sum to finish dot-product (fold and sum m256)
    int i = 0;
    // note: with a single pair of accumulators every FMA waits on the previous
    // one, so longer buffers use 2 pairs to keep the FMA unit busy, combining
    // them only if the 2 register loop runs (short buffers skip it)
    if(len > 15)
    {
      __m512 sumr2 = sumr;
      __m512 sumi2 = sumr;
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...]
        // accumulate the real and imaginary parts of the dot-product separately
        sumr = _mm512_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
        sumi = _mm512_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1);
        sumr2 = _mm512_fmadd_ps(ld1, ld2, sumr2);
        sumi2 = _mm512_fmadd_ps(sh, ld2, sumi2);
      }
      sumr = _mm512_add_ps(sumr, sumr2);
      sumi = _mm512_add_ps(sumi, sumi2);
    }
    // up to 1 full register remains if bit 3 of len is set
    if(len&8) // process 1 register (8 complex elements)
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      // accumulate the real and imaginary parts of the dot-product separately
      sumr = _mm512_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
      sumi = _mm512_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
      i += 8;
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
      // accumulate the real and imaginary parts of the dot-product separately
      sumr = _mm512_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
      sumi = _mm512_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
    }
    // horizontal reduction, interleaving the accumulators lines the 4 product
    // sums up as [ArBr, AiBr, AiBi, ArBi] in every 128-bit lane
    ld1 = _mm512_unpacklo_ps(sumr, sumi);// [r0, i0, r1, i1,.,. | r8, i8, ...]
    ld2 = _mm512_unpackhi_ps(sumr, sumi);// [r2, i2, r3, i3,... | r10, i10, ...]
    ld1 = _mm512_add_ps(ld1, ld2);  // [r0+r2, i0+i2, ... | r8+r10, i8+i10, ...]
    // add the upper half of ld1 to the lower half (cast through double to
    // avoid requiring the AVX512DQ instruction _mm512_extractf32x8_ps)
    hs2 = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(ld1), 1));
    hs2 = _mm256_add_ps(_mm512_castps512_ps256(ld1), hs2);
    // fold and sum the halves of hs2 into hs1
    hs1 = _mm256_extractf128_ps(hs2, 1); // extract upper half of register
    hs1 = _mm_add_ps(_mm256_castps256_ps128(hs2), hs1); // [ArBr,AiBr,AiBi,ArBi]
    // [ArBr-AiBi, AiBr+ArBi, ...], real and imaginary parts of the result
    hs1 = _mm_addsub_ps(hs1, _mm_permute_ps(hs1, 0x4e));
    // store both parts with one 8-byte store to avoid instruction stall with
    // certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs1);
    return;
  }
  #endif // end AVX512 complex*complex dotp
  #if !defined(DISABLE_AVX2) // AVX2 complex*complex dotp
  __attribute__((__target__("avx2,fma")))
  inline void dotp(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m256 ld1, ld2, sh;
    __m256 sumr = _mm256_set1_ps(0.0f); // real part of dot-product
    __m256 sumi = sumr;                 // imaginary part of dot-product
    __m128 hs;  // horizontal sum to finish dot-product
    int i = 0;
    // note: with a single pair of accumulators every FMA waits on the previous
    // one, so longer buffers use 2 pairs to keep the FMA units busy, combining
    // them only if the 2 register loop runs (short buffers skip it)
    if(len > 7)
    {
      __m256 sumr2 = sumr;
      __m256 sumi2 = sumr;
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...]
        // accumulate the real and imaginary parts of the dot-product separately
        sumr = _mm256_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
        sumi = _mm256_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1);
        sumr2 = _mm256_fmadd_ps(ld1, ld2, sumr2);
        sumi2 = _mm256_fmadd_ps(sh, ld2, sumi2);
      }
      sumr = _mm256_add_ps(sumr, sumr2);
      sumi = _mm256_add_ps(sumi, sumi2);
    }
    // up to 1 full register remains if bit 2 of len is set
    if(len&4) // process 1 register (4 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      // accumulate the real and imaginary parts of the dot-product separately
      sumr = _mm256_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
      sumi = _mm256_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
      i += 4;
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]), msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]), msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      // accumulate the real and imaginary parts of the dot-product separately
      sumr = _mm256_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
      sumi = _mm256_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
    }
    // horizontal reduction, interleaving the accumulators lines the 4 product
    // sums up as [ArBr, AiBr, AiBi, ArBi] in each 128-bit lane
    ld1 = _mm256_unpacklo_ps(sumr, sumi); // [r0, i0, r1, i1, r4, i4, r5, i5]
    ld2 = _mm256_unpackhi_ps(sumr, sumi); // [r2, i2, r3, i3, r6, i6, r7, i7]
    ld1 = _mm256_add_ps(ld1, ld2);  // [r0+r2, i0+i2, r1+r3, i1+i3, ...]
    // fold and sum the halves of ld1 into hs
    hs = _mm256_extractf128_ps(ld1, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(ld1), hs); // [ArBr,AiBr,AiBi,ArBi]
    // [ArBr-AiBi, AiBr+ArBi, ...], real and imaginary parts of the result
    hs = _mm_addsub_ps(hs, _mm_permute_ps(hs, 0x4e));
    // store both parts with one 8-byte store to avoid instruction stall with
    // certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs);
    return;
  }
  #endif // end AVX2 complex*complex dotp
  #if !defined(DISABLE_AVX) // AVX complex*complex dotp
  __attribute__((__target__("avx")))
  inline void dotp(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m256 ld1, ld2, sh;
    __m256 sumr = _mm256_set1_ps(0.0f); // real part of dot-product
    __m256 sumi = sumr;                 // imaginary part of dot-product
    __m128 hs;  // horizontal sum to finish dot-product
    int i = 0;
    // note: with a single pair of accumulators every mul/add waits on the
    // previous one, so longer buffers use 2 pairs to keep the pipeline full,
    // combining them only if the 2 register loop runs (short buffers skip it)
    if(len > 7)
    {
      __m256 sumr2 = sumr;
      __m256 sumi2 = sumr;
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...]
        // accumulate the real and imaginary parts of the dot-product separately
        ld1 = _mm256_mul_ps(ld1, ld2);  // [ArBr, AiBi]
        ld2 = _mm256_mul_ps(sh, ld2);   // [AiBr, ArBi]
        sumr = _mm256_add_ps(sumr, ld1);
        sumi = _mm256_add_ps(sumi, ld2);
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1);
        ld1 = _mm256_mul_ps(ld1, ld2);
        ld2 = _mm256_mul_ps(sh, ld2);
        sumr2 = _mm256_add_ps(sumr2, ld1);
        sumi2 = _mm256_add_ps(sumi2, ld2);
      }
      sumr = _mm256_add_ps(sumr, sumr2);
      sumi = _mm256_add_ps(sumi, sumi2);
    }
    // up to 1 full register remains if bit 2 of len is set
    if(len&4) // process 1 register (4 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      // accumulate the real and imaginary parts of the dot-product separately
      ld1 = _mm256_mul_ps(ld1, ld2);  // [ArBr, AiBi]
      ld2 = _mm256_mul_ps(sh, ld2);   // [AiBr, ArBi]
      sumr = _mm256_add_ps(sumr, ld1);
      sumi = _mm256_add_ps(sumi, ld2);
      i += 4;
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]), msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]), msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      // accumulate the real and imaginary parts of the dot-product separately
      ld1 = _mm256_mul_ps(ld1, ld2);  // [ArBr, AiBi]
      ld2 = _mm256_mul_ps(sh, ld2);   // [AiBr, ArBi]
      sumr = _mm256_add_ps(sumr, ld1);
      sumi = _mm256_add_ps(sumi, ld2);
    }
    // horizontal reduction, interleaving the accumulators lines the 4 product
    // sums up as [ArBr, AiBr, AiBi, ArBi] in each 128-bit lane
    ld1 = _mm256_unpacklo_ps(sumr, sumi); // [r0, i0, r1, i1, r4, i4, r5, i5]
    ld2 = _mm256_unpackhi_ps(sumr, sumi); // [r2, i2, r3, i3, r6, i6, r7, i7]
    ld1 = _mm256_add_ps(ld1, ld2);  // [r0+r2, i0+i2, r1+r3, i1+i3, ...]
    // fold and sum the halves of ld1 into hs
    hs = _mm256_extractf128_ps(ld1, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(ld1), hs); // [ArBr,AiBr,AiBi,ArBi]
    // [ArBr-AiBi, AiBr+ArBi, ...], real and imaginary parts of the result
    hs = _mm_addsub_ps(hs, _mm_permute_ps(hs, 0x4e));
    // store both parts with one 8-byte store to avoid instruction stall with
    // certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs);
    return;
  }
  #endif // end AVX complex*complex dotp
  __attribute__((__target__("default"))) // default complex*complex dotp
  inline void dotp(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out)
  {
    float reSum = 0.0f;
    float imSum = 0.0f;
    for(int i=0; i<len; ++i)
    {
      reSum += (in1[i].real()*in2[i].real()-in1[i].imag()*in2[i].imag());
      imSum += (in1[i].real()*in2[i].imag()+in1[i].imag()*in2[i].real()); 
    }
    out.real(reSum);
    out.imag(imSum);
    return;
  }
  // wrapper returning result by value, must be defined after every target
  // version so the calls go through the runtime dispatcher
  inline complex<float> dotp(complex<float> const * const in1,
      complex<float> const * const in2, const int len)
  {
    complex<float> out;
    dotp(in1, in2, len, out);
    return out;
  }

  // real*complex dot-product
  #if !defined(DISABLE_AVX512) // AVX512 real*complex dotp
  __attribute__((__target__("avx512f")))
  inline void dotp(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m512 ld1, ld2, ld3, sc;
    __m512 sum1 = _mm512_set1_ps(0.0f);
    __m512 sum2 = sum1;
    // constants for duplicating real elements
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // registers for horizontal reduction
    __m256 hs2;
    __m128 hs1;
    int i = 0;
    // note: with 2 accumulators every FMA waits on the previous one, so longer
    // buffers use 4 to keep the FMA unit busy
    if(len > 31)
    {
      __m512 sum3 = sum1;
      __m512 sum4 = sum1;
      for(; i<len-31; i+=32) // process 2 registers (32 real elements)
      {
        ld1 = _mm512_loadu_ps(&in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        sc = _mm512_permutexvar_ps(p1, ld1);
        ld1 = _mm512_permutexvar_ps(p2, ld1);
        sum1 = _mm512_fmadd_ps(ld2, sc, sum1);
        sum2 = _mm512_fmadd_ps(ld3, ld1, sum2);
        ld1 = _mm512_loadu_ps(&in1[i+16]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+16]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+24]));
        sc = _mm512_permutexvar_ps(p1, ld1);
        ld1 = _mm512_permutexvar_ps(p2, ld1);
        sum3 = _mm512_fmadd_ps(ld2, sc, sum3);
        sum4 = _mm512_fmadd_ps(ld3, ld1, sum4);
      }
      sum1 = _mm512_add_ps(sum1, sum3);
      sum2 = _mm512_add_ps(sum2, sum4);
    }
    // at least 1 full register of reals (2 of complex) remain if bit 4 of len
    // is set
    if(len&16) // process 1 real/2 complex registers (16 real elements)
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      sum1 = _mm512_fmadd_ps(ld2, sc, sum1);
      sum2 = _mm512_fmadd_ps(ld3, ld1, sum2);
      i += 16;
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
      sum1 = _mm512_fmadd_ps(ld2, sc, sum1);
      sum2 = _mm512_fmadd_ps(ld3, ld1, sum2);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      sum1 = _mm512_fmadd_ps(ld2, sc, sum1);
    }
    // combine both sums and reduce horizontally
    sum1 = _mm512_add_ps(sum1, sum2);
    // add the upper half of sum1 to the lower half (cast through double to
    // avoid requiring the AVX512DQ instruction _mm512_extractf32x8_ps)
    hs2 = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(sum1), 1));
    hs2 = _mm256_add_ps(_mm512_castps512_ps256(sum1), hs2);
    // add the upper half of hs2 to the lower half
    hs1 = _mm256_extractf128_ps(hs2, 1); // extract upper half of register
    hs1 = _mm_add_ps(_mm256_castps256_ps128(hs2), hs1);
    // add the upper 2 FP to the lower 2 in hs1
    hs1 = _mm_add_ps(hs1, _mm_permute_ps(hs1, 0x4e));
    // store both real and imaginary parts with one 8-byte store to avoid
    // instruction stall with certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs1);
    return;
  }
  #endif // end AVX512 real*complex dotp
  #if !defined(DISABLE_AVX2) // AVX2 real*complex dotp
  __attribute__((__target__("avx2,fma")))
  inline void dotp(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m256 ld1, ld2, ld3, sc;
    __m256 sum1 = _mm256_set1_ps(0.0f);
    __m256 sum2 = sum1;
    // constants for duplicating real elements
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    // register for horizontal reduction
    __m128 hs;
    int i = 0;
    // note: with 2 accumulators every FMA waits on the previous one, so longer
    // buffers use 4 to keep the FMA units busy
    if(len > 15)
    {
      __m256 sum3 = sum1;
      __m256 sum4 = sum1;
      for(; i<len-15; i+=16) // process 2 registers (16 real elements)
      {
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        sum1 = _mm256_fmadd_ps(ld2, sc, sum1);
        sum2 = _mm256_fmadd_ps(ld3, ld1, sum2);
        ld1 = _mm256_loadu_ps(&in1[i+8]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+12]));
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        sum3 = _mm256_fmadd_ps(ld2, sc, sum3);
        sum4 = _mm256_fmadd_ps(ld3, ld1, sum4);
      }
      sum1 = _mm256_add_ps(sum1, sum3);
      sum2 = _mm256_add_ps(sum2, sum4);
    }
    // at least 1 full register of reals (2 of complex) remain if bit 3 of len
    // is set
    if(len&8) // process 1 real/2 complex registers (8 real elements)
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      sum1 = _mm256_fmadd_ps(ld2, sc, sum1);
      sum2 = _mm256_fmadd_ps(ld3, ld1, sum2);
      i += 8;
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
        sum1 = _mm256_fmadd_ps(ld2, sc, sum1);
        sum2 = _mm256_fmadd_ps(ld3, ld1, sum2);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        sum1 = _mm256_fmadd_ps(ld2, sc, sum1);
      }
    }
    // combine both sums and reduce horizontally
    sum1 = _mm256_add_ps(sum1, sum2);
    // add the upper half and lower half of sum1 and combine into hs
    hs = _mm256_extractf128_ps(sum1, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(sum1), hs); // add to lower half
    // add the upper 2 FP to the lower 2 in hs
    hs = _mm_add_ps(hs, _mm_permute_ps(hs, 0x4e));
    // store both real and imaginary parts with one 8-byte store to avoid
    // instruction stall with certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs);
    return;
  }
  #endif // end AVX2 real*complex dotp
  #if !defined(DISABLE_AVX) // AVX real*complex dotp
  __attribute__((__target__("avx")))
  inline void dotp(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m256 ld1, ld2, ld3, sc1, sc2;
    __m256 sum1 = _mm256_set1_ps(0.0f);
    __m256 sum2 = sum1;
    // in-lane index for duplicating real elements, used on 4 reals broadcast
    // to both 128-bit lanes: [0,0,1,1] in the lower lane, [2,2,3,3] in upper
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    // register for horizontal reduction
    __m128 hs;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      // the broadcasts are plain loads, which leaves 2 shuffles per register
      // instead of 4 (improves performance with gcc, intel unaffected)
      sc1 = _mm256_permutevar_ps(_mm256_broadcast_ps(
          reinterpret_cast<__m128 const *>(&in1[i])), p1);  // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permutevar_ps(_mm256_broadcast_ps(
          reinterpret_cast<__m128 const *>(&in1[i+4])), p1);// [4,4,5,5,6,6,7,7]
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      sum1 = _mm256_add_ps(ld2, sum1);
      sum2 = _mm256_add_ps(ld3, sum2);
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
        sum1 = _mm256_add_ps(ld2, sum1);
        sum2 = _mm256_add_ps(ld3, sum2);
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
        sum1 = _mm256_add_ps(ld2, sum1);
      }
    }
    // combine both sums and reduce horizontally
    sum1 = _mm256_add_ps(sum1, sum2);
    // add the upper half and lower half of sum1 and combine into hs
    hs = _mm256_extractf128_ps(sum1, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(sum1), hs); // add to lower half
    // add the upper 2 FP to the lower 2 in hs
    hs = _mm_add_ps(hs, _mm_permute_ps(hs, 0x4e));
    // store both real and imaginary parts with one 8-byte store to avoid
    // instruction stall with certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs);
    return;
  }
  #endif // end AVX real*complex dotp
  __attribute__((__target__("default"))) // default real*complex dotp
  inline void dotp(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> &out)
  {
    float reSum = 0.0f;
    float imSum = 0.0f;
    for(int i=0; i<len; ++i)
    {
      reSum += in1[i]*in2[i].real();
      imSum += in1[i]*in2[i].imag(); 
    }
    out.real(reSum);
    out.imag(imSum);
    return;
  }
  // wrapper returning result by value, must be defined after every target
  // version so the calls go through the runtime dispatcher
  inline complex<float> dotp(float const * const in1,
      complex<float> const * const in2, const int len)
  {
    complex<float> out;
    dotp(in1, in2, len, out);
    return out;
  }

  // complex*complex conjugate dot-product
  #if !defined(DISABLE_AVX512) // AVX512 complex*conj(complex) dotp
  __attribute__((__target__("avx512f")))
  inline void dotpc(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m512 ld1, ld2, sh;
    __m512 sumr = _mm512_set1_ps(0.0f); // real part of dot-product
    __m512 sumi = sumr;                 // imaginary part of dot-product
    __m256 hs2;  // horizontal sum, 8 elements (fold and sum halves of m512)
    __m128 hs1;  // horizontal sum to finish dot-product (fold and sum m256)
    int i = 0;
    // note: with a single pair of accumulators every FMA waits on the previous
    // one, so longer buffers use 2 pairs to keep the FMA unit busy, combining
    // them only if the 2 register loop runs (short buffers skip it)
    if(len > 15)
    {
      __m512 sumr2 = sumr;
      __m512 sumi2 = sumr;
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...]
        // accumulate the real and imaginary parts of the dot-product separately
        sumr = _mm512_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
        sumi = _mm512_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i+8]));
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        sh = _mm512_shuffle_ps(ld1, ld1, 0xb1);
        sumr2 = _mm512_fmadd_ps(ld1, ld2, sumr2);
        sumi2 = _mm512_fmadd_ps(sh, ld2, sumi2);
      }
      sumr = _mm512_add_ps(sumr, sumr2);
      sumi = _mm512_add_ps(sumi, sumi2);
    }
    // up to 1 full register remains if bit 3 of len is set
    if(len&8) // process 1 register (8 complex elements)
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      // accumulate the real and imaginary parts of the dot-product separately
      sumr = _mm512_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
      sumi = _mm512_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
      i += 8;
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
      // accumulate the real and imaginary parts of the dot-product separately
      sumr = _mm512_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
      sumi = _mm512_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
    }
    // horizontal reduction, interleaving the accumulators lines the 4 product
    // sums up as [ArBr, AiBr, AiBi, ArBi] in every 128-bit lane
    ld1 = _mm512_unpacklo_ps(sumr, sumi);// [r0, i0, r1, i1,.,. | r8, i8, ...]
    ld2 = _mm512_unpackhi_ps(sumr, sumi);// [r2, i2, r3, i3,... | r10, i10, ...]
    ld1 = _mm512_add_ps(ld1, ld2);  // [r0+r2, i0+i2, ... | r8+r10, i8+i10, ...]
    // add the upper half of ld1 to the lower half (cast through double to
    // avoid requiring the AVX512DQ instruction _mm512_extractf32x8_ps)
    hs2 = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(ld1), 1));
    hs2 = _mm256_add_ps(_mm512_castps512_ps256(ld1), hs2);
    // fold and sum the halves of hs2 into hs1
    hs1 = _mm256_extractf128_ps(hs2, 1); // extract upper half of register
    hs1 = _mm_add_ps(_mm256_castps256_ps128(hs2), hs1); // [ArBr,AiBr,AiBi,ArBi]
    // swap upper and lower halves of hs1 and negate the ArBi term for conjugate
    // [ArBr+AiBi, AiBr-ArBi, ...], real and imaginary parts of the result
    hs1 = _mm_add_ps(hs1, _mm_xor_ps(_mm_permute_ps(hs1, 0x4e),
        _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f)));
    // store both parts with one 8-byte store to avoid instruction stall with
    // certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs1);
    return;
  }
  #endif // end AVX512 complex*conj(complex) dotp
  #if !defined(DISABLE_AVX2) // AVX2 complex*conj(complex) dotp
  __attribute__((__target__("avx2,fma")))
  inline void dotpc(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m256 ld1, ld2, sh;
    __m256 sumr = _mm256_set1_ps(0.0f); // real part of dot-product
    __m256 sumi = sumr;                 // imaginary part of dot-product
    __m128 hs;  // horizontal sum to finish dot-product
    int i = 0;
    // note: with a single pair of accumulators every FMA waits on the previous
    // one, so longer buffers use 2 pairs to keep the FMA units busy, combining
    // them only if the 2 register loop runs (short buffers skip it)
    if(len > 7)
    {
      __m256 sumr2 = sumr;
      __m256 sumi2 = sumr;
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...]
        // accumulate the real and imaginary parts of the dot-product separately
        sumr = _mm256_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
        sumi = _mm256_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1);
        sumr2 = _mm256_fmadd_ps(ld1, ld2, sumr2);
        sumi2 = _mm256_fmadd_ps(sh, ld2, sumi2);
      }
      sumr = _mm256_add_ps(sumr, sumr2);
      sumi = _mm256_add_ps(sumi, sumi2);
    }
    // up to 1 full register remains if bit 2 of len is set
    if(len&4) // process 1 register (4 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      // accumulate the real and imaginary parts of the dot-product separately
      sumr = _mm256_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
      sumi = _mm256_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
      i += 4;
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]), msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]), msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      // accumulate the real and imaginary parts of the dot-product separately
      sumr = _mm256_fmadd_ps(ld1, ld2, sumr); // [ArBr, AiBi]
      sumi = _mm256_fmadd_ps(sh, ld2, sumi);  // [AiBr, ArBi]
    }
    // horizontal reduction, interleaving the accumulators lines the 4 product
    // sums up as [ArBr, AiBr, AiBi, ArBi] in each 128-bit lane
    ld1 = _mm256_unpacklo_ps(sumr, sumi); // [r0, i0, r1, i1, r4, i4, r5, i5]
    ld2 = _mm256_unpackhi_ps(sumr, sumi); // [r2, i2, r3, i3, r6, i6, r7, i7]
    ld1 = _mm256_add_ps(ld1, ld2);  // [r0+r2, i0+i2, r1+r3, i1+i3, ...]
    // fold and sum the halves of ld1 into hs
    hs = _mm256_extractf128_ps(ld1, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(ld1), hs); // [ArBr,AiBr,AiBi,ArBi]
    // swap upper and lower halves of hs and negate the ArBi term for conjugate
    // [ArBr+AiBi, AiBr-ArBi, ...], real and imaginary parts of the result
    hs = _mm_add_ps(hs, _mm_xor_ps(_mm_permute_ps(hs, 0x4e),
        _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f)));
    // store both parts with one 8-byte store to avoid instruction stall with
    // certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs);
    return;
  }
  #endif // end AVX2 complex*conj(complex) dotp
  #if !defined(DISABLE_AVX) // AVX complex*conj(complex) dotp
  __attribute__((__target__("avx")))
  inline void dotpc(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m256 ld1, ld2, sh;
    __m256 sumr = _mm256_set1_ps(0.0f); // real part of dot-product
    __m256 sumi = sumr;                 // imaginary part of dot-product
    __m128 hs;  // horizontal sum to finish dot-product
    int i = 0;
    // note: with a single pair of accumulators every mul/add waits on the
    // previous one, so longer buffers use 2 pairs to keep the pipeline full,
    // combining them only if the 2 register loop runs (short buffers skip it)
    if(len > 7)
    {
      __m256 sumr2 = sumr;
      __m256 sumi2 = sumr;
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...]
        // accumulate the real and imaginary parts of the dot-product separately
        ld1 = _mm256_mul_ps(ld1, ld2);  // [ArBr, AiBi]
        ld2 = _mm256_mul_ps(sh, ld2);   // [AiBr, ArBi]
        sumr = _mm256_add_ps(sumr, ld1);
        sumi = _mm256_add_ps(sumi, ld2);
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i+4]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        sh = _mm256_shuffle_ps(ld1, ld1, 0xb1);
        ld1 = _mm256_mul_ps(ld1, ld2);
        ld2 = _mm256_mul_ps(sh, ld2);
        sumr2 = _mm256_add_ps(sumr2, ld1);
        sumi2 = _mm256_add_ps(sumi2, ld2);
      }
      sumr = _mm256_add_ps(sumr, sumr2);
      sumi = _mm256_add_ps(sumi, sumi2);
    }
    // up to 1 full register remains if bit 2 of len is set
    if(len&4) // process 1 register (4 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      // accumulate the real and imaginary parts of the dot-product separately
      ld1 = _mm256_mul_ps(ld1, ld2);  // [ArBr, AiBi]
      ld2 = _mm256_mul_ps(sh, ld2);   // [AiBr, ArBi]
      sumr = _mm256_add_ps(sumr, ld1);
      sumi = _mm256_add_ps(sumi, ld2);
      i += 4;
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in1[i]), msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]), msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      // accumulate the real and imaginary parts of the dot-product separately
      ld1 = _mm256_mul_ps(ld1, ld2);  // [ArBr, AiBi]
      ld2 = _mm256_mul_ps(sh, ld2);   // [AiBr, ArBi]
      sumr = _mm256_add_ps(sumr, ld1);
      sumi = _mm256_add_ps(sumi, ld2);
    }
    // horizontal reduction, interleaving the accumulators lines the 4 product
    // sums up as [ArBr, AiBr, AiBi, ArBi] in each 128-bit lane
    ld1 = _mm256_unpacklo_ps(sumr, sumi); // [r0, i0, r1, i1, r4, i4, r5, i5]
    ld2 = _mm256_unpackhi_ps(sumr, sumi); // [r2, i2, r3, i3, r6, i6, r7, i7]
    ld1 = _mm256_add_ps(ld1, ld2);  // [r0+r2, i0+i2, r1+r3, i1+i3, ...]
    // fold and sum the halves of ld1 into hs
    hs = _mm256_extractf128_ps(ld1, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(ld1), hs); // [ArBr,AiBr,AiBi,ArBi]
    // swap upper and lower halves of hs and negate the ArBi term for conjugate
    // [ArBr+AiBi, AiBr-ArBi, ...], real and imaginary parts of the result
    hs = _mm_add_ps(hs, _mm_xor_ps(_mm_permute_ps(hs, 0x4e),
        _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f)));
    // store both parts with one 8-byte store to avoid instruction stall with
    // certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs);
    return;
  }
  #endif // end AVX complex*conj(complex) dotp
  __attribute__((__target__("default"))) // default complex*conj(complex) dotp
  inline void dotpc(complex<float> const * const in1,
      complex<float> const * const in2, const int len, complex<float> &out)
  {
    float reSum = 0.0f;
    float imSum = 0.0f;
    for(int i=0; i<len; ++i)
    {
      reSum += (in1[i].real()*in2[i].real()+in1[i].imag()*in2[i].imag());
      imSum += (in1[i].imag()*in2[i].real()-in1[i].real()*in2[i].imag());
    }
    out.real(reSum);
    out.imag(imSum);
    return;
  }
  // wrapper returning result by value, must be defined after every target
  // version so the calls go through the runtime dispatcher
  inline complex<float> dotpc(complex<float>  const * const in1,
      complex<float> const * const in2, const int len)
  {
    complex<float> out;
    dotpc(in1, in2, len, out);
    return out;
  }
  
  // real*complex conjugate dot-product
  #if !defined(DISABLE_AVX512) // AVX512 real*conj(complex) dotp
  __attribute__((__target__("avx512f")))
  inline void dotpc(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m512 ld1, ld2, ld3, sc;
    __m512 sum1 = _mm512_set1_ps(0.0f);
    __m512 sum2 = sum1;
    // constants for duplicating real elements
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    // registers for horizontal reduction
    __m256 hs2;
    __m128 hs1;
    int i = 0;
    // note: with 2 accumulators every FMA waits on the previous one, so longer
    // buffers use 4 to keep the FMA unit busy
    if(len > 31)
    {
      __m512 sum3 = sum1;
      __m512 sum4 = sum1;
      for(; i<len-31; i+=32) // process 2 registers (32 real elements)
      {
        ld1 = _mm512_loadu_ps(&in1[i]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        sc = _mm512_permutexvar_ps(p1, ld1);
        ld1 = _mm512_permutexvar_ps(p2, ld1);
        sum1 = _mm512_fmadd_ps(ld2, sc, sum1);
        sum2 = _mm512_fmadd_ps(ld3, ld1, sum2);
        ld1 = _mm512_loadu_ps(&in1[i+16]);
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+16]));
        ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+24]));
        sc = _mm512_permutexvar_ps(p1, ld1);
        ld1 = _mm512_permutexvar_ps(p2, ld1);
        sum3 = _mm512_fmadd_ps(ld2, sc, sum3);
        sum4 = _mm512_fmadd_ps(ld3, ld1, sum4);
      }
      sum1 = _mm512_add_ps(sum1, sum3);
      sum2 = _mm512_add_ps(sum2, sum4);
    }
    // at least 1 full register of reals (2 of complex) remain if bit 4 of len
    // is set
    if(len&16) // process 1 real/2 complex registers (16 real elements)
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      sum1 = _mm512_fmadd_ps(ld2, sc, sum1);
      sum2 = _mm512_fmadd_ps(ld3, ld1, sum2);
      i += 16;
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
      sum1 = _mm512_fmadd_ps(ld2, sc, sum1);
      sum2 = _mm512_fmadd_ps(ld3, ld1, sum2);
    }
    else if(rem)
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps(mask16(rem), &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in2[i]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      sum1 = _mm512_fmadd_ps(ld2, sc, sum1);
    }
    // combine both sums and reduce horizontally
    sum1 = _mm512_add_ps(sum1, sum2);
    // add the upper half of sum1 to the lower half (cast through double to
    // avoid requiring the AVX512DQ instruction _mm512_extractf32x8_ps)
    hs2 = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(sum1), 1));
    hs2 = _mm256_add_ps(_mm512_castps512_ps256(sum1), hs2);
    // add the upper half of hs2 to the lower half
    hs1 = _mm256_extractf128_ps(hs2, 1); // extract upper half of register
    hs1 = _mm_add_ps(_mm256_castps256_ps128(hs2), hs1);
    // add the upper 2 FP to the lower 2 in hs1
    hs1 = _mm_add_ps(hs1, _mm_permute_ps(hs1, 0x4e));
    // flip the sign on the imaginary part for conjugate
    hs1 = _mm_xor_ps(hs1, _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f));
    // store both real and imaginary parts with one 8-byte store to avoid
    // instruction stall with certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs1);
    return;
  }
  #endif // end AVX512 real*conj(complex) dotp
  #if !defined(DISABLE_AVX2) // AVX2 real*conj(complex) dotp
  __attribute__((__target__("avx2,fma")))
  inline void dotpc(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m256 ld1, ld2, ld3, sc;
    __m256 sum1 = _mm256_set1_ps(0.0f);
    __m256 sum2 = sum1;
    // constants for duplicating real elements
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    // register for horizontal reduction
    __m128 hs;
    int i = 0;
    // note: with 2 accumulators every FMA waits on the previous one, so longer
    // buffers use 4 to keep the FMA units busy
    if(len > 15)
    {
      __m256 sum3 = sum1;
      __m256 sum4 = sum1;
      for(; i<len-15; i+=16) // process 2 registers (16 real elements)
      {
        ld1 = _mm256_loadu_ps(&in1[i]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        sum1 = _mm256_fmadd_ps(ld2, sc, sum1);
        sum2 = _mm256_fmadd_ps(ld3, ld1, sum2);
        ld1 = _mm256_loadu_ps(&in1[i+8]);
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+8]));
        ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+12]));
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        sum3 = _mm256_fmadd_ps(ld2, sc, sum3);
        sum4 = _mm256_fmadd_ps(ld3, ld1, sum4);
      }
      sum1 = _mm256_add_ps(sum1, sum3);
      sum2 = _mm256_add_ps(sum2, sum4);
    }
    // at least 1 full register of reals (2 of complex) remain if bit 3 of len
    // is set
    if(len&8) // process 1 real/2 complex registers (8 real elements)
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      sum1 = _mm256_fmadd_ps(ld2, sc, sum1);
      sum2 = _mm256_fmadd_ps(ld3, ld1, sum2);
      i += 8;
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
        sum1 = _mm256_fmadd_ps(ld2, sc, sum1);
        sum2 = _mm256_fmadd_ps(ld3, ld1, sum2);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const *>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in2[i]),
            msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        sum1 = _mm256_fmadd_ps(ld2, sc, sum1);
      }
    }
    // combine both sums and reduce horizontally
    sum1 = _mm256_add_ps(sum1, sum2);
    // add the upper half and lower half of sum1 and combine into hs
    hs = _mm256_extractf128_ps(sum1, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(sum1), hs); // add to lower half
    // add the upper 2 FP to the lower 2 in hs
    hs = _mm_add_ps(hs, _mm_permute_ps(hs, 0x4e));
    // flip the sign on the imaginary part for conjugate
    hs = _mm_xor_ps(hs, _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f));
    // store both real and imaginary parts with one 8-byte store to avoid
    // instruction stall with certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs);
    return;
  }
  #endif // end AVX2 real*conj(complex) dotp
  #if !defined(DISABLE_AVX) // AVX real*conj(complex) dotp
  __attribute__((__target__("avx")))
  inline void dotpc(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> &out)
  {
    if(len<=0)
    {
      out.real(0.0f);
      out.imag(0.0f);
      return;
    }
    __m256 ld1, ld2, ld3, sc1, sc2;
    __m256 sum1 = _mm256_set1_ps(0.0f);
    __m256 sum2 = sum1;
    // in-lane index for duplicating real elements, used on 4 reals broadcast
    // to both 128-bit lanes: [0,0,1,1] in the lower lane, [2,2,3,3] in upper
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    // register for horizontal reduction
    __m128 hs;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      // the broadcasts are plain loads, which leaves 2 shuffles per register
      // instead of 4 (improves performance with gcc, intel unaffected)
      sc1 = _mm256_permutevar_ps(_mm256_broadcast_ps(
          reinterpret_cast<__m128 const *>(&in1[i])), p1);  // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permutevar_ps(_mm256_broadcast_ps(
          reinterpret_cast<__m128 const *>(&in1[i+4])), p1);// [4,4,5,5,6,6,7,7]
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in2[i+4]));
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      sum1 = _mm256_add_ps(ld2, sum1);
      sum2 = _mm256_add_ps(ld3, sum2);
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
        sum1 = _mm256_add_ps(ld2, sum1);
        sum2 = _mm256_add_ps(ld3, sum2);
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
        sum1 = _mm256_add_ps(ld2, sum1);
      }
    }
    // combine both sums and reduce horizontally
    sum1 = _mm256_add_ps(sum1, sum2);
    // add the upper half and lower half of sum1 and combine into hs
    hs = _mm256_extractf128_ps(sum1, 1); // extract upper half of register
    hs = _mm_add_ps(_mm256_castps256_ps128(sum1), hs); // add to lower half
    // add the upper 2 FP to the lower 2 in hs
    hs = _mm_add_ps(hs, _mm_permute_ps(hs, 0x4e));
    // flip the sign on the imaginary part for conjugate
    hs = _mm_xor_ps(hs, _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f));
    // store both real and imaginary parts with one 8-byte store to avoid
    // instruction stall with certain compilers
    _mm_storel_pi(reinterpret_cast<__m64 *>(&out), hs);
    return;
  }
  #endif // end AVX real*conj(complex) dotp
  __attribute__((__target__("default"))) // default real*conj(complex) dotp
  inline void dotpc(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> &out)
  {
    float reSum = 0.0f;
    float imSum = 0.0f;
    for(int i=0; i<len; ++i)
    {
      reSum += in1[i]*in2[i].real();
      imSum -= in1[i]*in2[i].imag(); 
    }
    out.real(reSum);
    out.imag(imSum);
    return;
  }
  // wrapper returning result by value, must be defined after every target
  // version so the calls go through the runtime dispatcher
  inline complex<float> dotpc(float const * const in1,
      complex<float> const * const in2, const int len)
  {
    complex<float> out;
    dotpc(in1, in2, len, out);
    return out;
  }

  //========================================//
  // Raise to Power Function Implementation //
  //========================================//
  // squared real value
  #if !defined(DISABLE_AVX512) // AVX512 real^2
  __attribute__((__target__("avx512f")))
  inline void pow2(float const * const in, const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in[i]);
      ld1 = _mm512_mul_ps(ld1, ld1);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in[i]);
      ld1 = _mm512_mul_ps(ld1, ld1);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // end AVX512 real^2
  #if !defined(DISABLE_AVX) // AVX real^2
  __attribute__((__target__("avx")))
  inline void pow2(float const * const in, const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in[i]);
      ld1 = _mm256_mul_ps(ld1, ld1);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld1);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // end AVX real^2
  __attribute__((__target__("default"))) // default real^2
  inline void pow2(float const * const in, const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in[i]*in[i];
  }
  

  // squared complex value
  #if !defined(DISABLE_AVX512) // AVX512 complex^2
  __attribute__((__target__("avx512f")))
  inline void pow2(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, re, im;
    int i = 0;
    // note: 2 registers are split into separate real and imaginary registers
    // so the squaring needs no shuffles, then interleaved again for the store
    for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm512_unpacklo_ps(re, im); // A^2
      ld2 = _mm512_unpackhi_ps(re, im); // B^2
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld2);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem > 8) // one full and one partial register
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in[i+8])); // B
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm512_unpacklo_ps(re, im); // A^2
      ld2 = _mm512_unpackhi_ps(re, im); // B^2
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk, ld2);
    }
    else if(rem) // one partial register, split against itself
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in[i])); // A
      re = _mm512_shuffle_ps(ld1, ld1, 0x88); // [Ar0,Ar1,Ar0,Ar1,...]
      im = _mm512_shuffle_ps(ld1, ld1, 0xdd); // [Ai0,Ai1,Ai0,Ai1,...]
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm512_unpacklo_ps(re, im); // A^2
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // end AVX512 complex^2
  #if !defined(DISABLE_AVX2) // AVX2 complex^2
  __attribute__((__target__("avx2,fma")))
  inline void pow2(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, re, im;
    int i = 0;
    // note: 2 registers are split into separate real and imaginary registers
    // so the squaring needs no shuffles, then interleaved again for the store
    for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^2
      ld2 = _mm256_unpackhi_ps(re, im); // B^2
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem > 4) // one full and one partial register
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), msk);
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^2
      ld2 = _mm256_unpackhi_ps(re, im); // B^2
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk, ld2);
    }
    else if(rem) // one partial register, split against itself
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      re = _mm256_shuffle_ps(ld1, ld1, 0x88); // [Ar0,Ar1,Ar0,Ar1,...]
      im = _mm256_shuffle_ps(ld1, ld1, 0xdd); // [Ai0,Ai1,Ai0,Ai1,...]
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^2
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex^2
  #if !defined(DISABLE_AVX) // AVX complex^2
  __attribute__((__target__("avx")))
  inline void pow2(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, re, im;
    int i = 0;
    // note: 2 registers are split into separate real and imaginary registers
    // so the squaring needs no shuffles, then interleaved again for the store
    for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^2
      ld2 = _mm256_unpackhi_ps(re, im); // B^2
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem > 4) // one full and one partial register
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), msk);
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^2
      ld2 = _mm256_unpackhi_ps(re, im); // B^2
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk, ld2);
    }
    else if(rem) // one partial register, split against itself
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      re = _mm256_shuffle_ps(ld1, ld1, 0x88); // [Ar0,Ar1,Ar0,Ar1,...]
      im = _mm256_shuffle_ps(ld1, ld1, 0xdd); // [Ai0,Ai1,Ai0,Ai1,...]
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^2
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex^2
  __attribute__((__target__("default"))) // default complex^2
  inline void pow2(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in[i]*in[i];
    return;
  }
  
  // 4th power real value
  #if !defined(DISABLE_AVX512) // AVX512 real^4
  __attribute__((__target__("avx512f")))
  inline void pow4(float const * const in, const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in[i]);
      ld1 = _mm512_mul_ps(ld1, ld1);
      ld1 = _mm512_mul_ps(ld1, ld1);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in[i]);
      ld1 = _mm512_mul_ps(ld1, ld1);
      ld1 = _mm512_mul_ps(ld1, ld1);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // end AVX512 real^4
  #if !defined(DISABLE_AVX) // AVX real^4
  __attribute__((__target__("avx")))
  inline void pow4(float const * const in, const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in[i]);
      ld1 = _mm256_mul_ps(ld1, ld1);
      ld1 = _mm256_mul_ps(ld1, ld1);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld1);
      ld1 = _mm256_mul_ps(ld1, ld1);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // end AVX real^4
  __attribute__((__target__("default"))) // default real^4
  inline void pow4(float const * const in, const int len, float * const out)
  {
    for(int i=0; i<len; ++i)
    {
      out[i] = in[i]*in[i];
      out[i] = out[i]*out[i];
    }
  }
  
  // 4th power complex value
  #if !defined(DISABLE_AVX512) // AVX512 complex^4
  __attribute__((__target__("avx512f")))
  inline void pow4(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, re, im;
    int i = 0;
    // note: 2 registers are split into separate real and imaginary registers
    // so the squaring needs no shuffles, then interleaved again for the store
    for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm512_unpacklo_ps(re, im); // A^4
      ld2 = _mm512_unpackhi_ps(re, im); // B^4
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld2);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem > 8) // one full and one partial register
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in[i+8])); // B
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm512_unpacklo_ps(re, im); // A^4
      ld2 = _mm512_unpackhi_ps(re, im); // B^4
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk, ld2);
    }
    else if(rem) // one partial register, split against itself
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in[i])); // A
      re = _mm512_shuffle_ps(ld1, ld1, 0x88); // [Ar0,Ar1,Ar0,Ar1,...]
      im = _mm512_shuffle_ps(ld1, ld1, 0xdd); // [Ai0,Ai1,Ai0,Ai1,...]
      // square input
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm512_unpacklo_ps(re, im); // A^4
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // end AVX512 complex^4
  #if !defined(DISABLE_AVX2) // AVX2 complex^4
  __attribute__((__target__("avx2,fma")))
  inline void pow4(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, re, im;
    int i = 0;
    // note: 2 registers are split into separate real and imaginary registers
    // so the squaring needs no shuffles, then interleaved again for the store
    for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^4
      ld2 = _mm256_unpackhi_ps(re, im); // B^4
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem > 4) // one full and one partial register
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), msk);
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^4
      ld2 = _mm256_unpackhi_ps(re, im); // B^4
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk, ld2);
    }
    else if(rem) // one partial register, split against itself
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      re = _mm256_shuffle_ps(ld1, ld1, 0x88); // [Ar0,Ar1,Ar0,Ar1,...]
      im = _mm256_shuffle_ps(ld1, ld1, 0xdd); // [Ai0,Ai1,Ai0,Ai1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^4
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex^4
  #if !defined(DISABLE_AVX) // AVX complex^4
  __attribute__((__target__("avx")))
  inline void pow4(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, re, im;
    int i = 0;
    // note: 2 registers are split into separate real and imaginary registers
    // so the squaring needs no shuffles, then interleaved again for the store
    for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^4
      ld2 = _mm256_unpackhi_ps(re, im); // B^4
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem > 4) // one full and one partial register
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), msk);
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^4
      ld2 = _mm256_unpackhi_ps(re, im); // B^4
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk, ld2);
    }
    else if(rem) // one partial register, split against itself
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      re = _mm256_shuffle_ps(ld1, ld1, 0x88); // [Ar0,Ar1,Ar0,Ar1,...]
      im = _mm256_shuffle_ps(ld1, ld1, 0xdd); // [Ai0,Ai1,Ai0,Ai1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^4
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex^4
  __attribute__((__target__("default"))) // default complex^4
  inline void pow4(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i)
    {
      out[i] = in[i]*in[i];
      out[i] = out[i]*out[i];
    }
    return;
  }


  // 8th power real value
  #if !defined(DISABLE_AVX512) // AVX512 real^8
  __attribute__((__target__("avx512f")))
  inline void pow8(float const * const in, const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in[i]);
      ld1 = _mm512_mul_ps(ld1, ld1);
      ld1 = _mm512_mul_ps(ld1, ld1);
      ld1 = _mm512_mul_ps(ld1, ld1);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in[i]);
      ld1 = _mm512_mul_ps(ld1, ld1);
      ld1 = _mm512_mul_ps(ld1, ld1);
      ld1 = _mm512_mul_ps(ld1, ld1);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // end AVX512 real^8
  #if !defined(DISABLE_AVX) // AVX real^8
  __attribute__((__target__("avx")))
  inline void pow8(float const * const in, const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in[i]);
      ld1 = _mm256_mul_ps(ld1, ld1);
      ld1 = _mm256_mul_ps(ld1, ld1);
      ld1 = _mm256_mul_ps(ld1, ld1);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld1);
      ld1 = _mm256_mul_ps(ld1, ld1);
      ld1 = _mm256_mul_ps(ld1, ld1);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // end AVX real^8
  __attribute__((__target__("default"))) // default real^8
  inline void pow8(float const * const in, const int len, float * const out)
  {
    for(int i=0; i<len; ++i)
    {
      out[i] = in[i]*in[i];
      out[i] = out[i]*out[i];
      out[i] = out[i]*out[i];
    }
  }
  

  // 8th power complex value
  #if !defined(DISABLE_AVX512) // AVX512 complex^8
  __attribute__((__target__("avx512f")))
  inline void pow8(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, re, im;
    int i = 0;
    // note: 2 registers are split into separate real and imaginary registers
    // so the squaring needs no shuffles, then interleaved again for the store
    for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 8th power
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm512_unpacklo_ps(re, im); // A^8
      ld2 = _mm512_unpackhi_ps(re, im); // B^8
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld2);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem > 8) // one full and one partial register
    {
      const __mmask16 mk = mask16((rem-8)<<1); // 2 floats per complex
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in[i+8])); // B
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 8th power
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm512_unpacklo_ps(re, im); // A^8
      ld2 = _mm512_unpackhi_ps(re, im); // B^8
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]), mk, ld2);
    }
    else if(rem) // one partial register, split against itself
    {
      const __mmask16 mk = mask16(rem<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const *>(&in[i])); // A
      re = _mm512_shuffle_ps(ld1, ld1, 0x88); // [Ar0,Ar1,Ar0,Ar1,...]
      im = _mm512_shuffle_ps(ld1, ld1, 0xdd); // [Ai0,Ai1,Ai0,Ai1,...]
      // square input
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 8th power
      ld2 = _mm512_mul_ps(im, im);      // im^2
      ld1 = _mm512_add_ps(re, re);      // 2*re
      im = _mm512_mul_ps(ld1, im);      // 2*re*im
      re = _mm512_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm512_unpacklo_ps(re, im); // A^8
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // end AVX512 complex^8
  #if !defined(DISABLE_AVX2) // AVX2 complex^8
  __attribute__((__target__("avx2,fma")))
  inline void pow8(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, re, im;
    int i = 0;
    // note: 2 registers are split into separate real and imaginary registers
    // so the squaring needs no shuffles, then interleaved again for the store
    for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 8th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^8
      ld2 = _mm256_unpackhi_ps(re, im); // B^8
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem > 4) // one full and one partial register
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), msk);
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 8th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^8
      ld2 = _mm256_unpackhi_ps(re, im); // B^8
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk, ld2);
    }
    else if(rem) // one partial register, split against itself
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      re = _mm256_shuffle_ps(ld1, ld1, 0x88); // [Ar0,Ar1,Ar0,Ar1,...]
      im = _mm256_shuffle_ps(ld1, ld1, 0xdd); // [Ai0,Ai1,Ai0,Ai1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      // square again to get 8th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_fmsub_ps(re, re, ld2);// re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^8
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex^8
  #if !defined(DISABLE_AVX) // AVX complex^8
  __attribute__((__target__("avx")))
  inline void pow8(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, re, im;
    int i = 0;
    // note: 2 registers are split into separate real and imaginary registers
    // so the squaring needs no shuffles, then interleaved again for the store
    for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      // square again to get 8th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^8
      ld2 = _mm256_unpackhi_ps(re, im); // B^8
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem > 4) // one full and one partial register
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem-4)<<1]));
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), msk);
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      // square again to get 8th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^8
      ld2 = _mm256_unpackhi_ps(re, im); // B^8
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]), msk, ld2);
    }
    else if(rem) // one partial register, split against itself
    {
      // 2 floats per complex, so double the remainder for the mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      re = _mm256_shuffle_ps(ld1, ld1, 0x88); // [Ar0,Ar1,Ar0,Ar1,...]
      im = _mm256_shuffle_ps(ld1, ld1, 0xdd); // [Ai0,Ai1,Ai0,Ai1,...]
      // square input
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      // square again to get 4th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      // square again to get 8th power
      ld2 = _mm256_mul_ps(im, im);      // im^2
      ld1 = _mm256_add_ps(re, re);      // 2*re
      im = _mm256_mul_ps(ld1, im);      // 2*re*im
      re = _mm256_mul_ps(re, re);       // re^2
      re = _mm256_sub_ps(re, ld2);      // re^2 - im^2
      ld1 = _mm256_unpacklo_ps(re, im); // A^8
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex^8
  __attribute__((__target__("default"))) // default complex^8
  inline void pow8(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    for(int i=0; i<len; ++i)
    {
      out[i] = in[i]*in[i];
      out[i] = out[i]*out[i];
      out[i] = out[i]*out[i];
    }
    return;
  }


  //=====================================================//
  // Magnitude/Magnitude Squared Function Implementation //
  //=====================================================//
  // magnitude squared of complex
  #if !defined(DISABLE_AVX512) // AVX512 |complex|^2
  __attribute__((__target__("avx512f")))
  inline void magSq(complex<float> const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, re;
    __m512i p1 = _mm512_setr_epi32(0,1,4,5,8,9,12,13,2,3,6,7,10,11,14,15);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once then
    // compact the float result into a single register
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A and B (8 complex elements each)
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm512_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm512_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm512_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm512_permutexvar_ps(p1, ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      }
      if(rem > 8) // one full and one partial register
      {
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm512_maskz_loadu_ps(mask16((rem-8)<<1), // 2 floats per complex
            reinterpret_cast<float const *>(&in[i+8])); // B
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm512_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm512_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm512_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm512_permutexvar_ps(p1, ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]),
            mask16(rem), ld1);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    ld1 = _mm512_maskz_loadu_ps(mask16(rem<<1),
        reinterpret_cast<float const *>(&in[i])); // A
    ld1 = _mm512_mul_ps(ld1, ld1);
    ld2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0^2,Ar0^2,...,Ai7^2,Ar7^2]
    ld1 = _mm512_add_ps(ld1, ld2);
    ld1 = _mm512_permutexvar_ps(_mm512_setr_epi32(0,2,4,6,8,10,12,14,1,3,5,
        7,9,11,13,15), ld1);
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mask16(rem), ld1);
    return;
  }
  #endif // end AVX512 |complex|^2
  #if !defined(DISABLE_AVX2) // AVX2 |complex|^2
  __attribute__((__target__("avx2,fma")))
  inline void magSq(complex<float> const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, re;
    __m256i p1 = _mm256_setr_epi32(0,1,4,5,2,3,6,7);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once then
    // compact the float result into a single register
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A and B (4 complex elements each)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm256_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm256_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm256_permutevar8x32_ps(ld1, p1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      }
      if(rem > 4) // one full and one partial register
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), // B
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm256_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm256_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm256_permutevar8x32_ps(ld1, p1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[rem])),
            ld1);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), // A
        _mm256_load_si256(reinterpret_cast<__m256i const *>(
        masks[rem<<1])));
    ld1 = _mm256_mul_ps(ld1, ld1);
    ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0^2,Ar0^2,...,Ai7^2,Ar7^2]
    ld1 = _mm256_add_ps(ld1, ld2);
    ld1 = _mm256_permutevar8x32_ps(ld1, _mm256_setr_epi32(0,2,4,6,1,3,5,7));
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
        _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[(rem)])),
        ld1);
    return;
  }
  #endif // end AVX2 |complex|^2
  #if !defined(DISABLE_AVX) // AVX |complex|^2
  __attribute__((__target__("avx")))
  inline void magSq(complex<float> const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, re;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once then
    // compact the float result into a single register
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A and B (4 complex elements each)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
        // regroup the 128-bit halves as [A0,A1,B0,B1] and [A2,A3,B2,B3] so
        // the shuffles that follow generate packed real and imaginary registers
        re = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [A0,A1,B0,B1]
        ld2 = _mm256_permute2f128_ps(ld1, ld2, 0x31); // [A2,A3,B2,B3]
        // compute |A|^2 and |B|^2 already compressed into a single register
        ld1 = _mm256_shuffle_ps(re,ld2,0x88); // [Ar0,Ar1,Ar2,Ar3,Br0,...]
        ld2 = _mm256_shuffle_ps(re,ld2,0xdd); // [Ai0,Ai1,Ai2,Ai3,Bi0,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Ai2^2,Ai3^2,Bi0^2,...]
        ld1 = _mm256_mul_ps(ld1,ld1);
        ld1 = _mm256_add_ps(ld1,ld2); // [|A0|^2,...,|A3|^2,|B0|^2,...]
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      }
      if(rem > 4) // one full and one partial register
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), // B
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
        // regroup as in the loop above
        re = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [A0,A1,B0,B1]
        ld2 = _mm256_permute2f128_ps(ld1, ld2, 0x31); // [A2,A3,B2,B3]
        // compute |A|^2 and |B|^2 already compressed into a single register
        ld1 = _mm256_shuffle_ps(re,ld2,0x88); // [Ar0,Ar1,Ar2,Ar3,Br0,...]
        ld2 = _mm256_shuffle_ps(re,ld2,0xdd); // [Ai0,Ai1,Ai2,Ai3,Bi0,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Ai2^2,Ai3^2,Bi0^2,...]
        ld1 = _mm256_mul_ps(ld1,ld1);
        ld1 = _mm256_add_ps(ld1,ld2); // [|A0|^2,...,|A3|^2,|B0|^2,...]
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[rem])),
            ld1);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), // A
        _mm256_load_si256(reinterpret_cast<__m256i const *>(
        masks[rem<<1])));
    ld1 = _mm256_mul_ps(ld1, ld1);
    ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0^2,Ar0^2,...,Ai7^2,Ar7^2]
    ld1 = _mm256_add_ps(ld1, ld2);
    // re-arrange pairs of output values to put all A before B
    ld2 = _mm256_permute2f128_ps (ld1, ld1, 0x01); // swap halves of ld1
    ld1 = _mm256_shuffle_ps(ld1, ld2, 0x88);
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
        _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[(rem)])),
        ld1);
    return;
  }
  #endif // end AVX |complex|^2
  __attribute__((__target__("default"))) // default |complex|^2
  inline void magSq(complex<float> const * const in, const int len,
      float * const out)
  {
    for(int i=0; i<len; ++i)
    {
      out[i] = (in[i].real()*in[i].real()+in[i].imag()*in[i].imag());
    }
  }

  // wrapper for magnitude squared of real data, must be defined after every
  // pow2 target version so the call goes through the runtime dispatcher
  inline void magSq(float const * const in, const int len, float * const out)
  {
    pow2(in, len, out);
  }

  // absolute value of real
  #if !defined(DISABLE_AVX512) // AVX512 |real|
  __attribute__((__target__("avx512f")))
  inline void mag(float const * const in, const int len, float * const out)
  {
    if(len<=0) return;
    __m512 ld1;
    // constant to make the sign bit positive
    const __m512i pos = _mm512_set1_epi32(0x7fffffff);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in[i]);
      // force sign bit positive, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_and_ps (casting should be free)
      ld1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(ld1),pos));
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in[i]);
      // force sign bit positive, note the casting is done to avoid requiring
      // the AVX512DQ instruction _mm512_and_ps (casting should be free)
      ld1 = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(ld1),pos));
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // end AVX512 |real|
  #if !defined(DISABLE_AVX) // AVX |real|
  __attribute__((__target__("avx")))
  inline void mag(float const * const in, const int len, float * const out)
  {
    if(len<=0) return;
    __m256 ld1;
    // constant to make the sign bit positive
    const __m256 pos = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in[i]);
      ld1 = _mm256_and_ps(ld1, pos);
      _mm256_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in[i], msk);
      ld1 = _mm256_and_ps(ld1, pos);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // end AVX |real|
  __attribute__((__target__("default"))) // default |real|
  inline void mag(float const * const in, const int len, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = std::abs(in[i]);
    return;
  }
  // wrapper for magnitude of real data, must be defined after every mag
  // target version so the call goes through the runtime dispatcher
  inline void abs(float const * const in, const int len, float * const out)
  {
    mag(in, len, out);
  }
  
  // magnitude of complex (true square root)
  #if !defined(DISABLE_AVX512) // AVX512 |complex|
  __attribute__((__target__("avx512f")))
  inline void mag(complex<float> const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, re;
    __m512i p1 = _mm512_setr_epi32(0,1,4,5,8,9,12,13,2,3,6,7,10,11,14,15);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once then
    // compact the float result into a single register
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A and B (8 complex elements each)
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm512_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm512_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm512_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        ld1 = _mm512_sqrt_ps(ld1);    // [|A0|,|A1|,|B0|,|B1|,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm512_permutexvar_ps(p1, ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      }
      if(rem > 8) // one full and one partial register
      {
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm512_maskz_loadu_ps(mask16((rem-8)<<1), // 2 floats per complex
            reinterpret_cast<float const *>(&in[i+8])); // B
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm512_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm512_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm512_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        ld1 = _mm512_sqrt_ps(ld1);    // [|A0|,|A1|,|B0|,|B1|,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm512_permutexvar_ps(p1, ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]),
            mask16(rem), ld1);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    ld1 = _mm512_maskz_loadu_ps(mask16(rem<<1),
        reinterpret_cast<float const *>(&in[i])); // A
    ld1 = _mm512_mul_ps(ld1, ld1);
    ld2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0^2,Ar0^2,...,Ai7^2,Ar7^2]
    ld1 = _mm512_add_ps(ld1, ld2);
    ld1 = _mm512_sqrt_ps(ld1);
    ld1 = _mm512_permutexvar_ps(_mm512_setr_epi32(0,2,4,6,8,10,12,14,1,3,5,
        7,9,11,13,15), ld1);
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mask16(rem), ld1);
    return;
  }
  #endif // end AVX512 |complex|
  #if !defined(DISABLE_AVX2) // AVX2 |complex|
  __attribute__((__target__("avx2,fma")))
  inline void mag(complex<float> const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, re;
    __m256i p1 = _mm256_setr_epi32(0,1,4,5,2,3,6,7);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once then
    // compact the float result into a single register
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A and B (4 complex elements each)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm256_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm256_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        ld1 = _mm256_sqrt_ps(ld1);    // [|A0|,|A1|,|B0|,|B1|,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm256_permutevar8x32_ps(ld1, p1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      }
      if(rem > 4) // one full and one partial register
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), // B
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm256_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm256_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        ld1 = _mm256_sqrt_ps(ld1);    // [|A0|,|A1|,|B0|,|B1|,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm256_permutevar8x32_ps(ld1, p1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[rem])),
            ld1);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), // A
        _mm256_load_si256(reinterpret_cast<__m256i const *>(
        masks[rem<<1])));
    ld1 = _mm256_mul_ps(ld1, ld1);
    ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0^2,Ar0^2,...,Ai7^2,Ar7^2]
    ld1 = _mm256_add_ps(ld1, ld2);
    ld1 = _mm256_sqrt_ps(ld1);
    ld1 = _mm256_permutevar8x32_ps(ld1, _mm256_setr_epi32(0,2,4,6,1,3,5,7));
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
        _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[(rem)])),
        ld1);
    return;
  }
  #endif // end AVX2 |complex|
  #if !defined(DISABLE_AVX) // AVX |complex|
  __attribute__((__target__("avx")))
  inline void mag(complex<float> const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, re;
    int i = 0;
    // note: we process 2 registers worth of complex elements at once then
    // compact the float result into a single register
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A and B (4 complex elements each)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
        // regroup the 128-bit halves as [A0,A1,B0,B1] and [A2,A3,B2,B3] so
        // the shuffles that follow generate packed real and imaginary registers
        re = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [A0,A1,B0,B1]
        ld2 = _mm256_permute2f128_ps(ld1, ld2, 0x31); // [A2,A3,B2,B3]
        // compute |A|^2 and |B|^2 already compressed into a single register
        ld1 = _mm256_shuffle_ps(re,ld2,0x88); // [Ar0,Ar1,Ar2,Ar3,Br0,...]
        ld2 = _mm256_shuffle_ps(re,ld2,0xdd); // [Ai0,Ai1,Ai2,Ai3,Bi0,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Ai2^2,Ai3^2,Bi0^2,...]
        ld1 = _mm256_mul_ps(ld1,ld1);
        ld1 = _mm256_add_ps(ld1,ld2); // [|A0|^2,...,|A3|^2,|B0|^2,...]
        ld1 = _mm256_sqrt_ps(ld1);    // [|A0|,...,|A3|,|B0|,...]
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      }
      if(rem > 4) // one full and one partial register
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), // B
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
        // regroup as in the loop above
        re = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [A0,A1,B0,B1]
        ld2 = _mm256_permute2f128_ps(ld1, ld2, 0x31); // [A2,A3,B2,B3]
        // compute |A|^2 and |B|^2 already compressed into a single register
        ld1 = _mm256_shuffle_ps(re,ld2,0x88); // [Ar0,Ar1,Ar2,Ar3,Br0,...]
        ld2 = _mm256_shuffle_ps(re,ld2,0xdd); // [Ai0,Ai1,Ai2,Ai3,Bi0,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Ai2^2,Ai3^2,Bi0^2,...]
        ld1 = _mm256_mul_ps(ld1,ld1);
        ld1 = _mm256_add_ps(ld1,ld2); // [|A0|^2,...,|A3|^2,|B0|^2,...]
        ld1 = _mm256_sqrt_ps(ld1);    // [|A0|,...,|A3|,|B0|,...]
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[rem])),
            ld1);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), // A
        _mm256_load_si256(reinterpret_cast<__m256i const *>(
        masks[rem<<1])));
    ld1 = _mm256_mul_ps(ld1, ld1);
    ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0^2,Ar0^2,...,Ai7^2,Ar7^2]
    ld1 = _mm256_add_ps(ld1, ld2);
    ld1 = _mm256_sqrt_ps(ld1);
    // re-arrange pairs of output values to put all A before B
    ld2 = _mm256_permute2f128_ps (ld1, ld1, 0x01); // swap halves of ld1
    ld1 = _mm256_shuffle_ps(ld1, ld2, 0x88);
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
        _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[(rem)])),
        ld1);
    return;
  }
  #endif // end AVX |complex|
  __attribute__((__target__("default"))) // default |complex|
  inline void mag(complex<float> const * const in, const int len,
      float * const out)
  {
    for(int i=0; i<len; ++i)
    {
      out[i] = sqrtf(in[i].real()*in[i].real()+in[i].imag()*in[i].imag());
    }
  }
  
  // magnitude of complex (reciprocal square root approximation)
  #if !defined(DISABLE_AVX512) // AVX512 |complex| (rsqrt)
  __attribute__((__target__("avx512f")))
  inline void magr(complex<float> const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    const int rem = len&15; // remaining elements (note len&15 == len%16)
    __m512 ld1, ld2, re;
    // the rsqrt is zero-masked to the lanes where |in|^2 != 0, so |in| = 0
    // gives 0*0 = 0 instead of 0*rsqrt(0) = 0*inf = NaN
    __mmask16 nz;
    const __m512 zero = _mm512_setzero_ps();
    __m512i p1 = _mm512_setr_epi32(0,1,4,5,8,9,12,13,2,3,6,7,10,11,14,15);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once then
    // compact the float result into a single register
    if(len > 8)
    {
      for(; i<len-15; i+=16) // process 2 registers (16 complex elements)
      {
        // load A and B (8 complex elements each)
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm512_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm512_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm512_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        nz = _mm512_cmp_ps_mask(ld1, zero, _CMP_NEQ_UQ); // |in|^2 != 0
        ld2 = _mm512_maskz_rsqrt14_ps(nz, ld1); // [1/|A0|,1/|A1|,1/|B0|,...]
        ld1 = _mm512_mul_ps(ld1,ld2);     // [|A0|,|A1|,|B0|,|B1|,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm512_permutexvar_ps(p1, ld1);
        _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      }
      if(rem > 8) // one full and one partial register
      {
        ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm512_maskz_loadu_ps(mask16((rem-8)<<1), // 2 floats per complex
            reinterpret_cast<float const *>(&in[i+8])); // B
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm512_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm512_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm512_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm512_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        nz = _mm512_cmp_ps_mask(ld1, zero, _CMP_NEQ_UQ); // |in|^2 != 0
        ld2 = _mm512_maskz_rsqrt14_ps(nz, ld1); // [1/|A0|,1/|A1|,1/|B0|,...]
        ld1 = _mm512_mul_ps(ld1,ld2);     // [|A0|,|A1|,|B0|,|B1|,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm512_permutexvar_ps(p1, ld1);
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]),
            mask16(rem), ld1);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    ld1 = _mm512_maskz_loadu_ps(mask16(rem<<1),
        reinterpret_cast<float const *>(&in[i])); // A
    ld1 = _mm512_mul_ps(ld1, ld1);
    ld2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0^2,Ar0^2,...,Ai7^2,Ar7^2]
    ld1 = _mm512_add_ps(ld1, ld2);
    nz = _mm512_cmp_ps_mask(ld1, zero, _CMP_NEQ_UQ); // |in|^2 != 0
    ld2 = _mm512_maskz_rsqrt14_ps(nz, ld1);
    ld1 = _mm512_mul_ps(ld1,ld2);
    ld1 = _mm512_permutexvar_ps(_mm512_setr_epi32(0,2,4,6,8,10,12,14,1,3,5,
        7,9,11,13,15), ld1);
    _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mask16(rem), ld1);
    return;
  }
  #endif // end AVX512 |complex| (rsqrt)
  #if !defined(DISABLE_AVX2) // AVX2 |complex| (rsqrt)
  __attribute__((__target__("avx2,fma")))
  inline void magr(complex<float> const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, re;
    // |in| is computed as rcp(rsqrt(|in|^2)) rather than |in|^2*rsqrt(|in|^2)
    // so |in| = 0 gives rcp(inf) = 0 instead of 0*inf = NaN, at the same cost
    __m256i p1 = _mm256_setr_epi32(0,1,4,5,2,3,6,7);
    int i = 0;
    // note: we process 2 registers worth of complex elements at once then
    // compact the float result into a single register
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A and B (4 complex elements each)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm256_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm256_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        ld2 = _mm256_rsqrt_ps(ld1);   // [1/|A0|,1/|A1|,1/|B0|,1/|B1|,...]
        ld1 = _mm256_rcp_ps(ld2);     // [|A0|,|A1|,|B0|,|B1|,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm256_permutevar8x32_ps(ld1, p1);
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      }
      if(rem > 4) // one full and one partial register
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), // B
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
        // compute |A|^2 and |B|^2 and compress into 1 register
        re = _mm256_shuffle_ps(ld1,ld2,0x88); // [Ar0,Ar1,Br0,Br1,...]
        ld2 = _mm256_shuffle_ps(ld1,ld2,0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Bi0^2,Bi1^2,...]
        ld1 = _mm256_fmadd_ps(re,re,ld2); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
        ld2 = _mm256_rsqrt_ps(ld1);   // [1/|A0|,1/|A1|,1/|B0|,1/|B1|,...]
        ld1 = _mm256_rcp_ps(ld2);     // [|A0|,|A1|,|B0|,|B1|,...]
        // re-arrange pairs of output values to put all A before B
        ld1 = _mm256_permutevar8x32_ps(ld1, p1);
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[rem])),
            ld1);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), // A
        _mm256_load_si256(reinterpret_cast<__m256i const *>(
        masks[rem<<1])));
    ld1 = _mm256_mul_ps(ld1, ld1);
    ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0^2,Ar0^2,...,Ai7^2,Ar7^2]
    ld1 = _mm256_add_ps(ld1, ld2);
    ld2 = _mm256_rsqrt_ps(ld1);
    ld1 = _mm256_rcp_ps(ld2);
    ld1 = _mm256_permutevar8x32_ps(ld1, _mm256_setr_epi32(0,2,4,6,1,3,5,7));
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
        _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[(rem)])),
        ld1);
    return;
  }
  #endif // end AVX2 |complex| (rsqrt)
  #if !defined(DISABLE_AVX) // AVX |complex| (rsqrt)
  __attribute__((__target__("avx")))
  inline void magr(complex<float> const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    const int rem = len&7; // remaining elements (note len&7 == len%8)
    __m256 ld1, ld2, re;
    // |in| is computed as rcp(rsqrt(|in|^2)) rather than |in|^2*rsqrt(|in|^2)
    // so |in| = 0 gives rcp(inf) = 0 instead of 0*inf = NaN, at the same cost
    int i = 0;
    // note: we process 2 registers worth of complex elements at once then
    // compact the float result into a single register
    if(len > 4)
    {
      for(; i<len-7; i+=8) // process 2 registers (8 complex elements)
      {
        // load A and B (4 complex elements each)
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
        // regroup the 128-bit halves as [A0,A1,B0,B1] and [A2,A3,B2,B3] so
        // the shuffles that follow generate packed real and imaginary registers
        re = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [A0,A1,B0,B1]
        ld2 = _mm256_permute2f128_ps(ld1, ld2, 0x31); // [A2,A3,B2,B3]
        // compute |A|^2 and |B|^2 already compressed into a single register
        ld1 = _mm256_shuffle_ps(re,ld2,0x88); // [Ar0,Ar1,Ar2,Ar3,Br0,...]
        ld2 = _mm256_shuffle_ps(re,ld2,0xdd); // [Ai0,Ai1,Ai2,Ai3,Bi0,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Ai2^2,Ai3^2,Bi0^2,...]
        ld1 = _mm256_mul_ps(ld1,ld1);
        ld1 = _mm256_add_ps(ld1,ld2); // [|A0|^2,...,|A3|^2,|B0|^2,...]
        ld2 = _mm256_rsqrt_ps(ld1);   // [1/|A0|,...,1/|A3|,1/|B0|,...,1/|B3|]
        ld1 = _mm256_rcp_ps(ld2);     // [|A0|,...,|A3|,|B0|,...,|B3|]
        _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      }
      if(rem > 4) // one full and one partial register
      {
        ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]), // B
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
        // regroup as in the loop above
        re = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [A0,A1,B0,B1]
        ld2 = _mm256_permute2f128_ps(ld1, ld2, 0x31); // [A2,A3,B2,B3]
        // compute |A|^2 and |B|^2 already compressed into a single register
        ld1 = _mm256_shuffle_ps(re,ld2,0x88); // [Ar0,Ar1,Ar2,Ar3,Br0,...]
        ld2 = _mm256_shuffle_ps(re,ld2,0xdd); // [Ai0,Ai1,Ai2,Ai3,Bi0,...]
        ld2 = _mm256_mul_ps(ld2,ld2); // [Ai0^2,Ai1^2,Ai2^2,Ai3^2,Bi0^2,...]
        ld1 = _mm256_mul_ps(ld1,ld1);
        ld1 = _mm256_add_ps(ld1,ld2); // [|A0|^2,...,|A3|^2,|B0|^2,...]
        ld2 = _mm256_rsqrt_ps(ld1);   // [1/|A0|,...,1/|A3|,1/|B0|,...,1/|B3|]
        ld1 = _mm256_rcp_ps(ld2);     // [|A0|,...,|A3|,|B0|,...,|B3|]
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[rem])),
            ld1);
        return;
      }
      else if(!rem) return;
    }
    // if we got here, there is one partial register left to complete
    // each complex is 2 floats, so double rem
    ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), // A
        _mm256_load_si256(reinterpret_cast<__m256i const *>(
        masks[rem<<1])));
    ld1 = _mm256_mul_ps(ld1, ld1);
    ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0^2,Ar0^2,...,Ai7^2,Ar7^2]
    ld1 = _mm256_add_ps(ld1, ld2);
    ld2 = _mm256_rsqrt_ps(ld1);
    ld1 = _mm256_rcp_ps(ld2);
    // re-arrange pairs of output values to put all A before B
    ld2 = _mm256_permute2f128_ps (ld1, ld1, 0x01); // swap halves of ld1
    ld1 = _mm256_shuffle_ps(ld1, ld2, 0x88);
    _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]),
        _mm256_load_si256(reinterpret_cast<__m256i const *>(masks[(rem)])),
        ld1);
    return;
  }
  #endif // end AVX |complex| (rsqrt)
  __attribute__((__target__("default"))) // default |complex| (rsqrt)
  inline void magr(complex<float> const * const in, const int len,
      float * const out)
  {
    for(int i=0; i<len; ++i)
    {
      out[i] = sqrtf(in[i].real()*in[i].real()+in[i].imag()*in[i].imag());
    }
  }
  
  //============================================================//
  // Unit/Coordinate System Conversion Function Implementations //
  //============================================================//
  // 8-bit int to float conversion
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f,avx512bw,avx512vl")))
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    // note: the output is 4x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    __m512 cv;
    int i = len&~15; // start of the final partial register
    // handle remaining elements first (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi8_epi32(_mm_maskz_loadu_epi8(mk,
          &in[i])));
      _mm512_mask_storeu_ps(&out[i], mk, cv);
    }
    for(i-=16; i>=0; i-=16) // process 16 real elements per register
    {
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi8_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i]))));
      _mm512_storeu_ps(&out[i], cv);
    }
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX2) // AVX2
  __attribute__((__target__("avx2")))
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    // note: the output is 4x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    if(len < 8) // less than 1 register, convert one element at a time
    {
      for(int i=len-1; i>=0; --i) out[i] = static_cast<float>(in[i]);
      return;
    }
    __m256 cv;
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8), maskload can't
    // load 8-bit integers, instead load the last 8 elements (overlapping the
    // previous register) and only store the last rem
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(rmasks[rem])); // last rem floats
      cv = _mm256_cvtepi32_ps(_mm256_cvtepi8_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[len-8]))));
      _mm256_maskstore_ps(&out[len-8], msk, cv);
    }
    for(i-=8; i>=0; i-=8) // process 8 real elements per register
    {
      cv = _mm256_cvtepi32_ps(_mm256_cvtepi8_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i]))));
      _mm256_storeu_ps(&out[i], cv);
    }
    return;
  }
  #endif // end AVX2
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    // note: the output is 4x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    if(len < 8) // less than 1 register, convert one element at a time
    {
      for(int i=len-1; i>=0; --i) out[i] = static_cast<float>(in[i]);
      return;
    }
    // note: AVX has no 256-bit integer sign extension, so each half of the
    // register is extended separately with SSE4.1
    __m128i lo, hi;
    __m256 cv;
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8), maskload can't
    // load 8-bit integers, instead load the last 8 elements (overlapping the
    // previous register) and only store the last rem
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(rmasks[rem])); // last rem floats
      lo = _mm_cvtepi8_epi32(_mm_loadu_si32(&in[len-8]));
      hi = _mm_cvtepi8_epi32(_mm_loadu_si32(&in[len-4]));
      cv = _mm256_cvtepi32_ps(_mm256_insertf128_si256(
          _mm256_castsi128_si256(lo), hi, 1));
      _mm256_maskstore_ps(&out[len-8], msk, cv);
    }
    for(i-=8; i>=0; i-=8) // process 8 real elements per register
    {
      lo = _mm_cvtepi8_epi32(_mm_loadu_si32(&in[i]));   // [0,1,2,3]
      hi = _mm_cvtepi8_epi32(_mm_loadu_si32(&in[i+4])); // [4,5,6,7]
      cv = _mm256_cvtepi32_ps(_mm256_insertf128_si256(
          _mm256_castsi128_si256(lo), hi, 1));
      _mm256_storeu_ps(&out[i], cv);
    }
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out)
  {
    // note: going from the end to the start so this also works in place
    for(int i=len-1; i>=0; --i) out[i] = static_cast<float>(in[i]);
    return;
  }
  inline void cnvrtInt8ToFloat(void * const inout, const int len)
  {
    cnvrtInt8ToFloat(reinterpret_cast<int8_t const *>(inout), len,
        reinterpret_cast<float *>(inout));
  }
  
  // 8-bit integer to float conversion with min/max output
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f,avx512bw,avx512vl")))
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    // note: the output is 4x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    // note: the min/max is found on the 8-bit integers since converting
    // them to float doesn't change their order, that way each min and max
    // covers 64 elements (4 registers of output)
    __m512i ld;
    __m512i mn = _mm512_set1_epi8(in[0]);
    __m512i mx = mn;
    __m512 cv;
    int i = len&~63; // start of the remaining elements
    // handle remaining elements first (note len&63 == len%64), the min/max
    // is found on all of them at once, then they are converted one register
    // at a time
    const int rem = len&63;
    if(rem)
    {
      const __mmask64 mk = 0xffffffffffffffffULL >> (64-rem); // first rem
      ld = _mm512_maskz_loadu_epi8(mk, &in[i]);
      // only compare the loaded elements, the others keep their min/max
      mn = _mm512_mask_min_epi8(mn, mk, ld, mn);
      mx = _mm512_mask_max_epi8(mx, mk, ld, mx);
      // the last (partial) register first so that this works in place
      int j = i + (rem&~15); // start of the final partial register
      if(rem&15)
      {
        const __mmask16 mk16 = mask16(rem&15);
        cv = _mm512_cvtepi32_ps(_mm512_cvtepi8_epi32(_mm_maskz_loadu_epi8(mk16,
            &in[j])));
        _mm512_mask_storeu_ps(&out[j], mk16, cv);
      }
      for(j-=16; j>=i; j-=16) // full registers
      {
        cv = _mm512_cvtepi32_ps(_mm512_cvtepi8_epi32(_mm_loadu_si128(
            reinterpret_cast<__m128i const *>(&in[j]))));
        _mm512_storeu_ps(&out[j], cv);
      }
    }
    for(i-=64; i>=0; i-=64) // process 64 real elements (4 registers)
    {
      ld = _mm512_loadu_si512(&in[i]);
      mn = _mm512_min_epi8(ld, mn);
      mx = _mm512_max_epi8(ld, mx);
      // convert the last register first so that this works in place (the
      // integers are loaded again, which is faster than extracting them)
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi8_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i+48]))));
      _mm512_storeu_ps(&out[i+48], cv);
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi8_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i+32]))));
      _mm512_storeu_ps(&out[i+32], cv);
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi8_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i+16]))));
      _mm512_storeu_ps(&out[i+16], cv);
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi8_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i]))));
      _mm512_storeu_ps(&out[i], cv);
    }
    // sign extend the byte min/max 16 at a time to reduce them to 16
    // integers, then the horizontal min and max
    __m512i r = _mm512_min_epi32(_mm512_cvtepi8_epi32(
        _mm512_castsi512_si128(mn)),
        _mm512_cvtepi8_epi32(_mm512_extracti32x4_epi32(mn, 1)));
    r = _mm512_min_epi32(r, _mm512_cvtepi8_epi32(
        _mm512_extracti32x4_epi32(mn, 2)));
    r = _mm512_min_epi32(r, _mm512_cvtepi8_epi32(
        _mm512_extracti32x4_epi32(mn, 3)));
    const float fmn = static_cast<float>(_mm512_reduce_min_epi32(r));
    r = _mm512_max_epi32(_mm512_cvtepi8_epi32(
        _mm512_castsi512_si128(mx)),
        _mm512_cvtepi8_epi32(_mm512_extracti32x4_epi32(mx, 1)));
    r = _mm512_max_epi32(r, _mm512_cvtepi8_epi32(
        _mm512_extracti32x4_epi32(mx, 2)));
    r = _mm512_max_epi32(r, _mm512_cvtepi8_epi32(
        _mm512_extracti32x4_epi32(mx, 3)));
    const float fmx = static_cast<float>(_mm512_reduce_max_epi32(r));
    minOut = noInit? std::min(minOut, fmn) : fmn;
    maxOut = noInit? std::max(maxOut, fmx) : fmx;
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX2) // AVX2
  __attribute__((__target__("avx2")))
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    // note: the output is 4x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    if(len < 8) // less than 1 register, convert one element at a time
    {
      if(!noInit)
      {
        minOut = static_cast<float>(in[0]);
        maxOut = static_cast<float>(in[0]);
      }
      for(int i=len-1; i>=0; --i)
      {
        out[i] = static_cast<float>(in[i]);
        minOut = std::min(minOut, out[i]);
        maxOut = std::max(maxOut, out[i]);
      }
      return;
    }
    // note: the min/max is found on the 8-bit integers since converting
    // them to float doesn't change their order, that way each min and max
    // covers 32 elements (4 registers of output)
    __m256i ld;
    __m256i mn = _mm256_set1_epi8(in[0]);
    __m256i mx = mn;
    __m256 cv, cva, cvb, cvc, cvd;
    __m128i hs;
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8), maskload can't
    // load 8-bit integers, instead load the last 8 elements (overlapping the
    // previous register) and only store the last rem
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(rmasks[rem])); // last rem floats
      // for the min/max the 8 loaded bytes are repeated to fill the register
      ld = _mm256_broadcastq_epi64(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[len-8])));
      // the overlapping elements are compared twice, which doesn't matter
      mn = _mm256_min_epi8(ld, mn);
      mx = _mm256_max_epi8(ld, mx);
      cv = _mm256_cvtepi32_ps(_mm256_cvtepi8_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[len-8]))));
      _mm256_maskstore_ps(&out[len-8], msk, cv);
    }
    // then single registers until the rest is a multiple of 32 elements
    for(i-=8; i>=(len&~31); i-=8)
    {
      ld = _mm256_broadcastq_epi64(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i])));
      mn = _mm256_min_epi8(ld, mn);
      mx = _mm256_max_epi8(ld, mx);
      cv = _mm256_cvtepi32_ps(_mm256_cvtepi8_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i]))));
      _mm256_storeu_ps(&out[i], cv);
    }
    for(i=(len&~31)-32; i>=0; i-=32) // process 32 real elements (4 registers)
    {
      ld = _mm256_loadu_si256(reinterpret_cast<__m256i const *>(&in[i]));
      mn = _mm256_min_epi8(ld, mn);
      mx = _mm256_max_epi8(ld, mx);
      // load and convert every register before storing any of them so that
      // this works in place (the integers are loaded again, which is faster
      // than extracting them)
      cvd = _mm256_cvtepi32_ps(_mm256_cvtepi8_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i+24]))));
      cvc = _mm256_cvtepi32_ps(_mm256_cvtepi8_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i+16]))));
      cvb = _mm256_cvtepi32_ps(_mm256_cvtepi8_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i+8]))));
      cva = _mm256_cvtepi32_ps(_mm256_cvtepi8_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i]))));
      _mm256_storeu_ps(&out[i+24], cvd);
      _mm256_storeu_ps(&out[i+16], cvc);
      _mm256_storeu_ps(&out[i+8], cvb);
      _mm256_storeu_ps(&out[i], cva);
    }
    // sign extend the byte min/max 8 at a time to reduce them to 8
    // integers, then the horizontal min and max
    hs = _mm256_castsi256_si128(mn);
    ld = _mm256_min_epi32(_mm256_cvtepi8_epi32(hs),
        _mm256_cvtepi8_epi32(_mm_srli_si128(hs, 8)));
    hs = _mm256_extracti128_si256(mn, 1);
    ld = _mm256_min_epi32(ld, _mm256_cvtepi8_epi32(hs));
    mn = _mm256_min_epi32(ld, _mm256_cvtepi8_epi32(_mm_srli_si128(hs, 8)));
    hs = _mm256_castsi256_si128(mx);
    ld = _mm256_max_epi32(_mm256_cvtepi8_epi32(hs),
        _mm256_cvtepi8_epi32(_mm_srli_si128(hs, 8)));
    hs = _mm256_extracti128_si256(mx, 1);
    ld = _mm256_max_epi32(ld, _mm256_cvtepi8_epi32(hs));
    mx = _mm256_max_epi32(ld, _mm256_cvtepi8_epi32(_mm_srli_si128(hs, 8)));
    hs = _mm_min_epi32(_mm256_castsi256_si128(mn),
        _mm256_extracti128_si256(mn, 1)); // [min(0,4),min(1,5),...]
    hs = _mm_min_epi32(hs, _mm_shuffle_epi32(hs, 0x4e)); // [min(0,4,2,6),...]
    hs = _mm_min_epi32(hs, _mm_shuffle_epi32(hs, 0xb1)); // [min(0,...,7),...]
    const float fmn = static_cast<float>(_mm_cvtsi128_si32(hs));
    hs = _mm_max_epi32(_mm256_castsi256_si128(mx),
        _mm256_extracti128_si256(mx, 1)); // [max(0,4),max(1,5),...]
    hs = _mm_max_epi32(hs, _mm_shuffle_epi32(hs, 0x4e)); // [max(0,4,2,6),...]
    hs = _mm_max_epi32(hs, _mm_shuffle_epi32(hs, 0xb1)); // [max(0,...,7),...]
    const float fmx = static_cast<float>(_mm_cvtsi128_si32(hs));
    minOut = noInit? std::min(minOut, fmn) : fmn;
    maxOut = noInit? std::max(maxOut, fmx) : fmx;
    return;
  }
  #endif // end AVX2
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    // note: the output is 4x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    if(len < 8) // less than 1 register, convert one element at a time
    {
      if(!noInit)
      {
        minOut = static_cast<float>(in[0]);
        maxOut = static_cast<float>(in[0]);
      }
      for(int i=len-1; i>=0; --i)
      {
        out[i] = static_cast<float>(in[i]);
        minOut = std::min(minOut, out[i]);
        maxOut = std::max(maxOut, out[i]);
      }
      return;
    }
    // note: AVX has no 256-bit integer sign extension or min/max, so each
    // half of the register is extended and compared separately with SSE4.1
    // note: the min/max is found on the integers since converting them to
    // float doesn't change their order, integer min/max has a 1 cycle latency
    // (vs. 4 for float) so one min and max accumulator doesn't stall the loop
    __m128i lo, hi;
    __m128i mn = _mm_set1_epi32(in[0]);
    __m128i mx = mn;
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8), maskload can't
    // load 8-bit integers, instead load the last 8 elements (overlapping the
    // previous register) and only store the last rem
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(rmasks[rem])); // last rem floats
      lo = _mm_cvtepi8_epi32(_mm_loadu_si32(&in[len-8]));
      hi = _mm_cvtepi8_epi32(_mm_loadu_si32(&in[len-4]));
      // the overlapping elements are compared twice, which doesn't matter
      mn = _mm_min_epi32(mn, _mm_min_epi32(lo, hi));
      mx = _mm_max_epi32(mx, _mm_max_epi32(lo, hi));
      _mm256_maskstore_ps(&out[len-8], msk, _mm256_cvtepi32_ps(
          _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1)));
    }
    for(i-=8; i>=0; i-=8) // process 8 real elements per register
    {
      lo = _mm_cvtepi8_epi32(_mm_loadu_si32(&in[i]));   // [0,1,2,3]
      hi = _mm_cvtepi8_epi32(_mm_loadu_si32(&in[i+4])); // [4,5,6,7]
      mn = _mm_min_epi32(mn, _mm_min_epi32(lo, hi));
      mx = _mm_max_epi32(mx, _mm_max_epi32(lo, hi));
      _mm256_storeu_ps(&out[i], _mm256_cvtepi32_ps(
          _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1)));
    }
    // horizontal min and max
    mn = _mm_min_epi32(mn, _mm_shuffle_epi32(mn, 0x4e)); // [min(0,2),...]
    mn = _mm_min_epi32(mn, _mm_shuffle_epi32(mn, 0xb1)); // [min(0,...,3),...]
    mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0x4e)); // [max(0,2),...]
    mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0xb1)); // [max(0,...,3),...]
    const float fmn = static_cast<float>(_mm_cvtsi128_si32(mn));
    const float fmx = static_cast<float>(_mm_cvtsi128_si32(mx));
    minOut = noInit? std::min(minOut, fmn) : fmn;
    maxOut = noInit? std::max(maxOut, fmx) : fmx;
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void cnvrtInt8ToFloat(int8_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    if(!noInit)
    {
      minOut = static_cast<float>(in[0]);
      maxOut = static_cast<float>(in[0]);
    }
    // note: going from the end to the start so this also works in place
    for(int i=len-1; i>=0; --i)
    {
      out[i] = static_cast<float>(in[i]);
      minOut = std::min(minOut, out[i]);
      maxOut = std::max(maxOut, out[i]);
    }
    return;
  }


  // 16-bit int to float conversion
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f,avx512bw,avx512vl")))
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    // note: the output is 2x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    __m512 cv;
    int i = len&~15; // start of the final partial register
    // handle remaining elements first (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi16_epi32(_mm256_maskz_loadu_epi16(
          mk, &in[i])));
      _mm512_mask_storeu_ps(&out[i], mk, cv);
    }
    for(i-=16; i>=0; i-=16) // process 16 real elements per register
    {
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi16_epi32(_mm256_loadu_si256(
          reinterpret_cast<__m256i const *>(&in[i]))));
      _mm512_storeu_ps(&out[i], cv);
    }
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX2) // AVX2
  __attribute__((__target__("avx2")))
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    // note: the output is 2x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    if(len < 8) // less than 1 register, convert one element at a time
    {
      // note: the input is read with memcpy, in place the int16_t input and
      // float output alias and the compiler could otherwise assume they don't
      // and reorder the loop, the memcpy is actually free as it compiles out
      int16_t tmp;
      for(int i=len-1; i>=0; --i)
      {
        std::memcpy(&tmp, &in[i], sizeof(tmp));
        out[i] = static_cast<float>(tmp);
      }
      return;
    }
    __m256 cv;
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8), maskload can't
    // load 16-bit integers, instead load the last 8 elements (overlapping the
    // previous register) and only store the last rem
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(rmasks[rem])); // last rem floats
      cv = _mm256_cvtepi32_ps(_mm256_cvtepi16_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[len-8]))));
      _mm256_maskstore_ps(&out[len-8], msk, cv);
    }
    for(i-=8; i>=0; i-=8) // process 8 real elements per register
    {
      cv = _mm256_cvtepi32_ps(_mm256_cvtepi16_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i]))));
      _mm256_storeu_ps(&out[i], cv);
    }
    return;
  }
  #endif // end AVX2
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    // note: the output is 2x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    if(len < 8) // less than 1 register, convert one element at a time
    {
      // note: the input is read with memcpy, in place the int16_t input and
      // float output alias and the compiler could otherwise assume they don't
      // and reorder the loop, the memcpy is actually free as it compiles out
      int16_t tmp;
      for(int i=len-1; i>=0; --i)
      {
        std::memcpy(&tmp, &in[i], sizeof(tmp));
        out[i] = static_cast<float>(tmp);
      }
      return;
    }
    // note: AVX has no 256-bit integer sign extension, so each half of the
    // register is extended separately with SSE4.1
    __m128i lo, hi;
    __m256 cv;
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8), maskload can't
    // load 16-bit integers, instead load the last 8 elements (overlapping the
    // previous register) and only store the last rem
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(rmasks[rem])); // last rem floats
      lo = _mm_cvtepi16_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[len-8])));
      hi = _mm_cvtepi16_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[len-4])));
      cv = _mm256_cvtepi32_ps(_mm256_insertf128_si256(
          _mm256_castsi128_si256(lo), hi, 1));
      _mm256_maskstore_ps(&out[len-8], msk, cv);
    }
    for(i-=8; i>=0; i-=8) // process 8 real elements per register
    {
      lo = _mm_cvtepi16_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i])));   // [0,1,2,3]
      hi = _mm_cvtepi16_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i+4]))); // [4,5,6,7]
      cv = _mm256_cvtepi32_ps(_mm256_insertf128_si256(
          _mm256_castsi128_si256(lo), hi, 1));
      _mm256_storeu_ps(&out[i], cv);
    }
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out)
  {
    // note: going from the end to the start so this also works in place
    // note: the input is read with memcpy, in place the int16_t input and
    // float output alias and the compiler could otherwise assume they don't
    // and reorder the loop, the memcpy is actually free as it compiles out
    int16_t tmp;
    for(int i=len-1; i>=0; --i)
    {
      std::memcpy(&tmp, &in[i], sizeof(tmp));
      out[i] = static_cast<float>(tmp);
    }
    return;
  }
  inline void cnvrtInt16ToFloat(void * const inout, const int len)
  {
    cnvrtInt16ToFloat(reinterpret_cast<int16_t const *>(inout), len,
        reinterpret_cast<float *>(inout));
  }
  
  // 16-bit integer to float conversion with min/max output
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f,avx512bw,avx512vl")))
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    // note: the output is 2x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    // note: the min/max is found on the 16-bit integers since converting
    // them to float doesn't change their order, that way each min and max
    // covers 32 elements (2 registers of output)
    __m512i ld;
    __m512i mn = _mm512_set1_epi16(in[0]);
    __m512i mx = mn;
    __m512 cv;
    int i = len&~31; // start of the remaining elements
    // handle remaining elements first (note len&31 == len%32), the min/max
    // is found on all of them at once, then they are converted one register
    // at a time
    const int rem = len&31;
    if(rem)
    {
      const __mmask32 mk = 0xffffffffu >> (32-rem); // first rem
      ld = _mm512_maskz_loadu_epi16(mk, &in[i]);
      // only compare the loaded elements, the others keep their min/max
      mn = _mm512_mask_min_epi16(mn, mk, ld, mn);
      mx = _mm512_mask_max_epi16(mx, mk, ld, mx);
      // the last (partial) register first so that this works in place
      int j = i + (rem&~15); // start of the final partial register
      if(rem&15)
      {
        const __mmask16 mk16 = mask16(rem&15);
        cv = _mm512_cvtepi32_ps(_mm512_cvtepi16_epi32(
            _mm256_maskz_loadu_epi16(mk16, &in[j])));
        _mm512_mask_storeu_ps(&out[j], mk16, cv);
      }
      for(j-=16; j>=i; j-=16) // full registers
      {
        cv = _mm512_cvtepi32_ps(_mm512_cvtepi16_epi32(_mm256_loadu_si256(
            reinterpret_cast<__m256i const *>(&in[j]))));
        _mm512_storeu_ps(&out[j], cv);
      }
    }
    for(i-=32; i>=0; i-=32) // process 32 real elements (2 registers)
    {
      ld = _mm512_loadu_si512(&in[i]);
      mn = _mm512_min_epi16(ld, mn);
      mx = _mm512_max_epi16(ld, mx);
      // convert the last register first so that this works in place (the
      // integers are loaded again, which is faster than extracting them)
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi16_epi32(_mm256_loadu_si256(
          reinterpret_cast<__m256i const *>(&in[i+16]))));
      _mm512_storeu_ps(&out[i+16], cv);
      cv = _mm512_cvtepi32_ps(_mm512_cvtepi16_epi32(_mm256_loadu_si256(
          reinterpret_cast<__m256i const *>(&in[i]))));
      _mm512_storeu_ps(&out[i], cv);
    }
    // sign extend the word min/max 16 at a time to reduce them to 16
    // integers, then the horizontal min and max
    __m512i r = _mm512_min_epi32(_mm512_cvtepi16_epi32(
        _mm512_castsi512_si256(mn)),
        _mm512_cvtepi16_epi32(_mm512_extracti64x4_epi64(mn, 1)));
    const float fmn = static_cast<float>(_mm512_reduce_min_epi32(r));
    r = _mm512_max_epi32(_mm512_cvtepi16_epi32(
        _mm512_castsi512_si256(mx)),
        _mm512_cvtepi16_epi32(_mm512_extracti64x4_epi64(mx, 1)));
    const float fmx = static_cast<float>(_mm512_reduce_max_epi32(r));
    minOut = noInit? std::min(minOut, fmn) : fmn;
    maxOut = noInit? std::max(maxOut, fmx) : fmx;
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX2) // AVX2
  __attribute__((__target__("avx2")))
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    // note: the output is 2x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    if(len < 8) // less than 1 register, convert one element at a time
    {
      if(!noInit)
      {
        minOut = static_cast<float>(in[0]);
        maxOut = static_cast<float>(in[0]);
      }
      // note: the input is read with memcpy, in place the int16_t input and
      // float output alias and the compiler could otherwise assume they don't
      // and reorder the loop, the memcpy is actually free as it compiles out
      int16_t tmp;
      for(int i=len-1; i>=0; --i)
      {
        std::memcpy(&tmp, &in[i], sizeof(tmp));
        out[i] = static_cast<float>(tmp);
        minOut = std::min(minOut, out[i]);
        maxOut = std::max(maxOut, out[i]);
      }
      return;
    }
    // note: the min/max is found on the 16-bit integers since converting
    // them to float doesn't change their order, that way each min and max
    // covers 16 elements (2 registers of output)
    __m256i ld;
    __m256i mn = _mm256_set1_epi16(in[0]);
    __m256i mx = mn;
    __m256 cv, cva, cvb;
    __m128i hs;
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8), maskload can't
    // load 16-bit integers, instead load the last 8 elements (overlapping the
    // previous register) and only store the last rem
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(rmasks[rem])); // last rem floats
      // for the min/max the 8 loaded integers are repeated to fill the register
      ld = _mm256_broadcastsi128_si256(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[len-8])));
      // the overlapping elements are compared twice, which doesn't matter
      mn = _mm256_min_epi16(ld, mn);
      mx = _mm256_max_epi16(ld, mx);
      cv = _mm256_cvtepi32_ps(_mm256_cvtepi16_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[len-8]))));
      _mm256_maskstore_ps(&out[len-8], msk, cv);
    }
    // then a single register if needed to make the rest a multiple of 16
    for(i-=8; i>=(len&~15); i-=8)
    {
      ld = _mm256_broadcastsi128_si256(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i])));
      mn = _mm256_min_epi16(ld, mn);
      mx = _mm256_max_epi16(ld, mx);
      cv = _mm256_cvtepi32_ps(_mm256_cvtepi16_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i]))));
      _mm256_storeu_ps(&out[i], cv);
    }
    for(i=(len&~15)-16; i>=0; i-=16) // process 16 real elements (2 registers)
    {
      ld = _mm256_loadu_si256(reinterpret_cast<__m256i const *>(&in[i]));
      mn = _mm256_min_epi16(ld, mn);
      mx = _mm256_max_epi16(ld, mx);
      // load and convert every register before storing any of them so that
      // this works in place (the integers are loaded again, which is faster
      // than extracting them)
      cvb = _mm256_cvtepi32_ps(_mm256_cvtepi16_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i+8]))));
      cva = _mm256_cvtepi32_ps(_mm256_cvtepi16_epi32(_mm_loadu_si128(
          reinterpret_cast<__m128i const *>(&in[i]))));
      _mm256_storeu_ps(&out[i+8], cvb);
      _mm256_storeu_ps(&out[i], cva);
    }
    // sign extend the word min/max 8 at a time to reduce them to 8
    // integers, then the horizontal min and max
    mn = _mm256_min_epi32(_mm256_cvtepi16_epi32(_mm256_castsi256_si128(mn)),
        _mm256_cvtepi16_epi32(_mm256_extracti128_si256(mn, 1)));
    mx = _mm256_max_epi32(_mm256_cvtepi16_epi32(_mm256_castsi256_si128(mx)),
        _mm256_cvtepi16_epi32(_mm256_extracti128_si256(mx, 1)));
    hs = _mm_min_epi32(_mm256_castsi256_si128(mn),
        _mm256_extracti128_si256(mn, 1)); // [min(0,4),min(1,5),...]
    hs = _mm_min_epi32(hs, _mm_shuffle_epi32(hs, 0x4e)); // [min(0,4,2,6),...]
    hs = _mm_min_epi32(hs, _mm_shuffle_epi32(hs, 0xb1)); // [min(0,...,7),...]
    const float fmn = static_cast<float>(_mm_cvtsi128_si32(hs));
    hs = _mm_max_epi32(_mm256_castsi256_si128(mx),
        _mm256_extracti128_si256(mx, 1)); // [max(0,4),max(1,5),...]
    hs = _mm_max_epi32(hs, _mm_shuffle_epi32(hs, 0x4e)); // [max(0,4,2,6),...]
    hs = _mm_max_epi32(hs, _mm_shuffle_epi32(hs, 0xb1)); // [max(0,...,7),...]
    const float fmx = static_cast<float>(_mm_cvtsi128_si32(hs));
    minOut = noInit? std::min(minOut, fmn) : fmn;
    maxOut = noInit? std::max(maxOut, fmx) : fmx;
    return;
  }
  #endif // end AVX2
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    // note: the output is 2x the size of the input, so to work in place the
    // buffers are processed from the end to the start, that way each store
    // only overwrites input that has already been loaded
    if(len < 8) // less than 1 register, convert one element at a time
    {
      if(!noInit)
      {
        minOut = static_cast<float>(in[0]);
        maxOut = static_cast<float>(in[0]);
      }
      // note: the input is read with memcpy, in place the int16_t input and
      // float output alias and the compiler could otherwise assume they don't
      // and reorder the loop, the memcpy is actually free as it compiles out
      int16_t tmp;
      for(int i=len-1; i>=0; --i)
      {
        std::memcpy(&tmp, &in[i], sizeof(tmp));
        out[i] = static_cast<float>(tmp);
        minOut = std::min(minOut, out[i]);
        maxOut = std::max(maxOut, out[i]);
      }
      return;
    }
    // note: AVX has no 256-bit integer sign extension or min/max, so each
    // half of the register is extended and compared separately with SSE4.1
    // note: the min/max is found on the integers since converting them to
    // float doesn't change their order, integer min/max has a 1 cycle latency
    // (vs. 4 for float) so one min and max accumulator doesn't stall the loop
    __m128i lo, hi;
    __m128i mn = _mm_set1_epi32(in[0]);
    __m128i mx = mn;
    int i = len&~7; // start of the final partial register
    // handle remaining elements first (note len&7 == len%8), maskload can't
    // load 16-bit integers, instead load the last 8 elements (overlapping the
    // previous register) and only store the last rem
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(rmasks[rem])); // last rem floats
      lo = _mm_cvtepi16_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[len-8])));
      hi = _mm_cvtepi16_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[len-4])));
      // the overlapping elements are compared twice, which doesn't matter
      mn = _mm_min_epi32(mn, _mm_min_epi32(lo, hi));
      mx = _mm_max_epi32(mx, _mm_max_epi32(lo, hi));
      _mm256_maskstore_ps(&out[len-8], msk, _mm256_cvtepi32_ps(
          _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1)));
    }
    for(i-=8; i>=0; i-=8) // process 8 real elements per register
    {
      lo = _mm_cvtepi16_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i])));   // [0,1,2,3]
      hi = _mm_cvtepi16_epi32(_mm_loadl_epi64(
          reinterpret_cast<__m128i const *>(&in[i+4]))); // [4,5,6,7]
      mn = _mm_min_epi32(mn, _mm_min_epi32(lo, hi));
      mx = _mm_max_epi32(mx, _mm_max_epi32(lo, hi));
      _mm256_storeu_ps(&out[i], _mm256_cvtepi32_ps(
          _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1)));
    }
    // horizontal min and max
    mn = _mm_min_epi32(mn, _mm_shuffle_epi32(mn, 0x4e)); // [min(0,2),...]
    mn = _mm_min_epi32(mn, _mm_shuffle_epi32(mn, 0xb1)); // [min(0,...,3),...]
    mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0x4e)); // [max(0,2),...]
    mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0xb1)); // [max(0,...,3),...]
    const float fmn = static_cast<float>(_mm_cvtsi128_si32(mn));
    const float fmx = static_cast<float>(_mm_cvtsi128_si32(mx));
    minOut = noInit? std::min(minOut, fmn) : fmn;
    maxOut = noInit? std::max(maxOut, fmx) : fmx;
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void cnvrtInt16ToFloat(int16_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    if(!noInit)
    {
      minOut = static_cast<float>(in[0]);
      maxOut = static_cast<float>(in[0]);
    }
    // note: going from the end to the start so this also works in place
    // note: the input is read with memcpy, in place the int16_t input and
    // float output alias and the compiler could otherwise assume they don't
    // and reorder the loop, the memcpy is actually free as it compiles out
    int16_t tmp;
    for(int i=len-1; i>=0; --i)
    {
      std::memcpy(&tmp, &in[i], sizeof(tmp));
      out[i] = static_cast<float>(tmp);
      minOut = std::min(minOut, out[i]);
      maxOut = std::max(maxOut, out[i]);
    }
    return;
  }

  // 32-bit int to float conversion
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f")))
  inline void cnvrtInt32ToFloat(int32_t const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    __m512 cv;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      cv = _mm512_cvtepi32_ps(_mm512_loadu_si512(&in[i]));
      _mm512_storeu_ps(&out[i], cv);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      cv = _mm512_cvtepi32_ps(_mm512_maskz_loadu_epi32(mk, &in[i]));
      _mm512_mask_storeu_ps(&out[i], mk, cv);
    }
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void cnvrtInt32ToFloat(int32_t const * const in, const int len,
      float * const out)
  {
    if(len<=0) return;
    __m256 cv;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      cv = _mm256_cvtepi32_ps(_mm256_loadu_si256(
          reinterpret_cast<__m256i const *>(&in[i])));
      _mm256_storeu_ps(&out[i], cv);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      // integer maskload needs AVX2, so load the integers as floats
      cv = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      cv = _mm256_cvtepi32_ps(_mm256_castps_si256(cv));
      _mm256_maskstore_ps(&out[i], msk, cv);
    }
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void cnvrtInt32ToFloat(int32_t const * const in, const int len,
      float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = static_cast<float>(in[i]);
    return;
  }
  inline void cnvrtInt32ToFloat(void * const inout, const int len)
  {
    cnvrtInt32ToFloat(reinterpret_cast<int32_t const *>(inout), len,
        reinterpret_cast<float *>(inout));
  }

  // 32-bit integer to float conversion with min/max output
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f")))
  inline void cnvrtInt32ToFloat(int32_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    // note: the min/max is found on the integers since converting them to
    // float doesn't change their order, integer min/max has a 1 cycle latency
    // (vs. 4 for float) so one min and max accumulator doesn't stall the loop
    __m512i ld;
    __m512i mn = _mm512_set1_epi32(in[0]);
    __m512i mx = mn;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld = _mm512_loadu_si512(&in[i]);
      mn = _mm512_min_epi32(ld, mn);
      mx = _mm512_max_epi32(ld, mx);
      _mm512_storeu_ps(&out[i], _mm512_cvtepi32_ps(ld));
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16(rem);
      ld = _mm512_maskz_loadu_epi32(mk, &in[i]);
      // only compare the loaded elements, the others keep their min/max
      mn = _mm512_mask_min_epi32(mn, mk, ld, mn);
      mx = _mm512_mask_max_epi32(mx, mk, ld, mx);
      _mm512_mask_storeu_ps(&out[i], mk, _mm512_cvtepi32_ps(ld));
    }
    const float fmn = static_cast<float>(_mm512_reduce_min_epi32(mn));
    const float fmx = static_cast<float>(_mm512_reduce_max_epi32(mx));
    minOut = noInit? std::min(minOut, fmn) : fmn;
    maxOut = noInit? std::max(maxOut, fmx) : fmx;
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX2) // AVX2
  __attribute__((__target__("avx2")))
  inline void cnvrtInt32ToFloat(int32_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    // note: the min/max is found on the integers since converting them to
    // float doesn't change their order, integer min/max has a 1 cycle latency
    // (vs. 4 for float) so one min and max accumulator doesn't stall the loop
    __m256i ld;
    __m256i mn = _mm256_set1_epi32(in[0]);
    __m256i mx = mn;
    __m128i hs;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld = _mm256_loadu_si256(reinterpret_cast<__m256i const *>(&in[i]));
      mn = _mm256_min_epi32(ld, mn);
      mx = _mm256_max_epi32(ld, mx);
      _mm256_storeu_ps(&out[i], _mm256_cvtepi32_ps(ld));
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      ld = _mm256_maskload_epi32(&in[i], msk);
      // maskload sets the unused elements to 0, which could become the min or
      // max, so set them to in[i] instead (always loaded)
      ld = _mm256_or_si256(ld, _mm256_andnot_si256(msk,
          _mm256_set1_epi32(in[i])));
      mn = _mm256_min_epi32(ld, mn);
      mx = _mm256_max_epi32(ld, mx);
      _mm256_maskstore_ps(&out[i], msk, _mm256_cvtepi32_ps(ld));
    }
    // horizontal min and max
    hs = _mm_min_epi32(_mm256_castsi256_si128(mn),
        _mm256_extracti128_si256(mn, 1)); // [min(0,4),min(1,5),...]
    hs = _mm_min_epi32(hs, _mm_shuffle_epi32(hs, 0x4e)); // [min(0,4,2,6),...]
    hs = _mm_min_epi32(hs, _mm_shuffle_epi32(hs, 0xb1)); // [min(0,...,7),...]
    const float fmn = static_cast<float>(_mm_cvtsi128_si32(hs));
    hs = _mm_max_epi32(_mm256_castsi256_si128(mx),
        _mm256_extracti128_si256(mx, 1)); // [max(0,4),max(1,5),...]
    hs = _mm_max_epi32(hs, _mm_shuffle_epi32(hs, 0x4e)); // [max(0,4,2,6),...]
    hs = _mm_max_epi32(hs, _mm_shuffle_epi32(hs, 0xb1)); // [max(0,...,7),...]
    const float fmx = static_cast<float>(_mm_cvtsi128_si32(hs));
    minOut = noInit? std::min(minOut, fmn) : fmn;
    maxOut = noInit? std::max(maxOut, fmx) : fmx;
    return;
  }
  #endif // end AVX2
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void cnvrtInt32ToFloat(int32_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    // note: AVX has no 256-bit integer min/max, so the two halves of each
    // register are compared separately with SSE4.1
    // note: the min/max is found on the integers since converting them to
    // float doesn't change their order, integer min/max has a 1 cycle latency
    // (vs. 4 for float) so one min and max accumulator doesn't stall the loop
    __m256i ld;
    __m128i lo, hi;
    __m128i mn = _mm_set1_epi32(in[0]);
    __m128i mx = mn;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld = _mm256_loadu_si256(reinterpret_cast<__m256i const *>(&in[i]));
      lo = _mm256_castsi256_si128(ld);      // [0,1,2,3]
      hi = _mm256_extractf128_si256(ld, 1); // [4,5,6,7]
      mn = _mm_min_epi32(mn, _mm_min_epi32(lo, hi));
      mx = _mm_max_epi32(mx, _mm_max_epi32(lo, hi));
      _mm256_storeu_ps(&out[i], _mm256_cvtepi32_ps(ld));
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[rem]));
      // integer maskload needs AVX2, so load the integers as floats
      __m256 fl = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]),
          msk);
      // maskload sets the unused elements to 0, which could become the min or
      // max, so set them to in[i] instead (always loaded)
      fl = _mm256_or_ps(fl, _mm256_andnot_ps(_mm256_castsi256_ps(msk),
          _mm256_broadcast_ss(reinterpret_cast<float const *>(&in[i]))));
      ld = _mm256_castps_si256(fl);
      lo = _mm256_castsi256_si128(ld);      // [0,1,2,3]
      hi = _mm256_extractf128_si256(ld, 1); // [4,5,6,7]
      mn = _mm_min_epi32(mn, _mm_min_epi32(lo, hi));
      mx = _mm_max_epi32(mx, _mm_max_epi32(lo, hi));
      _mm256_maskstore_ps(&out[i], msk, _mm256_cvtepi32_ps(ld));
    }
    // horizontal min and max
    mn = _mm_min_epi32(mn, _mm_shuffle_epi32(mn, 0x4e)); // [min(0,2),...]
    mn = _mm_min_epi32(mn, _mm_shuffle_epi32(mn, 0xb1)); // [min(0,...,3),...]
    mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0x4e)); // [max(0,2),...]
    mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0xb1)); // [max(0,...,3),...]
    const float fmn = static_cast<float>(_mm_cvtsi128_si32(mn));
    const float fmx = static_cast<float>(_mm_cvtsi128_si32(mx));
    minOut = noInit? std::min(minOut, fmn) : fmn;
    maxOut = noInit? std::max(maxOut, fmx) : fmx;
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void cnvrtInt32ToFloat(int32_t const * const in, const int len,
      float * const out, float &minOut, float &maxOut, const bool noInit)
  {
    if(len<=0) return;
    if(!noInit)
    {
      minOut = static_cast<float>(in[0]);
      maxOut = static_cast<float>(in[0]);
    }
    for(int i=0; i<len; ++i)
    {
      out[i] = static_cast<float>(in[i]);
      minOut = std::min(minOut, out[i]);
      maxOut = std::max(maxOut, out[i]);
    }
    return;
  }

  // Polar to cartesian conversion with phase in cycles
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f")))
  inline void polToCart(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, mag, ph, sn, cs;
    int i = 0;
    // note: the magnitudes and phases of 2 registers are separated so that
    // each sincos call computes 16 phases (m = magnitude, p = phase)
    for(; i<len-15; i+=16) // process 16 complex elements (2 registers)
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
      mag = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm512_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      TRIG::sincosCyc512(ph, sn, cs);          // sin(2*pi*ph), cos(2*pi*ph)
      cs = _mm512_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm512_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm512_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm512_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld2);
    }
    // handle remaining elements (note len&15 == len%16), the sincos call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 8)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16((rem>8? 8 : rem)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(&in[i]));
      ld2 = _mm512_setzero_ps();
      if(rem>8)
        ld2 = _mm512_maskz_loadu_ps(mask16((rem-8)<<1),
            reinterpret_cast<float const *>(&in[i+8]));
      mag = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm512_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      TRIG::sincosCyc512(ph, sn, cs);          // sin(2*pi*ph), cos(2*pi*ph)
      cs = _mm512_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm512_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm512_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm512_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
      if(rem>8)
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]),
            mask16((rem-8)<<1), ld2);
    }
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX2) // AVX2
  __attribute__((__target__("avx2,fma")))
  inline void polToCart(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    // note: the same as the AVX version, but evm_trig.h's sin/cos uses FMA
    __m256 ld1, ld2, mag, ph, sn, cs;
    int i = 0;
    // note: the magnitudes and phases of 2 registers are separated so that
    // each sincos call computes 8 phases (m = magnitude, p = phase)
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      mag = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm256_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      TRIG::sincosCycAvx2(ph, sn, cs);         // sin(2*pi*ph), cos(2*pi*ph)
      cs = _mm256_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm256_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm256_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm256_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8), the sincos call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 4)
    const int rem = len&7;
    if(rem)
    {
      // note msk accounts for 2 reals per element for the complex buffers
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem>4? 4 : rem)<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      ld2 = _mm256_setzero_ps();
      if(rem>4)
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
      mag = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm256_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      TRIG::sincosCycAvx2(ph, sn, cs);         // sin(2*pi*ph), cos(2*pi*ph)
      cs = _mm256_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm256_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm256_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm256_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
      if(rem>4)
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])), ld2);
    }
    return;
  }
  #endif // end AVX2
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void polToCart(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, mag, ph, sn, cs;
    int i = 0;
    // note: the magnitudes and phases of 2 registers are separated so that
    // each sincos call computes 8 phases (m = magnitude, p = phase)
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      mag = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm256_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      TRIG::sincosCycAvx(ph, sn, cs);          // sin(2*pi*ph), cos(2*pi*ph)
      cs = _mm256_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm256_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm256_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm256_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8), the sincos call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 4)
    const int rem = len&7;
    if(rem)
    {
      // note msk accounts for 2 reals per element for the complex buffers
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem>4? 4 : rem)<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      ld2 = _mm256_setzero_ps();
      if(rem>4)
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
      mag = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm256_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      TRIG::sincosCycAvx(ph, sn, cs);          // sin(2*pi*ph), cos(2*pi*ph)
      cs = _mm256_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm256_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm256_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm256_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
      if(rem>4)
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])), ld2);
    }
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void polToCart(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    float tmpx, tmpy, mag;
    // note: whole cycles are removed from the phase first, so large phases
    // keep their accuracy
    for(int i=0; i<len; ++i)
    {
      mag = in[i].real();
      sincosf(M_TWO_PI*(in[i].imag() - rintf(in[i].imag())), &tmpy, &tmpx);
      out[i].real(mag*tmpx);
      out[i].imag(mag*tmpy);
    }
    return;
  }

  // Polar to cartesian conversion with phase in radians
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f")))
  inline void polToCartRad(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, mag, ph, sn, cs;
    int i = 0;
    // note: the magnitudes and phases of 2 registers are separated so that
    // each sincos call computes 16 phases (m = magnitude, p = phase)
    for(; i<len-15; i+=16) // process 16 complex elements (2 registers)
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
      mag = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm512_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      #if defined(EVM_SVML)
      sn = _mm512_sincos_ps(&cs, ph);          // sin(ph), cos(ph) in cs
      #else
      TRIG::sincosRad512(ph, sn, cs);          // sin(ph), cos(ph) in cs
      #endif
      cs = _mm512_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm512_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm512_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm512_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld2);
    }
    // handle remaining elements (note len&15 == len%16), the sincos call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 8)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16((rem>8? 8 : rem)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(&in[i]));
      ld2 = _mm512_setzero_ps();
      if(rem>8)
        ld2 = _mm512_maskz_loadu_ps(mask16((rem-8)<<1),
            reinterpret_cast<float const *>(&in[i+8]));
      mag = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm512_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      #if defined(EVM_SVML)
      sn = _mm512_sincos_ps(&cs, ph);          // sin(ph), cos(ph) in cs
      #else
      TRIG::sincosRad512(ph, sn, cs);          // sin(ph), cos(ph) in cs
      #endif
      cs = _mm512_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm512_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm512_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm512_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
      if(rem>8)
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]),
            mask16((rem-8)<<1), ld2);
    }
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX2) // AVX2
  __attribute__((__target__("avx2,fma")))
  inline void polToCartRad(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    // note: the same as the AVX version, but SVML's AVX2 sincos is about 2.4x
    // faster than its AVX one, and evm_trig.h's uses FMA
    __m256 ld1, ld2, mag, ph, sn, cs;
    int i = 0;
    // note: the magnitudes and phases of 2 registers are separated so that
    // each sincos call computes 8 phases (m = magnitude, p = phase)
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      mag = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm256_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      #if defined(EVM_SVML)
      sn = _mm256_sincos_ps(&cs, ph);          // sin(ph), cos(ph) in cs
      #else
      TRIG::sincosRadAvx2(ph, sn, cs);         // sin(ph), cos(ph) in cs
      #endif
      cs = _mm256_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm256_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm256_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm256_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8), the sincos call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 4)
    const int rem = len&7;
    if(rem)
    {
      // note msk accounts for 2 reals per element for the complex buffers
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem>4? 4 : rem)<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      ld2 = _mm256_setzero_ps();
      if(rem>4)
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
      mag = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm256_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      #if defined(EVM_SVML)
      sn = _mm256_sincos_ps(&cs, ph);          // sin(ph), cos(ph) in cs
      #else
      TRIG::sincosRadAvx2(ph, sn, cs);         // sin(ph), cos(ph) in cs
      #endif
      cs = _mm256_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm256_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm256_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm256_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
      if(rem>4)
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])), ld2);
    }
    return;
  }
  #endif // end AVX2
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void polToCartRad(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, mag, ph, sn, cs;
    int i = 0;
    // note: the magnitudes and phases of 2 registers are separated so that
    // each sincos call computes 8 phases (m = magnitude, p = phase)
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      mag = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm256_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      #if defined(EVM_SVML)
      sn = _mm256_sincos_ps(&cs, ph);          // sin(ph), cos(ph) in cs
      #else
      TRIG::sincosRadAvx(ph, sn, cs);          // sin(ph), cos(ph) in cs
      #endif
      cs = _mm256_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm256_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm256_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm256_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8), the sincos call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 4)
    const int rem = len&7;
    if(rem)
    {
      // note msk accounts for 2 reals per element for the complex buffers
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem>4? 4 : rem)<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      ld2 = _mm256_setzero_ps();
      if(rem>4)
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
      mag = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Am0,Am1,Bm0,Bm1,...]
      ph = _mm256_shuffle_ps(ld1, ld2, 0xdd);  // [Ap0,Ap1,Bp0,Bp1,...]
      #if defined(EVM_SVML)
      sn = _mm256_sincos_ps(&cs, ph);          // sin(ph), cos(ph) in cs
      #else
      TRIG::sincosRadAvx(ph, sn, cs);          // sin(ph), cos(ph) in cs
      #endif
      cs = _mm256_mul_ps(mag, cs);             // [Ax0,Ax1,Bx0,Bx1,...]
      sn = _mm256_mul_ps(mag, sn);             // [Ay0,Ay1,By0,By1,...]
      ld1 = _mm256_unpacklo_ps(cs, sn);        // [Ax0,Ay0,Ax1,Ay1,...]
      ld2 = _mm256_unpackhi_ps(cs, sn);        // [Bx0,By0,Bx1,By1,...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
      if(rem>4)
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])), ld2);
    }
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void polToCartRad(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    float tmpx, tmpy, mag;
    for(int i=0; i<len; ++i)
    {
      mag = in[i].real();
      sincosf(in[i].imag(), &tmpy, &tmpx);
      out[i].real(mag*tmpx);
      out[i].imag(mag*tmpy);
    }
    return;
  }

  // Cartesian to polar conversion with phase in cycles
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f")))
  inline void cartToPol(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, re, im, mag, ang;
    #if defined(EVM_SVML)
    // SVML's atan2 takes a much slower path for any register with a (0,0)
    // input (up to 5x slower on zero padded data), and gives NaN when both
    // are denormal and denormals are treated as zero (DAZ, icpx's default).
    // atan2(+-0,+-1) gives the same result as atan2(+-0,+-0), so re is set
    // to +-1 where re|im compares equal to 0 (both zero, or both denormal
    // under DAZ) before the atan2 (or'ing in the bits of 1.0 keeps the
    // sign), the magnitude still uses the original re
    const __m512i one = _mm512_castps_si512(_mm512_set1_ps(1.0f));
    __mmask16 zr;
    const __m512 sc = _mm512_set1_ps(static_cast<float>(RAD_TO_CYC));
    #endif
    int i = 0;
    // note: the real and imaginary parts of 2 registers are separated so
    // that each atan2 call computes 16 angles
    for(; i<len-15; i+=16) // process 16 complex elements (2 registers)
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm512_mul_ps(im, im);
      mag = _mm512_fmadd_ps(re, re, mag); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
      mag = _mm512_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      zr = _mm512_cmp_ps_mask(_mm512_castsi512_ps(_mm512_or_epi32(
          _mm512_castps_si512(re), _mm512_castps_si512(im))),
          _mm512_setzero_ps(), _CMP_EQ_OQ); // both 0 (or denormal under DAZ)
      ang = _mm512_castsi512_ps(_mm512_mask_or_epi32(_mm512_castps_si512(re),
          zr, _mm512_castps_si512(re), one)); // +-1 where both are +-0
      ang = _mm512_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      ang = _mm512_mul_ps(ang, sc);       // radians to cycles
      #else
      ang = TRIG::atan2Cyc512(im, re);    // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm512_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm512_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld2);
    }
    // handle remaining elements (note len&15 == len%16), the atan2 call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 8)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16((rem>8? 8 : rem)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(&in[i]));
      ld2 = _mm512_setzero_ps();
      if(rem>8)
        ld2 = _mm512_maskz_loadu_ps(mask16((rem-8)<<1),
            reinterpret_cast<float const *>(&in[i+8]));
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm512_mul_ps(im, im);
      mag = _mm512_fmadd_ps(re, re, mag); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
      mag = _mm512_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      zr = _mm512_cmp_ps_mask(_mm512_castsi512_ps(_mm512_or_epi32(
          _mm512_castps_si512(re), _mm512_castps_si512(im))),
          _mm512_setzero_ps(), _CMP_EQ_OQ); // both 0 (or denormal under DAZ)
      ang = _mm512_castsi512_ps(_mm512_mask_or_epi32(_mm512_castps_si512(re),
          zr, _mm512_castps_si512(re), one)); // +-1 where both are +-0
      ang = _mm512_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      ang = _mm512_mul_ps(ang, sc);       // radians to cycles
      #else
      ang = TRIG::atan2Cyc512(im, re);    // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm512_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm512_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
      if(rem>8)
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]),
            mask16((rem-8)<<1), ld2);
    }
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX2) // AVX2
  __attribute__((__target__("avx2,fma")))
  inline void cartToPol(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    // note: the same as the AVX version except for the FMA in the magnitude,
    // SVML's AVX2 atan2 is 1.2-1.4x faster than its AVX one
    __m256 ld1, ld2, re, im, mag, ang;
    #if defined(EVM_SVML)
    // SVML's atan2 gives NaN when re and im are both denormal and denormals
    // are treated as zero (DAZ, icpx's default), atan2(+-0,+-1) gives the same
    // result as atan2(+-0,+-0), so re is set to +-1 where re|im compares
    // equal to 0 (or'ing in the bits of 1.0 keeps the sign)
    const __m256 one = _mm256_set1_ps(1.0f);
    const __m256 sc = _mm256_set1_ps(static_cast<float>(RAD_TO_CYC));
    #endif
    int i = 0;
    // note: the real and imaginary parts of 2 registers are separated so
    // that each atan2 call computes 8 angles
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm256_mul_ps(im, im);
      mag = _mm256_fmadd_ps(re, re, mag); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
      mag = _mm256_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      ang = _mm256_cmp_ps(_mm256_or_ps(re, im), _mm256_setzero_ps(),
          _CMP_EQ_OQ); // re and im both 0 (or denormal under DAZ)
      ang = _mm256_or_ps(re, _mm256_and_ps(ang, one)); // +-1 where both 0
      ang = _mm256_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      ang = _mm256_mul_ps(ang, sc);       // radians to cycles
      #else
      ang = TRIG::atan2CycAvx2(im, re);   // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm256_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm256_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8), the atan2 call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 4)
    const int rem = len&7;
    if(rem)
    {
      // note msk accounts for 2 reals per element for the complex buffers
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem>4? 4 : rem)<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      ld2 = _mm256_setzero_ps();
      if(rem>4)
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm256_mul_ps(im, im);
      mag = _mm256_fmadd_ps(re, re, mag); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
      mag = _mm256_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      ang = _mm256_cmp_ps(_mm256_or_ps(re, im), _mm256_setzero_ps(),
          _CMP_EQ_OQ); // re and im both 0 (or denormal under DAZ)
      ang = _mm256_or_ps(re, _mm256_and_ps(ang, one)); // +-1 where both 0
      ang = _mm256_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      ang = _mm256_mul_ps(ang, sc);       // radians to cycles
      #else
      ang = TRIG::atan2CycAvx2(im, re);   // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm256_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm256_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
      if(rem>4)
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])), ld2);
    }
    return;
  }
  #endif // end AVX2
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void cartToPol(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, re, im, mag, ang;
    #if defined(EVM_SVML)
    // SVML's atan2 gives NaN when re and im are both denormal and denormals
    // are treated as zero (DAZ, icpx's default), atan2(+-0,+-1) gives the same
    // result as atan2(+-0,+-0), so re is set to +-1 where re|im compares
    // equal to 0 (or'ing in the bits of 1.0 keeps the sign)
    const __m256 one = _mm256_set1_ps(1.0f);
    const __m256 sc = _mm256_set1_ps(static_cast<float>(RAD_TO_CYC));
    #endif
    int i = 0;
    // note: the real and imaginary parts of 2 registers are separated so
    // that each atan2 call computes 8 angles
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm256_add_ps(_mm256_mul_ps(re, re), _mm256_mul_ps(im, im));
      mag = _mm256_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      ang = _mm256_cmp_ps(_mm256_or_ps(re, im), _mm256_setzero_ps(),
          _CMP_EQ_OQ); // re and im both 0 (or denormal under DAZ)
      ang = _mm256_or_ps(re, _mm256_and_ps(ang, one)); // +-1 where both 0
      ang = _mm256_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      ang = _mm256_mul_ps(ang, sc);       // radians to cycles
      #else
      ang = TRIG::atan2CycAvx(im, re);    // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm256_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm256_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8), the atan2 call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 4)
    const int rem = len&7;
    if(rem)
    {
      // note msk accounts for 2 reals per element for the complex buffers
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem>4? 4 : rem)<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      ld2 = _mm256_setzero_ps();
      if(rem>4)
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm256_add_ps(_mm256_mul_ps(re, re), _mm256_mul_ps(im, im));
      mag = _mm256_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      ang = _mm256_cmp_ps(_mm256_or_ps(re, im), _mm256_setzero_ps(),
          _CMP_EQ_OQ); // re and im both 0 (or denormal under DAZ)
      ang = _mm256_or_ps(re, _mm256_and_ps(ang, one)); // +-1 where both 0
      ang = _mm256_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      ang = _mm256_mul_ps(ang, sc);       // radians to cycles
      #else
      ang = TRIG::atan2CycAvx(im, re);    // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm256_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm256_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
      if(rem>4)
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])), ld2);
    }
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void cartToPol(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    // note: where im is 0 the angle only depends on the sign of re, so re is
    // moved away from 0 by +-1 there (keeping its sign, without a branch):
    // icpx vectorizes this loop with SVML, whose atan2 gives NaN when re and
    // im are both denormal and denormals are treated as zero (DAZ, icpx's
    // default)
    float ang, mag;
    for(int i=0; i<len; ++i)
    {
      mag = sqrtf(in[i].real()*in[i].real()+in[i].imag()*in[i].imag());
      ang = atan2f(in[i].imag(), in[i].real() +
          copysignf(in[i].imag() == 0.0f, in[i].real()));
      out[i].real(mag);
      out[i].imag(ang*RAD_TO_CYC);
    }
    return;
  }

  // Cartesian to polar conversion with phase in radians
  #if !defined(DISABLE_AVX512) // AVX512 
  __attribute__((__target__("avx512f")))
  inline void cartToPolRad(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m512 ld1, ld2, re, im, mag, ang;
    #if defined(EVM_SVML)
    // SVML's atan2 takes a much slower path for any register with a (0,0)
    // input (up to 5x slower on zero padded data), and gives NaN when both
    // are denormal and denormals are treated as zero (DAZ, icpx's default).
    // atan2(+-0,+-1) gives the same result as atan2(+-0,+-0), so re is set
    // to +-1 where re|im compares equal to 0 (both zero, or both denormal
    // under DAZ) before the atan2 (or'ing in the bits of 1.0 keeps the
    // sign), the magnitude still uses the original re
    const __m512i one = _mm512_castps_si512(_mm512_set1_ps(1.0f));
    __mmask16 zr;
    #endif
    int i = 0;
    // note: the real and imaginary parts of 2 registers are separated so
    // that each atan2 call computes 16 angles
    for(; i<len-15; i+=16) // process 16 complex elements (2 registers)
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const *>(&in[i+8])); // B
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm512_mul_ps(im, im);
      mag = _mm512_fmadd_ps(re, re, mag); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
      mag = _mm512_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      zr = _mm512_cmp_ps_mask(_mm512_castsi512_ps(_mm512_or_epi32(
          _mm512_castps_si512(re), _mm512_castps_si512(im))),
          _mm512_setzero_ps(), _CMP_EQ_OQ); // both 0 (or denormal under DAZ)
      ang = _mm512_castsi512_ps(_mm512_mask_or_epi32(_mm512_castps_si512(re),
          zr, _mm512_castps_si512(re), one)); // +-1 where both are +-0
      ang = _mm512_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      #else
      ang = TRIG::atan2Rad512(im, re);    // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm512_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm512_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm512_storeu_ps(reinterpret_cast<float *>(&out[i+8]), ld2);
    }
    // handle remaining elements (note len&15 == len%16), the atan2 call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 8)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = mask16((rem>8? 8 : rem)<<1); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const *>(&in[i]));
      ld2 = _mm512_setzero_ps();
      if(rem>8)
        ld2 = _mm512_maskz_loadu_ps(mask16((rem-8)<<1),
            reinterpret_cast<float const *>(&in[i+8]));
      re = _mm512_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm512_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm512_mul_ps(im, im);
      mag = _mm512_fmadd_ps(re, re, mag); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
      mag = _mm512_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      zr = _mm512_cmp_ps_mask(_mm512_castsi512_ps(_mm512_or_epi32(
          _mm512_castps_si512(re), _mm512_castps_si512(im))),
          _mm512_setzero_ps(), _CMP_EQ_OQ); // both 0 (or denormal under DAZ)
      ang = _mm512_castsi512_ps(_mm512_mask_or_epi32(_mm512_castps_si512(re),
          zr, _mm512_castps_si512(re), one)); // +-1 where both are +-0
      ang = _mm512_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      #else
      ang = TRIG::atan2Rad512(im, re);    // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm512_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm512_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i]), mk, ld1);
      if(rem>8)
        _mm512_mask_storeu_ps(reinterpret_cast<float *>(&out[i+8]),
            mask16((rem-8)<<1), ld2);
    }
    return;
  }
  #endif // end AVX512
  #if !defined(DISABLE_AVX2) // AVX2
  __attribute__((__target__("avx2,fma")))
  inline void cartToPolRad(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    // note: the same as the AVX version except for the FMA in the magnitude,
    // SVML's AVX2 atan2 is 1.2-1.4x faster than its AVX one
    __m256 ld1, ld2, re, im, mag, ang;
    #if defined(EVM_SVML)
    // SVML's atan2 gives NaN when re and im are both denormal and denormals
    // are treated as zero (DAZ, icpx's default), atan2(+-0,+-1) gives the same
    // result as atan2(+-0,+-0), so re is set to +-1 where re|im compares
    // equal to 0 (or'ing in the bits of 1.0 keeps the sign)
    const __m256 one = _mm256_set1_ps(1.0f);
    #endif
    int i = 0;
    // note: the real and imaginary parts of 2 registers are separated so
    // that each atan2 call computes 8 angles
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm256_mul_ps(im, im);
      mag = _mm256_fmadd_ps(re, re, mag); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
      mag = _mm256_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      ang = _mm256_cmp_ps(_mm256_or_ps(re, im), _mm256_setzero_ps(),
          _CMP_EQ_OQ); // re and im both 0 (or denormal under DAZ)
      ang = _mm256_or_ps(re, _mm256_and_ps(ang, one)); // +-1 where both 0
      ang = _mm256_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      #else
      ang = TRIG::atan2RadAvx2(im, re);   // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm256_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm256_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8), the atan2 call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 4)
    const int rem = len&7;
    if(rem)
    {
      // note msk accounts for 2 reals per element for the complex buffers
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem>4? 4 : rem)<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      ld2 = _mm256_setzero_ps();
      if(rem>4)
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm256_mul_ps(im, im);
      mag = _mm256_fmadd_ps(re, re, mag); // [|A0|^2,|A1|^2,|B0|^2,|B1|^2,...]
      mag = _mm256_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      ang = _mm256_cmp_ps(_mm256_or_ps(re, im), _mm256_setzero_ps(),
          _CMP_EQ_OQ); // re and im both 0 (or denormal under DAZ)
      ang = _mm256_or_ps(re, _mm256_and_ps(ang, one)); // +-1 where both 0
      ang = _mm256_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      #else
      ang = TRIG::atan2RadAvx2(im, re);   // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm256_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm256_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
      if(rem>4)
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])), ld2);
    }
    return;
  }
  #endif // end AVX2
  #if !defined(DISABLE_AVX) // AVX 
  __attribute__((__target__("avx")))
  inline void cartToPolRad(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    if(len<=0) return;
    __m256 ld1, ld2, re, im, mag, ang;
    #if defined(EVM_SVML)
    // SVML's atan2 gives NaN when re and im are both denormal and denormals
    // are treated as zero (DAZ, icpx's default), atan2(+-0,+-1) gives the same
    // result as atan2(+-0,+-0), so re is set to +-1 where re|im compares
    // equal to 0 (or'ing in the bits of 1.0 keeps the sign)
    const __m256 one = _mm256_set1_ps(1.0f);
    #endif
    int i = 0;
    // note: the real and imaginary parts of 2 registers are separated so
    // that each atan2 call computes 8 angles
    for(; i<len-7; i+=8) // process 8 complex elements (2 registers)
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i])); // A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const *>(&in[i+4])); // B
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm256_add_ps(_mm256_mul_ps(re, re), _mm256_mul_ps(im, im));
      mag = _mm256_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      ang = _mm256_cmp_ps(_mm256_or_ps(re, im), _mm256_setzero_ps(),
          _CMP_EQ_OQ); // re and im both 0 (or denormal under DAZ)
      ang = _mm256_or_ps(re, _mm256_and_ps(ang, one)); // +-1 where both 0
      ang = _mm256_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      #else
      ang = TRIG::atan2RadAvx(im, re);    // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm256_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm256_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i]), ld1);
      _mm256_storeu_ps(reinterpret_cast<float *>(&out[i+4]), ld2);
    }
    // handle remaining elements (note len&7 == len%8), the atan2 call takes
    // as long for 1 register as for 2, so the remainder is always processed
    // as 2 registers (the second is zero if rem <= 4)
    const int rem = len&7;
    if(rem)
    {
      // note msk accounts for 2 reals per element for the complex buffers
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const *>(masks[(rem>4? 4 : rem)<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i]), msk);
      ld2 = _mm256_setzero_ps();
      if(rem>4)
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const *>(&in[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])));
      re = _mm256_shuffle_ps(ld1, ld2, 0x88); // [Ar0,Ar1,Br0,Br1,...]
      im = _mm256_shuffle_ps(ld1, ld2, 0xdd); // [Ai0,Ai1,Bi0,Bi1,...]
      mag = _mm256_add_ps(_mm256_mul_ps(re, re), _mm256_mul_ps(im, im));
      mag = _mm256_sqrt_ps(mag);          // [|A0|,|A1|,|B0|,|B1|,...]
      #if defined(EVM_SVML)
      ang = _mm256_cmp_ps(_mm256_or_ps(re, im), _mm256_setzero_ps(),
          _CMP_EQ_OQ); // re and im both 0 (or denormal under DAZ)
      ang = _mm256_or_ps(re, _mm256_and_ps(ang, one)); // +-1 where both 0
      ang = _mm256_atan2_ps(im, ang);     // [arg(A0),arg(A1),arg(B0),...]
      #else
      ang = TRIG::atan2RadAvx(im, re);    // [arg(A0),arg(A1),arg(B0),...]
      #endif
      ld1 = _mm256_unpacklo_ps(mag, ang); // [|A0|,arg(A0),|A1|,arg(A1),...]
      ld2 = _mm256_unpackhi_ps(mag, ang); // [|B0|,arg(B0),|B1|,arg(B1),...]
      _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i]), msk, ld1);
      if(rem>4)
        _mm256_maskstore_ps(reinterpret_cast<float *>(&out[i+4]),
            _mm256_load_si256(reinterpret_cast<__m256i const *>(
            masks[(rem-4)<<1])), ld2);
    }
    return;
  }
  #endif // end AVX 
  __attribute__((__target__("default")))
  inline void cartToPolRad(complex<float> const * const in, const int len,
      complex<float> * const out)
  {
    // note: where im is 0 the angle only depends on the sign of re, so re is
    // moved away from 0 by +-1 there (keeping its sign, without a branch):
    // icpx vectorizes this loop with SVML, whose atan2 gives NaN when re and
    // im are both denormal and denormals are treated as zero (DAZ, icpx's
    // default)
    float ang, mag;
    for(int i=0; i<len; ++i)
    {
      mag = sqrtf(in[i].real()*in[i].real()+in[i].imag()*in[i].imag());
      ang = atan2f(in[i].imag(), in[i].real() +
          copysignf(in[i].imag() == 0.0f, in[i].real()));
      out[i].real(mag);
      out[i].imag(ang);
    }
    return;
  }

  //===============================================//
  // Specialized Functions Function Implementation //
  //===============================================//

} // namespace EVM
#endif // EFFICIENT_VECTOR_MATH_H