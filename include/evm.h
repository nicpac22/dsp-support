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
//  Implements some common efficient vector math functions using SIMD (single
//  instruction multiple data) intrinsics.  These functions use runtime
//  dispatch to automatically choose the most advanced instruction set based
//  on the architecture being run on.  For example, if a function supports both
//  an AVX and an AVX512 branch, it will take the AVX512 branch if executed on
//  an Icelake cpu or the AVX branch if executed on a Haswell cpu, regardless
//  of compiler optimization flags.  In some cases (unit/regression testing,
//  etc.) it may be desireable to disable certain branches.  This can be done
//  by specifying one or more of the following macros at compile time:
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
//  To help mitigate the cost of load operations, this class offers an aligned
//  allocator for stl vectors (aligned_allocator.h) which both aligns memory
//  on boundaries that are a multiple of the SIMD register size and reserves
//  additional space in the buffer equivalent to one full SIMD register.  There
//  is also a helper function to return the number of elements per SIMD
//  register based on data type (float, complex<float>, etc.).  Below is an
//  example of creating a SIMD aligned vector of complex single precision
//  floating point elements, sized to a multiple of the SIMD register size:
//    ALIGNED_VECTOR(complex<float> ) myvec(EVM::getSize<complex<float> >(29));
//
//  The following mathematical operations are currently implemented:
//
//  Multiply/Divide:
//    scale = multiply an input buffer by a constant scale factor
//    mult = point-by-point vector multiply
//    multc = point-by-point complex conjugate multiply, second input conjugated
//    div = point-by-point vector divide
//    divr = point-by-point vector divide using fast reciprocal approximation
//    divc = point-by-point complex conjugate divide, second input conjugated
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
//  Note that for each method there is an "internal" fast version that starts
//  with "_" (i.e. _mult) that avoids all bounds checking (both load and store)
//  but requires the buffer length being operated on to be an integer multiple
//  of the SIMD register size for a given architecture.  The getSize(<type>)
//  method can be used to assist in computing the number of elements to
//  operate on.
//
//  IMPORTANT: unless otherwise specified in the function declaration, when
//  using the internal "_" methods, assume you must set length to
//  getSize<complex<float> >(len) if all arguments are complex or
//  getSize<float>(len) if ANY argument is real.
//
//==============================================================================
#ifndef EFFICIENT_VECTOR_MATH_H
#define EFFICIENT_VECTOR_MATH_H

#include <math.h>
#include <complex.h>
#include <string.h>
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

using std::vector;

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
#define ALIGNED_VECTOR(type) vector<type, aligned_allocator<type, evm_align> >
#endif
#ifndef ALIGNED_ITERATOR
#define ALIGNED_ITERATOR(type) vector<type, aligned_allocator<type, evm_align> >::iterator
#endif

// define some basic macros for doing masked load/store
// macro to generate mask for AVX512 maskz_load and mask_store, if output bit
// is 0, that float does not get loaded/stored, so take the number you want
// loaded/stored, and shift right which will put 16-X zeroes in the msb's and
// X ones in the lsb's
#define MASK16(X) (unsigned short)(0xffff>>(16-X))

// mask table for AVX instructions, rows can be cast as __m256i* directly
alignas(32) static const int32_t masks[8][8] = {
  { 0, 0, 0, 0, 0, 0, 0, 0},   // remainder 0 (unused)
  {-1, 0, 0, 0, 0, 0, 0, 0},   // remainder 1
  {-1,-1, 0, 0, 0, 0, 0, 0},   // remainder 2
  {-1,-1,-1, 0, 0, 0, 0, 0},   // remainder 3
  {-1,-1,-1,-1, 0, 0, 0, 0},   // remainder 4
  {-1,-1,-1,-1,-1, 0, 0, 0},   // remainder 5
  {-1,-1,-1,-1,-1,-1, 0, 0},   // remainder 6
  {-1,-1,-1,-1,-1,-1,-1, 0},   // remainder 7
};


using std::complex;

namespace EVM
{
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
  
  // Returns the highest optimization level as a string
  //    EVM::getOptLevel = string representation of highest optimization level
  std::string getOptLevel();

  //==============//
  // Vector Scale //
  //==============//
  // Multiply each element of a vector by a constant scale factor.  The output
  // buffer may be the same as the input buffer.
  //    in1 = buffer to apply scale factor to
  //    scaleFactor = scale factor to apply to in1
  //    len = number of elements to scale
  //    out = output buffer, product of in1*scaleFactor
  //    overindex_on_read = whether or not input buffers can be over-indexed
  //        when loading into SIMD registers
  void scale(float const * const in1, const float scaleFactor, const int len,
      float * const out);
  void scale(complex<float> const * const in1, const float scaleFactor,
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
  //    overindex_on_read = whether or not input buffers can be over-indexed
  //        when loading into SIMD registers
  void mult(float const * const in1, float const * const in2, const int len,
      float * const out);
  void mult(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  void mult(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  
  // Point by point complex-conjugate multiply of two input buffers, second
  // input conjugated.  The output buffer may be the same as one of the input
  // buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply, conjugated
  //    len = number of elements to multiply
  //    out = output buffer to recieve result of complex multiply
  //    overindex_on_read = whether or not input buffers can be over-indexed
  //        when loading into SIMD registers
  void multc(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  void multc(float const * const in1, complex<float> const * const in2,
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
  //    overindex_on_read = whether or not input buffers can be over-indexed
  //        when loading into SIMD registers
  void mults(float const * const in1, float const * const in2, const int len,
      const float scale, float * const out);
  void mults(complex<float> const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out);
  void mults(float const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out);
  // Point by point complex conjugate multiply of two input buffers with a
  // constant scale factor.  The output buffer may be the same as one of the
  // input buffers.
  //    in1 = first input buffer for multiply
  //    in2 = second input buffer for multiply, conjugated
  //    len = number of elements to multiply
  //    scale = scale factor
  //    out = output buffer to recieve result of complex multiply
  //    overindex_on_read = whether or not input buffers can be over-indexed
  //        when loading into SIMD registers
  void multcs(complex<float> const * const in1,
      complex<float> const * const in2, const int len, const float scale,
      complex<float> * const out);

  //===============//
  // Vector Divide //
  //===============//
  // Point by point multiply of two input buffers.  The output buffer may be
  // the same as one of the input buffers.
  //    in1 = first input buffer for divide (dividend)
  //    in2 = second input buffer for divide (divisor)
  //    len = number of elements to divide
  //    out = output buffer to recieve result of complex division
  void div(float const * const in1, float const * const in2, const int len,
      float * const out);
  void div(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  void div(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  void div(complex<float> const * const in1, float const * const in2,
      const int len, complex<float> * const out);
    
  // Point by point complex-conjugate divide of two input buffers, second
  // input conjugated.  The output buffer may be the same as one of the input
  // buffers.
  //    in1 = first input buffer for divide (dividend)
  //    in2 = second input buffer for divide (divisor), conjugated
  //    len = number of elements to divide
  //    out = output buffer to receive result of complex division
  void divc(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);
  void divc(float const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out);

  //=============//
  // Dot-Product //
  //=============//

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
        std::is_same<T, double>::value ||
        std::is_same<T, float>::value ||
        std::is_same<T, complex<float> >::value ||
        std::is_same<T, complex<double> >::value,
        "[EVM::getSize] Cannot use evm.h with types other than float, double, "
        "complex<float>, or complex<double>");
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
  // implementation for each target if it wasn't explicitly disabled
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  std::string getOptLevel() {return std::string("AVX512");}
  #endif
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2")))
  std::string getOptLevel() {return std::string("AVX2");}
  #endif
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  std::string getOptLevel() {return std::string("AVX");}
  #endif
  #if !defined(DISABLE_SSE4)
  __attribute__((__target__("sse4.2")))
  std::string getOptLevel() {return std::string("SSE4.2");}
  __attribute__((__target__("sse4.1")))
  std::string getOptLevel() {return std::string("SSE4.1");}
  #endif
  #if !defined(DISABLE_SSSE3)
  __attribute__((__target__("ssse3")))
  std::string getOptLevel() {return std::string("SSSE3");}
  #endif
  #if !defined(DISABLE_SSE3)
  __attribute__((__target__("sse3")))
  std::string getOptLevel() {return std::string("SSE3");}
  #endif
  #if !defined(DISABLE_SSE2)
  __attribute__((__target__("sse2")))
  std::string getOptLevel() {return std::string("SSE2");}
  #endif
  #if !defined(DISABLE_SSE)
  __attribute__((__target__("sse")))
  std::string getOptLevel() {return std::string("SSE");}
  #endif
  __attribute__((__target__("default")))
  std::string getOptLevel() {return std::string("SIMD Disabled");}

  //======================================//
  // Vector Scale Function Implementation //
  //======================================//
  // real data, real scale factor
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void scale(float const * const in1, const float scaleFactor, const int len,
      float * const out)
  {
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
      const __mmask16 mk = MASK16(rem);
      ld = _mm512_maskz_loadu_ps (mk, &in1[i]);
      ld = _mm512_mul_ps(ld, sc);
      _mm512_mask_storeu_ps(&out[i], mk, ld);
    }
    return;
  }
  #endif // AVX512 real x real scale
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  void scale(float const * const in1, const float scaleFactor, const int len,
      float * const out)
  {
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
          reinterpret_cast<__m256i const * const>(masks[rem]));
      ld = _mm256_maskload_ps(&in1[i], msk);
      ld = _mm256_mul_ps(ld, sc);
      _mm256_maskstore_ps(&out[i], msk, ld);
    }
    return;
  }
  #endif // AVX real x real scale
  __attribute__((__target__("default")))
  void scale(float const * const in1, const float scaleFactor, const int len,
      float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = scaleFactor*in1[i];
    return;
  }

  // Complex data, real scale factor
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void scale(complex<float> const * const in1, const float scaleFactor,
      const int len, complex<float> * const out)
  {
    __m512 ld;
    __m512 sc = _mm512_set1_ps(scaleFactor);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));
      ld = _mm512_mul_ps(ld, sc);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = MASK16((rem<<1));
      ld = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
        &in1[i]));
      ld = _mm512_mul_ps(ld, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld);
    }
    return;
  }
  #endif // end AVX512 complex x real scale
  #if !defined(DISABLE_AVX) // AVX complex x real scale
  __attribute__((__target__("avx")))
  void scale(complex<float> const * const in1, const float scaleFactor,
      const int len, complex<float> * const out)
  {
    __m256 ld;
    __m256 sc = _mm256_set1_ps(scaleFactor);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld = _mm256_loadu_ps(reinterpret_cast<float const  * const>(&in1[i]));
      ld = _mm256_mul_ps(ld, sc);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld = _mm256_mul_ps(ld, sc);
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld);
    }
    return;
  }
  #endif // AVX512 complex x real scale
  __attribute__((__target__("default")))
  void scale(complex<float> const * const in1, const float scaleFactor,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*scaleFactor;
  }

  // real data, complex scale factor
  #if !defined(DISABLE_AVX512) // AVX512 real data x complex scale
  __attribute__((__target__("avx512f")))
  void scale(float const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out)
  {
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
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i+8]), ld2);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      ld1 = _mm512_maskz_loadu_ps (MASK16(rem), &in1[i]);
      ld2 = _mm512_mul_ps(ld1, re);
      ld3 = _mm512_mul_ps(ld1, im);
      ld1 = _mm512_permutex2var_ps(ld2, p1, ld3);
      ld2 = _mm512_permutex2var_ps(ld2, p2, ld3);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i+8]),
          MASK16(((rem-8)<<1)), ld2);
    }
    else if(rem)
    {
      ld1 = _mm512_maskz_loadu_ps (MASK16(rem), &in1[i]);
      ld2 = _mm512_mul_ps(ld1, re);
      ld3 = _mm512_mul_ps(ld1, im);
      ld1 = _mm512_permutex2var_ps(ld2, p1, ld3);
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]),
          MASK16((rem<<1)), ld1);
    }
    return;
  }
  #endif // end AVX512 real x complex scale
  #if !defined(DISABLE_AVX) // AVX real x complex scale
  __attribute__((__target__("avx")))
  void scale(float const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out)
  {
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
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld3);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i+8]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem>4) // if remainder is > 4, need 2 registers worth
    {
      // 2 floats per complex, so double rem for msk2
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem]));
      const __m256i msk2 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[(rem-4)<<1]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      ld2 = _mm256_mul_ps(ld1, re); // [r0,r1,r2,r3,r4,r5,r6,r7]
      ld3 = _mm256_mul_ps(ld1, im); // [i0,i1,i2,i3,i4.i5,i6,i7]
      ld1 = _mm256_unpacklo_ps(ld2, ld3); // [r0,i0,r1,i1,r4,i4,r5,i5]
      ld2 = _mm256_unpackhi_ps(ld2, ld3); // [r2,i2,r3,i3,r6,i6,r7,i7]
      ld3 = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [r0,i0,r1,i1,r2,i2,r3,i3]
      ld1 = _mm256_permute2f128_ps(ld1, ld2, 0x31); // [r4,i4,r5,i5,r6,i6,r7,i7]
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld3);
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i+4]), msk2,
          ld1);
    }
    else if(rem)
    {
      // 2 floats per complex, so double rem for msk2
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem]));
      const __m256i msk2 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      ld2 = _mm256_mul_ps(ld1, re); // [r0,r1,r2,r3,r4,r5,r6,r7]
      ld3 = _mm256_mul_ps(ld1, im); // [i0,i1,i2,i3,i4.i5,i6,i7]
      ld1 = _mm256_unpacklo_ps(ld2, ld3); // [r0,i0,r1,i1,r4,i4,r5,i5]
      ld2 = _mm256_unpackhi_ps(ld2, ld3); // [r2,i2,r3,i3,r6,i6,r7,i7]
      ld3 = _mm256_permute2f128_ps(ld1, ld2, 0x20); // [r0,i0,r1,i1,r2,i2,r3,i3]
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk2,
          ld3);
    }
    return;
  }
  #endif // end AVX real x complex scale
  __attribute__((__target__("default"))) // default real x complex scale
  void scale(float const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = scaleFactor*in1[i];
    return;
  }

  // Complex data, complex scale factor
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void scale(complex<float> const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out)
  {
    __m512 ld1, ld2;
    __m512 re = _mm512_set1_ps(scaleFactor.real());
    __m512 im = _mm512_set1_ps(scaleFactor.imag());
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));
      ld2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      ld2 = _mm512_mul_ps(ld2, im); // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = MASK16((rem<<1));
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
        &in1[i]));
      ld2 = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      ld2 = _mm512_mul_ps(ld2, im); // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x complex scale
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  void scale(complex<float> const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out)
  {
    __m256 ld1, ld2;
    __m256 re = _mm256_set1_ps(scaleFactor.real());
    __m256 im = _mm256_set1_ps(scaleFactor.imag());
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const  * const>(&in1[i]));
      ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1);
      ld1 = _mm256_mul_ps(ld1, re);
      ld2 = _mm256_mul_ps(ld2, im);
      ld1 = _mm256_addsub_ps(ld1, ld2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld2 = _mm256_shuffle_ps(ld1, ld1, 0xb1);
      ld1 = _mm256_mul_ps(ld1, re);
      ld2 = _mm256_mul_ps(ld2, im);
      ld1 = _mm256_addsub_ps(ld1, ld2);
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX complex x complex scale
  __attribute__((__target__("default")))
  void scale(complex<float> const * const in1, const complex<float> scaleFactor,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*scaleFactor;
    return;
  }

  //=================//
  // Vector Multiply //
  //=================//
  // real x real
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void mult(float const * const in1, float const * const in2, const int len,
      float * const out)
  {
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
      const __mmask16 mk = MASK16(rem);
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
  void mult(float const * const in1, float const * const in2, const int len,
      float * const out)
  {
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
          reinterpret_cast<__m256i const * const>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      ld1 = _mm256_mul_ps(ld1, ld2);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX real x real
  __attribute__((__target__("default")))
  void mult(float const * const in1, float const * const in2, const int len,
    float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i];
    return;
  }

  // complex x complex
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void mult(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    __m512 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = MASK16((rem<<1));
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in1[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in2[i])); // B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x complex
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  void mult(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    __m256 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX complex x complex
  __attribute__((__target__("default")))
  void mult(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i];
    return;
  }
  
  // real x complex
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void mult(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = MASK16(((rem-8)<<1)); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps (MASK16(rem), &in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const * const>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = MASK16((rem<<1)); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps (MASK16(rem), &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const * const>(&in2[i]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 real x complex
  #if !defined(DISABLE_AVX2) // AVX2 real x complex 
  __attribute__((__target__("avx2")))
  void mult(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i+4]));
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const * const>(
            &in2[i+4]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
            msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX2 real x complex
  #if !defined(DISABLE_AVX) // AVX real x complex
  __attribute__((__target__("avx")))
  void mult(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    __m256 ld1, ld2, ld3, sc1, sc2;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i+4]));
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const * const>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const * const>(
            &in2[i+4]), msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
            msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX real x complex
  __attribute__((__target__("default"))) // default real x complex
  void mult(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i];
  } // end default real x complex

  // complex x conj(complex)
  #if !defined(DISABLE_AVX512) // AVX512 complex x conj(complex)
  __attribute__((__target__("avx512f")))
  void multc(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    __m512 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = MASK16((rem<<1));
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in1[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in2[i])); // B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x conj(complex)
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  void multc(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    __m256 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex x conj(complex)
  #if !defined(DISABLE_AVX) // AVX complex x conj(complex)
  __attribute__((__target__("avx")))
  void multc(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    __m256 ld1, ld2, sh, re, im;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
          msk);
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex x conj(complex)
  __attribute__((__target__("default")))
  void multc(complex<float> const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*conj(in2[i]);
    return;
  }
  
  // real x conj(complex)
  #if !defined(DISABLE_AVX512) // AVX512 real x conj(complex)
  __attribute__((__target__("avx512f")))
  void multc(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    __m512 ld1, ld2, ld3, sc;
    const __m512i p1 = _mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7);
    const __m512i p2 = _mm512_setr_epi32(8,8,9,9,10,10,11,11,12,12,13,13,14,14,
        15,15);
    const __m512 neg = _mm512_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0);
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      sc = _mm512_xor_ps(sc, neg);    // negate every other element
      ld1 = _mm512_xor_ps(ld1, neg);  // negate every other element
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = MASK16(((rem-8)<<1)); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps (MASK16(rem), &in1[i]);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const * const>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      sc = _mm512_xor_ps(sc, neg);    // negate every other element
      ld1 = _mm512_xor_ps(ld1, neg);  // negate every other element
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = MASK16((rem<<1)); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps (MASK16(rem), &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const * const>(&in2[i]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      sc = _mm512_xor_ps(sc, neg);    // negate every other element
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 real x conj(complex)
  #if !defined(DISABLE_AVX2) // AVX2 real x conj(complex)
  __attribute__((__target__("avx2")))
  void multc(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    __m256 ld1, ld2, ld3, sc;
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i+4]));
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      sc = _mm256_xor_ps(sc, neg);    // negate every other element
      ld1 = _mm256_xor_ps(ld1, neg);  // negate every other element
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const * const>(
            &in2[i+4]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        sc = _mm256_xor_ps(sc, neg);    // negate every other element
        ld1 = _mm256_xor_ps(ld1, neg);  // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
            msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        sc = _mm256_xor_ps(sc, neg);    // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // AVX2 real x conj(complex)
  #if !defined(DISABLE_AVX) // AVX real x conj(complex)
  __attribute__((__target__("avx")))
  void multc(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    __m256 ld1, ld2, ld3, sc1, sc2;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i+4]));
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      sc1 = _mm256_xor_ps(sc1, neg);      // negate every other element
      sc2 = _mm256_xor_ps(sc2, neg);      // negate every other element
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const * const>(
            &in2[i+4]), msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        sc1 = _mm256_xor_ps(sc1, neg);      // negate every other element
        sc2 = _mm256_xor_ps(sc2, neg);      // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
            msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc1 = _mm256_xor_ps(sc1, neg);        // negate every other element
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX real x conj(complex)
  __attribute__((__target__("default"))) // default real x conj(complex)
  void multc(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*conj(in2[i]);
    return;
  }
  
  //============================//
  // Vectory Multiply and Scale //
  //============================//
  // real x real x scale
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void mults(float const * const in1, float const * const in2, const int len,
      const float scale, float * const out)
  {
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
      const __mmask16 mk = MASK16(rem);
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
  void mults(float const * const in1, float const * const in2, const int len,
      const float scale, float * const out)
  {
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
          reinterpret_cast<__m256i const * const>(masks[rem]));
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
  void mults(float const * const in1, float const * const in2, const int len,
    const float scale, float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i]*scale;
  }

  // complex x complex x scale
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void mults(complex<float> const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    __m512 ld1, ld2, sh, re, im;
    const __m512 sc = _mm512_set1_ps(scale);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      ld2 = _mm512_mul_ps(ld2, sc);
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = MASK16((rem<<1));
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in1[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in2[i])); // B
      ld2 = _mm512_mul_ps(ld2, sc);
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmaddsub_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x complex x scale
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  void mults(complex<float> const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    __m256 ld1, ld2, sh, re, im;
    const __m256 sc = _mm256_set1_ps(scale);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      ld2 = _mm256_mul_ps(ld2, sc);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
          msk);
      ld2 = _mm256_mul_ps(ld2, sc);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // AVX complex x complex x scale
  __attribute__((__target__("default"))) // default complex x complex x scale
  void mults(complex<float> const * const in1, complex<float> const * const in2,
    const int len, const float scale, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i]*scale;
  }

  // real x complex x scale
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void mults(float const * const in1, complex<float> const * const in2,
    const int len, const float scale, complex<float> * const out)
  {
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
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i+8]), ld3);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem>8) // if remainder is > 8, need 2 registers worth
    {
      const __mmask16 mk = MASK16(((rem-8)<<1)); // 2 floats per complex
      ld1 = _mm512_maskz_loadu_ps (MASK16(rem), &in1[i]);
      ld1 = _mm512_mul_ps(ld1, sca);
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const * const>(&in2[i+8]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld1 = _mm512_permutexvar_ps(p2, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      ld3 = _mm512_mul_ps(ld3, ld1);
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i+8]), mk,
          ld3);
    }
    else if(rem)
    {
      const __mmask16 mk = MASK16((rem<<1)); // 2 floats per complex 
      ld1 = _mm512_maskz_loadu_ps (MASK16(rem), &in1[i]);
      ld1 = _mm512_mul_ps(ld1, sca);
      ld2 = _mm512_maskz_loadu_ps(mk,
          reinterpret_cast<float const * const>(&in2[i]));
      sc = _mm512_permutexvar_ps(p1, ld1);
      ld2 = _mm512_mul_ps(ld2, sc);
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld2);
    }
    return;
  }
  #endif // end AVX512 real x complex x scale
  #if !defined(DISABLE_AVX2) // AVX2 real x complex x scale
  __attribute__((__target__("avx2")))
  void mults(float const * const in1, complex<float> const * const in2,
    const int len, const float scale, complex<float> * const out)
  {
    __m256 ld1, ld2, ld3, sc;
    const __m256 sca = _mm256_set1_ps(scale);
    const __m256i p1 = _mm256_setr_epi32(0,0,1,1,2,2,3,3);
    const __m256i p2 = _mm256_setr_epi32(4,4,5,5,6,6,7,7);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld1 = _mm256_mul_ps(ld1, sca);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i+4]));
      sc = _mm256_permutevar8x32_ps(ld1, p1);
      ld1 = _mm256_permutevar8x32_ps(ld1, p2);
      ld2 = _mm256_mul_ps(ld2, sc);
      ld3 = _mm256_mul_ps(ld3, ld1);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      ld1 = _mm256_mul_ps(ld1, sca);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const * const>(
            &in2[i+4]), msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld1 = _mm256_permutevar8x32_ps(ld1, p2);
        ld2 = _mm256_mul_ps(ld2, sc);
        ld3 = _mm256_mul_ps(ld3, ld1);
        _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
            msk2);
        sc = _mm256_permutevar8x32_ps(ld1, p1);
        ld2 = _mm256_mul_ps(ld2, sc);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX2 real x complex x scale
  #if !defined(DISABLE_AVX) // AVX real x complex x scale
  __attribute__((__target__("avx")))
  void mults(float const * const in1, complex<float> const * const in2,
    const int len, const float scale, complex<float> * const out)
  {
    __m256 ld1, ld2, ld3, sc1, sc2;
    const __m256 sca = _mm256_set1_ps(scale);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 real elements per register
    {
      ld1 = _mm256_loadu_ps(&in1[i]);
      ld1 = _mm256_mul_ps(ld1, sca);
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
      ld3 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i+4]));
      sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
      ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
      sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
      sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
      ld2 = _mm256_mul_ps(ld2, sc1);
      ld3 = _mm256_mul_ps(ld3, sc2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i+4]), ld3);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // note msk2 accounts for 2 reals per element for the complex buffer
      const __m256i msk1 = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk1);
      ld1 = _mm256_mul_ps(ld1, sca);
      if(rem>4) // if remainder is > 4, need 2 registers worth
      {
        const __m256i msk2 = _mm256_load_si256(
           reinterpret_cast<__m256i const * const>(masks[(rem-4)<<1]));
        ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));
        ld3 = _mm256_maskload_ps(reinterpret_cast<float const * const>(
            &in2[i+4]), msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        sc2 = _mm256_permute2f128_ps(sc2, ld1, 0x31); // [4,4,5,5,6,6,7,7]
        ld2 = _mm256_mul_ps(ld2, sc1);
        ld3 = _mm256_mul_ps(ld3, sc2);
        _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld2);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i+4]), msk2,
            ld3);
      }
      else
      {
        const __m256i msk2 = _mm256_load_si256(
            reinterpret_cast<__m256i const * const>(masks[rem<<1]));
        ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
            msk2);
        sc2 = _mm256_unpacklo_ps(ld1, ld1); // [0,0,1,1,4,4,5,5]
        ld1 = _mm256_unpackhi_ps(ld1, ld1); // [2,2,3,3,6,6,7,7]
        sc1 = _mm256_permute2f128_ps(sc2, ld1, 0x20); // [0,0,1,1,2,2,3,3]
        ld2 = _mm256_mul_ps(ld2, sc1);
        _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk2,
            ld2);
      }
    }
    return;
  }
  #endif // end AVX real x complex x scale
  __attribute__((__target__("default"))) // default real x complex x scale
  void mults(float const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*in2[i]*scale;
  }

  // complex x conj(complex) x scale
  #if !defined(DISABLE_AVX512) // AVX512 complex x conj(complex) x scale
  __attribute__((__target__("avx512f")))
  void multcs(complex<float> const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    __m512 ld1, ld2, sh, re, im;
    const __m512 sc = _mm512_set1_ps(scale);
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      ld2 = _mm512_mul_ps(ld2, sc);
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = MASK16((rem<<1));
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in1[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in2[i])); // B
      ld2 = _mm512_mul_ps(ld2, sc);
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex x conj(complex) x scale
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  void multcs(complex<float> const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    __m256 ld1, ld2, sh, re, im;
    const __m256 sc = _mm256_set1_ps(scale);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      ld2 = _mm256_mul_ps(ld2, sc);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
          msk);
      ld2 = _mm256_mul_ps(ld2, sc);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex x conj(complex) x scale
  #if !defined(DISABLE_AVX) // AVX complex x conj(complex) x scale
  __attribute__((__target__("avx")))
  void multcs(complex<float> const * const in1, complex<float> const * const in2,
      const int len, const float scale, complex<float> * const out)
  {
    __m256 ld1, ld2, sh, re, im;
    const __m256 sc = _mm256_setr_ps(scale, -scale, scale, -scale, scale,
        -scale, scale, -scale);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm256_mul_ps(ld2, sc);// conj(B) x scale
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
          msk);
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld2 = _mm256_mul_ps(ld2, sc);// conj(B) x scale
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex x conj(complex) x scale
  // default complex x conj(complex) x scale
  __attribute__((__target__("default")))
  void multcs(float const * const in1, complex<float> const * const in2,
    const int len, const float scale, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]*conj(in2[i])*scale;
  }

  //===============//
  // Vector Divide //
  //===============//
  // real / real
  #if !defined(DISABLE_AVX512)
  __attribute__((__target__("avx512f")))
  void div(float const * const in1, float const * const in2, const int len,
      float * const out)
  {
    __m512 ld1, ld2;
    int i = 0;
    for(; i<len-15; i+=16) // process 16 real elements per register
    {
      ld1 = _mm512_loadu_ps(&in1[i]);
      ld2 = _mm512_loadu_ps(&in2[i]);
      ld1 = _mm512_div_ps(ld1, ld2);
      _mm512_storeu_ps(&out[i], ld1);
    }
    // handle remaining elements (note len&15 == len%16)
    const int rem = len&15;
    if(rem)
    {
      const __mmask16 mk = MASK16(rem);
      ld1 = _mm512_maskz_loadu_ps(mk, &in1[i]);
      ld2 = _mm512_maskz_loadu_ps(mk, &in2[i]);
      ld1 = _mm512_div_ps(ld1, ld2);
      _mm512_mask_storeu_ps(&out[i], mk, ld1);
    }
    return;
  }
  #endif // AVX512 real / real
  #if !defined(DISABLE_AVX)
  __attribute__((__target__("avx")))
  void div(float const * const in1, float const * const in2, const int len,
      float * const out)
  {
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
          reinterpret_cast<__m256i const * const>(masks[rem]));
      ld1 = _mm256_maskload_ps(&in1[i], msk);
      ld2 = _mm256_maskload_ps(&in2[i], msk);
      ld1 = _mm256_div_ps(ld1, ld2);
      _mm256_maskstore_ps(&out[i], msk, ld1);
    }
    return;
  }
  #endif // AVX real / real
  __attribute__((__target__("default"))) // default real / real
  void div(float const * const in1, float const * const in2, const int len,
    float * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }

  // complex / complex
  #if !defined(DISABLE_AVX512) // AVX512 complex / complex
  __attribute__((__target__("avx512f")))
  void div(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    __m512 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-7; i+=8) // process 8 complex elements per register
    {
      ld1 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm512_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&7 == len%8)
    const int rem = len&7;
    if(rem)
    {
      // each complex is 2 floats, so double rem
      const __mmask16 mk = MASK16((rem<<1));
      ld1 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in1[i])); // A
      ld2 = _mm512_maskz_loadu_ps(mk, reinterpret_cast<float const * const>(
          &in2[i])); // B
      sh = _mm512_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm512_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm512_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm512_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm512_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm512_mask_storeu_ps(reinterpret_cast<float * const>(&out[i]), mk, ld1);
    }
    return;
  }
  #endif // AVX512 complex / complex
  #if !defined(DISABLE_AVX2)
  __attribute__((__target__("avx2,fma")))
  void div(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    __m256 ld1, ld2, sh, re, im;
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
          msk);
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai7,Ar7]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi7,Bi7]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br7,Br7]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai7*Bi7,Ar7*Bi7]
      ld1 = _mm256_fmsubadd_ps(re, ld1, ld2);// [Br0*Ar0-Bi0Ai0,Br0*Ai0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX2 complex / complex
  #if !defined(DISABLE_AVX) // AVX complex / complex
  __attribute__((__target__("avx")))
  void div(complex<float> const * const in1, complex<float> const * const in2,
      const int len, complex<float> * const out)
  {
    __m256 ld1, ld2, sh, re, im;
    const __m256 neg = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f,
        0.0f, -0.0f);
    int i = 0;
    for(; i<len-3; i+=4) // process 4 complex elements per register
    {
      ld2 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in2[i]));// B
      ld1 = _mm256_loadu_ps(reinterpret_cast<float const * const>(&in1[i]));// A
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_storeu_ps(reinterpret_cast<float * const>(&out[i]), ld1);
    }
    // handle remaining elements (note len&3 == len%4)
    const int rem = len&3;
    if(rem)
    {
      // 2 floats per complex, so double rem for mask
      const __m256i msk = _mm256_load_si256(
          reinterpret_cast<__m256i const * const>(masks[rem<<1]));
      ld2 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in2[i]),
          msk);
      ld1 = _mm256_maskload_ps(reinterpret_cast<float const * const>(&in1[i]),
          msk);
      ld2 = _mm256_xor_ps(ld2, neg);// conj(B)
      sh = _mm256_shuffle_ps(ld1, ld1, 0xb1); // [Ai0,Ar0,Ai1,Ar1,...,Ai3,Ar3]
      im = _mm256_movehdup_ps(ld2); // [Bi0,Bi0,Bi1,Bi1,...,Bi3,Bi3]
      re = _mm256_moveldup_ps(ld2); // [Br0,Br0,Br1,Br1,...,Br3,Br3]
      ld2 = _mm256_mul_ps(sh, im);  // [Ai0*Bi0,Ar0*Bi0,...,Ai3*Bi3,Ar3*Bi3]
      ld1 = _mm256_mul_ps(ld1, re); // [Ar0*Br0,Ai0*Br0,...,Ar3*Br3,Ai3*Br3]
      ld1 = _mm256_addsub_ps(ld1, ld2);// [Ar0*Br0-Ai0*Bi0,Ai0*Br0+Ar0*Bi0]
      _mm256_maskstore_ps(reinterpret_cast<float * const>(&out[i]), msk, ld1);
    }
    return;
  }
  #endif // end AVX complex / complex
  __attribute__((__target__("default"))) // default complex / complex
  void div(complex<float> const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {
    for(int i=0; i<len; ++i) out[i] = in1[i]/in2[i];
  }

  
  __attribute__((__target__("default"))) // default real / complex
  void div(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {

  }

  __attribute__((__target__("default"))) // default complex / real
  void div(complex<float> const * const in1, float const * const in2,
    const int len, complex<float> * const out)
  {

  }
  
  __attribute__((__target__("default"))) // default complex / conj(complex)
  void divc(complex<float> const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {

  }

  __attribute__((__target__("default"))) // default real / conj(complex)
  void divc(float const * const in1, complex<float> const * const in2,
    const int len, complex<float> * const out)
  {

  }

  
} // namespace EVM
#endif // EFFICIENT_VECTOR_MATH_H