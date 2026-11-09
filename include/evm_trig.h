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
// Name:     evm_trig.h
//
// Purpose:  SIMD sine, cosine and atan2 for evm.h that don't need SVML
//
// Created:  2026/11/05
//
// Description:
//  Inline helpers that compute the sine and cosine, or the atan2, of a whole
//  register at once, for the AVX512 (avx512f), AVX2 (avx2,fma) and AVX (avx)
//  versions of the evm.h functions.  Each helper has the target attribute of
//  the instruction set it is written for, so it can only be called from a
//  function with the same (or a larger) target:
//    sincosCyc  sine and cosine of 2*pi*ph, phase in cycles
//    sincosRad  sine and cosine of x, x in radians
//    atan2Cyc   atan2(y, x) in cycles
//    atan2Rad   atan2(y, x) in radians
//  with the suffix 512, Avx2 or Avx.
//
//  The polynomials are the Cephes single precision minimax ones (Stephen
//  Moshier), the cycles versions fold 2*pi (sin/cos) or 1/(2*pi) (atan) into
//  the coefficients.  Accuracy (maximum error measured over random and edge
//  case inputs):
//    sincosCyc  1.7 ulp at any phase: whole cycles and quarter cycles are
//               removed exactly, unlike sincos(2*pi*ph), which loses accuracy
//               as the phase grows (abs error 8e-7 by 2 cycles, 4e-3 by 1e4);
//               the AVX2/AVX versions give NaN past |ph| = 8.5e37
//    sincosRad  1.6 ulp (AVX512/AVX2) or 3.3 ulp (AVX) for |x| <= 65536,
//               larger |x| fall back to scalar sincosf for those lanes
//    atan2      2.5 ulp (radians), 2.3 ulp (cycles)
//  Zeros (with their signs), infinities and NaN give the same results as the
//  C library's atan2f, sinf and cosf.
//==============================================================================
#ifndef EVM_TRIG_H
#define EVM_TRIG_H

#include <cmath>     // for sincosf
#include <math.h>
#include <x86intrin.h>

namespace EVM
{
namespace TRIG
{
  // the radians reduction stays accurate up to this |x|, larger values fall
  // back to scalar sincosf
  constexpr float SINCOS_RAD_LIMIT = 65536.0f;

  #define TA5 __attribute__((__target__("avx512f"), __always_inline__))
  #define TA2 __attribute__((__target__("avx2,fma"), __always_inline__))
  #define TA1 __attribute__((__target__("avx"), __always_inline__))

  // sincosf for the lanes set in k, the other lanes keep sn and cs
  __attribute__((__target__("avx512f"), __noinline__))
  inline void sincosRadFar512(const __m512 x, const __mmask16 k, __m512 &sn,
      __m512 &cs)
  {
    alignas(64) float xv[16], sv[16], cv[16];
    _mm512_store_ps(xv, x);
    _mm512_store_ps(sv, sn);
    _mm512_store_ps(cv, cs);
    for(int i=0; i<16; ++i) if((k>>i)&1) sincosf(xv[i], &sv[i], &cv[i]);
    sn = _mm512_load_ps(sv);
    cs = _mm512_load_ps(cv);
  }
  __attribute__((__target__("avx"), __noinline__))
  inline void sincosRadFar256(const __m256 x, const int k, __m256 &sn,
      __m256 &cs)
  {
    alignas(32) float xv[8], sv[8], cv[8];
    _mm256_store_ps(xv, x);
    _mm256_store_ps(sv, sn);
    _mm256_store_ps(cv, cs);
    for(int i=0; i<8; ++i) if((k>>i)&1) sincosf(xv[i], &sv[i], &cv[i]);
    sn = _mm256_load_ps(sv);
    cs = _mm256_load_ps(cv);
  }

  //=========================================================== AVX-512 ====//
  // quadrant fix-up: q odd swaps sin and cos, sin is negated for q = 2,3 and
  // cos for q = 1,2 (mod 4), i.e. bit 1 of q and of q+1
  TA5 inline void quadrant512(const __m512i qi, const __m512 s,
      const __m512 c, __m512 &sn, __m512 &cs)
  {
    const __m512i sgn = _mm512_set1_epi32(0x80000000);
    const __mmask16 odd = _mm512_test_epi32_mask(qi, _mm512_set1_epi32(1));
    const __m512i fs = _mm512_and_epi32(_mm512_slli_epi32(qi, 30), sgn);
    const __m512i fc = _mm512_and_epi32(_mm512_slli_epi32(
        _mm512_add_epi32(qi, _mm512_set1_epi32(1)), 30), sgn);
    sn = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(
        _mm512_mask_blend_ps(odd, s, c)), fs));
    cs = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(
        _mm512_mask_blend_ps(odd, c, s)), fc));
  }

  TA5 inline void sincosCyc512(const __m512 ph, __m512 &sn, __m512 &cs)
  {
    // q4 = ph rounded to the nearest quarter cycle (roundscale keeping 2
    // fraction bits, so nothing can overflow), ph - q4 is exact (|t| <= 1/8),
    // so the reduction adds no error however large ph is
    const __m512 q4 = _mm512_roundscale_ps(ph, (2<<4) |
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    const __m512 t = _mm512_sub_ps(ph, q4);
    const __m512 z = _mm512_mul_ps(t, t);
    __m512 s = _mm512_fmadd_ps(z, _mm512_set1_ps(-75.4456558f),
        _mm512_set1_ps(81.5937653f));
    s = _mm512_fmadd_ps(z, s, _mm512_set1_ps(-41.341671f));
    s = _mm512_fmadd_ps(z, s, _mm512_set1_ps(6.28318548f));
    s = _mm512_mul_ps(t, s);                         // sin(2*pi*t)
    __m512 c = _mm512_fmadd_ps(z, _mm512_set1_ps(59.3497009f),
        _mm512_set1_ps(-85.4471436f));
    c = _mm512_fmadd_ps(z, c, _mm512_set1_ps(64.9393616f));
    c = _mm512_fmadd_ps(z, c, _mm512_set1_ps(-19.7392082f));
    c = _mm512_fmadd_ps(z, c, _mm512_set1_ps(1.0f)); // cos(2*pi*t)
    // quadrant from 4*q4 (for |ph| >= 2^29 the conversion overflows to
    // 0x80000000, quadrant 0, which is right for a whole number of cycles)
    quadrant512(_mm512_cvtps_epi32(_mm512_mul_ps(q4, _mm512_set1_ps(4.0f))),
        s, c, sn, cs);
  }

  TA5 inline void sincosRad512(const __m512 x, __m512 &sn, __m512 &cs)
  {
    // r = x - q*pi/2 with pi/2 in three pieces (Cody-Waite), |r| <~ pi/4
    const __m512 q = _mm512_roundscale_ps(_mm512_mul_ps(x,
        _mm512_set1_ps(0.636619747f)),
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m512 r = _mm512_fnmadd_ps(q, _mm512_set1_ps(1.57079637f), x);
    r = _mm512_fnmadd_ps(q, _mm512_set1_ps(-4.37113883e-08f), r);
    r = _mm512_fnmadd_ps(q, _mm512_set1_ps(-1.71512451e-15f), r);
    const __m512 z = _mm512_mul_ps(r, r);
    __m512 s = _mm512_fmadd_ps(z, _mm512_set1_ps(-0.000195152956f),
        _mm512_set1_ps(0.00833216123f));
    s = _mm512_fmadd_ps(z, s, _mm512_set1_ps(-0.166666552f));
    s = _mm512_fmadd_ps(_mm512_mul_ps(z, r), s, r);  // r + r^3*P(z)
    __m512 c = _mm512_fmadd_ps(z, _mm512_set1_ps(2.44331568e-05f),
        _mm512_set1_ps(-0.00138873165f));
    c = _mm512_fmadd_ps(z, c, _mm512_set1_ps(0.0416666456f));
    c = _mm512_fmadd_ps(_mm512_mul_ps(z, z), c,      // 1 - z/2 + z^2*Q(z)
        _mm512_fnmadd_ps(z, _mm512_set1_ps(0.5f), _mm512_set1_ps(1.0f)));
    quadrant512(_mm512_cvtps_epi32(q), s, c, sn, cs);
    // lanes past the accurate reduction range (rare) go through sincosf
    const __mmask16 far = _mm512_cmp_ps_mask(_mm512_abs_ps(x),
        _mm512_set1_ps(SINCOS_RAD_LIMIT), _CMP_GT_OQ);
    if(far) sincosRadFar512(x, far, sn, cs);
  }

  // atan2 with one division: with mn = min(|x|,|y|), mx = max(|x|,|y|),
  // u = mn/mx, or (mn-mx)/(mn+mx) plus pi/4 when mn/mx > tan(pi/8), so the
  // polynomial only sees |u| <= tan(pi/8)
  template<bool CYC>
  TA5 inline __m512 atan2_512(const __m512 y, const __m512 x)
  {
    const __m512i noSign = _mm512_set1_epi32(0x7fffffff);
    const __m512 ax = _mm512_castsi512_ps(_mm512_and_epi32(
        _mm512_castps_si512(x), noSign));
    const __m512 ay = _mm512_castsi512_ps(_mm512_and_epi32(
        _mm512_castps_si512(y), noSign));
    const __mmask16 sw = _mm512_cmp_ps_mask(ay, ax, _CMP_GT_OQ);
    const __m512 mn = _mm512_mask_blend_ps(sw, ay, ax);
    const __m512 mx = _mm512_mask_blend_ps(sw, ax, ay);
    __mmask16 big = _mm512_cmp_ps_mask(mn, _mm512_mul_ps(mx,
        _mm512_set1_ps(0.414213568f)), _CMP_GT_OQ);
    // both infinite: angle is the pi/4 base, both zero: angle is 0
    const __mmask16 inf = _mm512_cmp_ps_mask(mn,
        _mm512_set1_ps(__builtin_inff()), _CMP_EQ_OQ);
    const __mmask16 zero = _mm512_cmp_ps_mask(mx, _mm512_setzero_ps(),
        _CMP_EQ_OQ);
    __m512 n = _mm512_mask_sub_ps(mn, big, mn, mx);
    __m512 d = _mm512_mask_add_ps(mx, big, mn, mx);
    n = _mm512_mask_mov_ps(n, inf, _mm512_setzero_ps());
    d = _mm512_mask_mov_ps(d, inf | zero, _mm512_set1_ps(1.0f));
    big |= inf;
    const __m512 u = _mm512_div_ps(n, d);
    const __m512 z = _mm512_mul_ps(u, u);
    __m512 a;
    if(CYC) // polynomial scaled by 1/(2*pi), quadrant constants exact
    {
      __m512 p = _mm512_fmadd_ps(z, _mm512_set1_ps(0.0128179323f),
          _mm512_set1_ps(-0.0220870227f));
      p = _mm512_fmadd_ps(z, p, _mm512_set1_ps(0.0317955129f));
      p = _mm512_fmadd_ps(z, p, _mm512_set1_ps(-0.0530510359f));
      // u/(2*pi) with 1/(2*pi) as hi + lo
      a = _mm512_fmadd_ps(u, _mm512_set1_ps(0.159154937f), _mm512_fmadd_ps(u,
          _mm512_set1_ps(6.42063824e-09f), _mm512_mul_ps(_mm512_mul_ps(u, z),
          p)));
      a = _mm512_mask_add_ps(a, big, a, _mm512_set1_ps(0.125f));
      a = _mm512_mask_sub_ps(a, sw, _mm512_set1_ps(0.25f), a);
    }
    else
    {
      __m512 p = _mm512_fmadd_ps(z, _mm512_set1_ps(0.0805374458f),
          _mm512_set1_ps(-0.138776854f));
      p = _mm512_fmadd_ps(z, p, _mm512_set1_ps(0.199777111f));
      p = _mm512_fmadd_ps(z, p, _mm512_set1_ps(-0.333329499f));
      a = _mm512_fmadd_ps(_mm512_mul_ps(u, z), p, u);
      // pi/4, pi/2 and pi as hi + lo so the constants add no rounding error
      a = _mm512_mask_add_ps(a, big, _mm512_add_ps(a,
          _mm512_set1_ps(-2.18556941e-08f)), _mm512_set1_ps(0.785398185f));
      a = _mm512_mask_add_ps(a, sw, _mm512_sub_ps(_mm512_set1_ps(1.57079637f),
          a), _mm512_set1_ps(-4.37113883e-08f));
    }
    // x negative (sign bit, so -0 counts): pi - a, then the sign of y
    const __mmask16 xneg = _mm512_cmplt_epi32_mask(_mm512_castps_si512(x),
        _mm512_setzero_si512());
    if(CYC) a = _mm512_mask_sub_ps(a, xneg, _mm512_set1_ps(0.5f), a);
    else a = _mm512_mask_add_ps(a, xneg, _mm512_sub_ps(
        _mm512_set1_ps(3.14159274f), a), _mm512_set1_ps(-8.74227766e-08f));
    a = _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(a),
        _mm512_andnot_epi32(noSign, _mm512_castps_si512(y))));
    // NaN in, NaN out
    return _mm512_mask_add_ps(a, _mm512_cmp_ps_mask(x, y, _CMP_UNORD_Q), x, y);
  }

  //============================================================== AVX2 ====//
  // m ? b : a without blendv (GCC 12 can turn blendv into branches)
  __attribute__((__target__("avx"), __always_inline__))
  inline __m256 sel256(const __m256 m, const __m256 a, const __m256 b)
  {
    return _mm256_or_ps(_mm256_and_ps(m, b), _mm256_andnot_ps(m, a));
  }

  TA2 inline void quadrantAvx2(const __m256i qi, const __m256 s,
      const __m256 c, __m256 &sn, __m256 &cs)
  {
    const __m256i sgn = _mm256_set1_epi32(0x80000000);
    const __m256i one = _mm256_set1_epi32(1);
    const __m256 odd = _mm256_castsi256_ps(_mm256_cmpeq_epi32(
        _mm256_and_si256(qi, one), one));
    const __m256 fs = _mm256_castsi256_ps(_mm256_and_si256(
        _mm256_slli_epi32(qi, 30), sgn));
    const __m256 fc = _mm256_castsi256_ps(_mm256_and_si256(
        _mm256_slli_epi32(_mm256_add_epi32(qi, one), 30), sgn));
    sn = _mm256_xor_ps(sel256(odd, s, c), fs);
    cs = _mm256_xor_ps(sel256(odd, c, s), fc);
  }

  TA2 inline void sincosCycAvx2(const __m256 ph, __m256 &sn, __m256 &cs)
  {
    // q = nearest quarter cycle, ph - q/4 is exact (|t| <= 1/8), so the
    // reduction adds no error however large ph is, up to |ph| = 8.5e37 where
    // 4*ph overflows and the result is NaN (removing whole cycles first would
    // avoid that, but costs 13-20%)
    const __m256 q = _mm256_round_ps(_mm256_mul_ps(ph, _mm256_set1_ps(4.0f)),
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    const __m256 t = _mm256_fnmadd_ps(q, _mm256_set1_ps(0.25f), ph);
    const __m256 z = _mm256_mul_ps(t, t);
    __m256 s = _mm256_fmadd_ps(z, _mm256_set1_ps(-75.4456558f),
        _mm256_set1_ps(81.5937653f));
    s = _mm256_fmadd_ps(z, s, _mm256_set1_ps(-41.341671f));
    s = _mm256_fmadd_ps(z, s, _mm256_set1_ps(6.28318548f));
    s = _mm256_mul_ps(t, s);
    __m256 c = _mm256_fmadd_ps(z, _mm256_set1_ps(59.3497009f),
        _mm256_set1_ps(-85.4471436f));
    c = _mm256_fmadd_ps(z, c, _mm256_set1_ps(64.9393616f));
    c = _mm256_fmadd_ps(z, c, _mm256_set1_ps(-19.7392082f));
    c = _mm256_fmadd_ps(z, c, _mm256_set1_ps(1.0f));
    quadrantAvx2(_mm256_cvtps_epi32(q), s, c, sn, cs);
  }

  TA2 inline void sincosRadAvx2(const __m256 x, __m256 &sn, __m256 &cs)
  {
    const __m256 q = _mm256_round_ps(_mm256_mul_ps(x,
        _mm256_set1_ps(0.636619747f)),
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m256 r = _mm256_fnmadd_ps(q, _mm256_set1_ps(1.57079637f), x);
    r = _mm256_fnmadd_ps(q, _mm256_set1_ps(-4.37113883e-08f), r);
    r = _mm256_fnmadd_ps(q, _mm256_set1_ps(-1.71512451e-15f), r);
    const __m256 z = _mm256_mul_ps(r, r);
    __m256 s = _mm256_fmadd_ps(z, _mm256_set1_ps(-0.000195152956f),
        _mm256_set1_ps(0.00833216123f));
    s = _mm256_fmadd_ps(z, s, _mm256_set1_ps(-0.166666552f));
    s = _mm256_fmadd_ps(_mm256_mul_ps(z, r), s, r);
    __m256 c = _mm256_fmadd_ps(z, _mm256_set1_ps(2.44331568e-05f),
        _mm256_set1_ps(-0.00138873165f));
    c = _mm256_fmadd_ps(z, c, _mm256_set1_ps(0.0416666456f));
    c = _mm256_fmadd_ps(_mm256_mul_ps(z, z), c,
        _mm256_fnmadd_ps(z, _mm256_set1_ps(0.5f), _mm256_set1_ps(1.0f)));
    quadrantAvx2(_mm256_cvtps_epi32(q), s, c, sn, cs);
    // lanes past the accurate reduction range (rare) go through sincosf
    const int far = _mm256_movemask_ps(_mm256_cmp_ps(_mm256_andnot_ps(
        _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)), x),
        _mm256_set1_ps(SINCOS_RAD_LIMIT), _CMP_GT_OQ));
    if(far) sincosRadFar256(x, far, sn, cs);
  }

  template<bool CYC>
  TA2 inline __m256 atan2Avx2(const __m256 y, const __m256 x)
  {
    const __m256 noSign = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    const __m256 ax = _mm256_and_ps(x, noSign);
    const __m256 ay = _mm256_and_ps(y, noSign);
    const __m256 sw = _mm256_cmp_ps(ay, ax, _CMP_GT_OQ);
    const __m256 mn = sel256(sw, ay, ax);
    const __m256 mx = sel256(sw, ax, ay);
    __m256 big = _mm256_cmp_ps(mn, _mm256_mul_ps(mx,
        _mm256_set1_ps(0.414213568f)), _CMP_GT_OQ);
    const __m256 inf = _mm256_cmp_ps(mn, _mm256_set1_ps(__builtin_inff()),
        _CMP_EQ_OQ);
    const __m256 zero = _mm256_cmp_ps(mx, _mm256_setzero_ps(), _CMP_EQ_OQ);
    __m256 n = sel256(big, mn, _mm256_sub_ps(mn, mx));
    __m256 d = sel256(big, mx, _mm256_add_ps(mn, mx));
    n = _mm256_andnot_ps(inf, n);
    d = sel256(_mm256_or_ps(inf, zero), d, _mm256_set1_ps(1.0f));
    big = _mm256_or_ps(big, inf);
    const __m256 u = _mm256_div_ps(n, d);
    const __m256 z = _mm256_mul_ps(u, u);
    __m256 a;
    if(CYC)
    {
      __m256 p = _mm256_fmadd_ps(z, _mm256_set1_ps(0.0128179323f),
          _mm256_set1_ps(-0.0220870227f));
      p = _mm256_fmadd_ps(z, p, _mm256_set1_ps(0.0317955129f));
      p = _mm256_fmadd_ps(z, p, _mm256_set1_ps(-0.0530510359f));
      a = _mm256_fmadd_ps(u, _mm256_set1_ps(0.159154937f), _mm256_fmadd_ps(u,
          _mm256_set1_ps(6.42063824e-09f), _mm256_mul_ps(_mm256_mul_ps(u, z),
          p)));
      a = _mm256_add_ps(a, _mm256_and_ps(big, _mm256_set1_ps(0.125f)));
      a = sel256(sw, a, _mm256_sub_ps(_mm256_set1_ps(0.25f), a));
      a = sel256(_mm256_castsi256_ps(_mm256_srai_epi32(
          _mm256_castps_si256(x), 31)), a, _mm256_sub_ps(
          _mm256_set1_ps(0.5f), a));
    }
    else
    {
      __m256 p = _mm256_fmadd_ps(z, _mm256_set1_ps(0.0805374458f),
          _mm256_set1_ps(-0.138776854f));
      p = _mm256_fmadd_ps(z, p, _mm256_set1_ps(0.199777111f));
      p = _mm256_fmadd_ps(z, p, _mm256_set1_ps(-0.333329499f));
      a = _mm256_fmadd_ps(_mm256_mul_ps(u, z), p, u);
      // pi/4, pi/2 and pi as hi + lo so the constants add no rounding error
      a = _mm256_add_ps(_mm256_add_ps(a, _mm256_and_ps(big,
          _mm256_set1_ps(-2.18556941e-08f))), _mm256_and_ps(big,
          _mm256_set1_ps(0.785398185f)));
      a = sel256(sw, a, _mm256_add_ps(_mm256_sub_ps(
          _mm256_set1_ps(1.57079637f), a), _mm256_set1_ps(-4.37113883e-08f)));
      a = sel256(_mm256_castsi256_ps(_mm256_srai_epi32(
          _mm256_castps_si256(x), 31)), a, _mm256_add_ps(_mm256_sub_ps(
          _mm256_set1_ps(3.14159274f), a), _mm256_set1_ps(-8.74227766e-08f)));
    }
    a = _mm256_or_ps(a, _mm256_andnot_ps(noSign, y));
    const __m256 nan = _mm256_cmp_ps(x, y, _CMP_UNORD_Q);
    return sel256(nan, a, _mm256_add_ps(x, y));
  }

  //=============================================================== AVX ====//
  // no FMA and no 256-bit integer instructions: the quadrant logic works on
  // q as a float, and the radians reduction uses pi/2 in five pieces whose
  // first four have 8 bits, so q*piece is exact for |q| < 65536
  TA1 inline void quadrantAvx(const __m256 q, const __m256 s, const __m256 c,
      __m256 &sn, __m256 &cs)
  {
    const __m256 sgn = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
    // q mod 4 in -2..2 (exact for any integer valued q)
    const __m256 m = _mm256_sub_ps(q, _mm256_mul_ps(_mm256_round_ps(
        _mm256_mul_ps(q, _mm256_set1_ps(0.25f)),
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC), _mm256_set1_ps(4.0f)));
    const __m256 odd = _mm256_cmp_ps(_mm256_andnot_ps(sgn, m),
        _mm256_set1_ps(1.0f), _CMP_EQ_OQ);
    const __m256 fs = _mm256_and_ps(sgn, _mm256_or_ps(       // m = 2,-1,-2
        _mm256_cmp_ps(m, _mm256_set1_ps(1.5f), _CMP_GT_OQ),
        _mm256_cmp_ps(m, _mm256_set1_ps(-0.5f), _CMP_LT_OQ)));
    const __m256 fc = _mm256_and_ps(sgn, _mm256_or_ps(       // m = 1,2,-2
        _mm256_cmp_ps(m, _mm256_set1_ps(0.5f), _CMP_GT_OQ),
        _mm256_cmp_ps(m, _mm256_set1_ps(-1.5f), _CMP_LT_OQ)));
    const __m256 d = _mm256_and_ps(_mm256_xor_ps(s, c), odd);
    sn = _mm256_xor_ps(_mm256_xor_ps(s, d), fs);
    cs = _mm256_xor_ps(_mm256_xor_ps(c, d), fc);
  }

  TA1 inline void sincosCycAvx(const __m256 ph, __m256 &sn, __m256 &cs)
  {
    // as in the AVX2 version (q*0.25 is exact, so is ph - q/4)
    const __m256 q = _mm256_round_ps(_mm256_mul_ps(ph, _mm256_set1_ps(4.0f)),
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    const __m256 t = _mm256_sub_ps(ph, _mm256_mul_ps(q, _mm256_set1_ps(0.25f)));
    const __m256 z = _mm256_mul_ps(t, t);
    __m256 s = _mm256_add_ps(_mm256_mul_ps(z, _mm256_set1_ps(-75.4456558f)),
        _mm256_set1_ps(81.5937653f));
    s = _mm256_add_ps(_mm256_mul_ps(z, s), _mm256_set1_ps(-41.341671f));
    s = _mm256_add_ps(_mm256_mul_ps(z, s), _mm256_set1_ps(6.28318548f));
    s = _mm256_mul_ps(t, s);
    __m256 c = _mm256_add_ps(_mm256_mul_ps(z, _mm256_set1_ps(59.3497009f)),
        _mm256_set1_ps(-85.4471436f));
    c = _mm256_add_ps(_mm256_mul_ps(z, c), _mm256_set1_ps(64.9393616f));
    c = _mm256_add_ps(_mm256_mul_ps(z, c), _mm256_set1_ps(-19.7392082f));
    c = _mm256_add_ps(_mm256_mul_ps(z, c), _mm256_set1_ps(1.0f));
    quadrantAvx(q, s, c, sn, cs);
  }

  TA1 inline void sincosRadAvx(const __m256 x, __m256 &sn, __m256 &cs)
  {
    const __m256 q = _mm256_round_ps(_mm256_mul_ps(x,
        _mm256_set1_ps(0.636619747f)),
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m256 r = _mm256_sub_ps(x, _mm256_mul_ps(q, _mm256_set1_ps(1.5703125f)));
    r = _mm256_sub_ps(r, _mm256_mul_ps(q, _mm256_set1_ps(0.000482559204f)));
    r = _mm256_sub_ps(r, _mm256_mul_ps(q, _mm256_set1_ps(1.2665987e-06f)));
    r = _mm256_sub_ps(r, _mm256_mul_ps(q, _mm256_set1_ps(9.89530236e-10f)));
    r = _mm256_sub_ps(r, _mm256_mul_ps(q, _mm256_set1_ps(2.56334407e-12f)));
    const __m256 z = _mm256_mul_ps(r, r);
    __m256 s = _mm256_add_ps(_mm256_mul_ps(z, _mm256_set1_ps(-0.000195152956f)),
        _mm256_set1_ps(0.00833216123f));
    s = _mm256_add_ps(_mm256_mul_ps(z, s), _mm256_set1_ps(-0.166666552f));
    s = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(z, r), s), r);
    __m256 c = _mm256_add_ps(_mm256_mul_ps(z, _mm256_set1_ps(2.44331568e-05f)),
        _mm256_set1_ps(-0.00138873165f));
    c = _mm256_add_ps(_mm256_mul_ps(z, c), _mm256_set1_ps(0.0416666456f));
    c = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(z, z), c), _mm256_sub_ps(
        _mm256_set1_ps(1.0f), _mm256_mul_ps(z, _mm256_set1_ps(0.5f))));
    quadrantAvx(q, s, c, sn, cs);
    // lanes past the accurate reduction range (rare) go through sincosf
    const int far = _mm256_movemask_ps(_mm256_cmp_ps(_mm256_andnot_ps(
        _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)), x),
        _mm256_set1_ps(SINCOS_RAD_LIMIT), _CMP_GT_OQ));
    if(far) sincosRadFar256(x, far, sn, cs);
  }

  template<bool CYC>
  TA1 inline __m256 atan2Avx(const __m256 y, const __m256 x)
  {
    const __m256 noSign = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    const __m256 one = _mm256_set1_ps(1.0f);
    const __m256 ax = _mm256_and_ps(x, noSign);
    const __m256 ay = _mm256_and_ps(y, noSign);
    const __m256 sw = _mm256_cmp_ps(ay, ax, _CMP_GT_OQ);
    const __m256 mn = sel256(sw, ay, ax);
    const __m256 mx = sel256(sw, ax, ay);
    __m256 big = _mm256_cmp_ps(mn, _mm256_mul_ps(mx,
        _mm256_set1_ps(0.414213568f)), _CMP_GT_OQ);
    const __m256 inf = _mm256_cmp_ps(mn, _mm256_set1_ps(__builtin_inff()),
        _CMP_EQ_OQ);
    const __m256 zero = _mm256_cmp_ps(mx, _mm256_setzero_ps(), _CMP_EQ_OQ);
    __m256 n = sel256(big, mn, _mm256_sub_ps(mn, mx));
    __m256 d = sel256(big, mx, _mm256_add_ps(mn, mx));
    n = _mm256_andnot_ps(inf, n);
    d = sel256(_mm256_or_ps(inf, zero), d, one);
    big = _mm256_or_ps(big, inf);
    const __m256 u = _mm256_div_ps(n, d);
    const __m256 z = _mm256_mul_ps(u, u);
    // x negative by its sign bit (so -0 counts): +-1 has x's sign
    const __m256 xneg = _mm256_cmp_ps(_mm256_or_ps(_mm256_andnot_ps(noSign, x),
        one), _mm256_setzero_ps(), _CMP_LT_OQ);
    __m256 a;
    if(CYC)
    {
      __m256 p = _mm256_add_ps(_mm256_mul_ps(z, _mm256_set1_ps(0.0128179323f)),
          _mm256_set1_ps(-0.0220870227f));
      p = _mm256_add_ps(_mm256_mul_ps(z, p), _mm256_set1_ps(0.0317955129f));
      p = _mm256_add_ps(_mm256_mul_ps(z, p), _mm256_set1_ps(-0.0530510359f));
      a = _mm256_add_ps(_mm256_mul_ps(u, _mm256_set1_ps(0.159154937f)),
          _mm256_add_ps(_mm256_mul_ps(u, _mm256_set1_ps(6.42063824e-09f)),
          _mm256_mul_ps(_mm256_mul_ps(u, z), p)));
      a = _mm256_add_ps(a, _mm256_and_ps(big, _mm256_set1_ps(0.125f)));
      a = sel256(sw, a, _mm256_sub_ps(_mm256_set1_ps(0.25f), a));
      a = sel256(xneg, a, _mm256_sub_ps(_mm256_set1_ps(0.5f), a));
    }
    else
    {
      __m256 p = _mm256_add_ps(_mm256_mul_ps(z, _mm256_set1_ps(0.0805374458f)),
          _mm256_set1_ps(-0.138776854f));
      p = _mm256_add_ps(_mm256_mul_ps(z, p), _mm256_set1_ps(0.199777111f));
      p = _mm256_add_ps(_mm256_mul_ps(z, p), _mm256_set1_ps(-0.333329499f));
      a = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(u, z), p), u);
      a = _mm256_add_ps(_mm256_add_ps(a, _mm256_and_ps(big,
          _mm256_set1_ps(-2.18556941e-08f))), _mm256_and_ps(big,
          _mm256_set1_ps(0.785398185f)));
      a = sel256(sw, a, _mm256_add_ps(_mm256_sub_ps(
          _mm256_set1_ps(1.57079637f), a), _mm256_set1_ps(-4.37113883e-08f)));
      a = sel256(xneg, a, _mm256_add_ps(_mm256_sub_ps(
          _mm256_set1_ps(3.14159274f), a), _mm256_set1_ps(-8.74227766e-08f)));
    }
    a = _mm256_or_ps(a, _mm256_andnot_ps(noSign, y));
    const __m256 nan = _mm256_cmp_ps(x, y, _CMP_UNORD_Q);
    return sel256(nan, a, _mm256_add_ps(x, y));
  }

  //=============================================================== atan2 ===//
  // atan2(y, x) in cycles (Cyc) or radians (Rad)
  TA5 inline __m512 atan2Cyc512(const __m512 y, const __m512 x)
  { return atan2_512<true>(y, x); }
  TA5 inline __m512 atan2Rad512(const __m512 y, const __m512 x)
  { return atan2_512<false>(y, x); }
  TA2 inline __m256 atan2CycAvx2(const __m256 y, const __m256 x)
  { return atan2Avx2<true>(y, x); }
  TA2 inline __m256 atan2RadAvx2(const __m256 y, const __m256 x)
  { return atan2Avx2<false>(y, x); }
  TA1 inline __m256 atan2CycAvx(const __m256 y, const __m256 x)
  { return atan2Avx<true>(y, x); }
  TA1 inline __m256 atan2RadAvx(const __m256 y, const __m256 x)
  { return atan2Avx<false>(y, x); }

  #undef TA5
  #undef TA2
  #undef TA1
} // namespace TRIG
} // namespace EVM
#endif // EVM_TRIG_H
