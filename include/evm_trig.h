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
// Purpose:  SIMD sine, cosine, atan2, exp and log for evm.h that don't need
//           SVML
//
// Created:  2026/11/05
//
// Description:
//  Inline helpers that compute the sine and cosine, the atan2, the exp or the
//  log of a whole register at once, for the AVX512 (avx512f), AVX2
//  (avx2,fma) and AVX (avx) versions of the evm.h functions.  Each helper has
//  the target attribute of the instruction set it is written for, so it can
//  only be called from a function with the same (or a larger) target:
//    sincosCyc  sine and cosine of 2*pi*ph, phase in cycles
//    sincosCycDbl  the same for phases given as doubles (two registers),
//               reduced in double first (for the tune kernels)
//    sincosRad  sine and cosine of x, x in radians
//    atan2Cyc   atan2(y, x) in cycles
//    atan2Rad   atan2(y, x) in radians
//    exp        e^x
//    log        natural logarithm of x
//  with the suffix 512, Avx2 or Avx.
//
//  The polynomials are the Cephes single precision minimax ones (Stephen
//  Moshier), the cycles versions fold 2*pi (sin/cos) or 1/(2*pi) (atan) into
//  the coefficients.  exp and log are Cephes only on AVX: on AVX512 and AVX2
//  they are table-driven (8 to 32 entry tables read with permutes and
//  shorter near-minimax polynomials), 1.6-2x faster and somewhat less
//  accurate.  Accuracy (maximum error measured over random and edge case
//  inputs, for exp and log over every float input):
//    sincosCyc  1.7 ulp at any phase: whole cycles and quarter cycles are
//               removed exactly, unlike sincos(2*pi*ph), which loses accuracy
//               as the phase grows (abs error 8e-7 by 2 cycles, 4e-3 by 1e4);
//               the AVX2/AVX versions give NaN past |ph| = 8.5e37
//    sincosRad  1.6 ulp (AVX512/AVX2) or 3.3 ulp (AVX) for |x| <= 65536,
//               larger |x| fall back to scalar sincosf for those lanes
//    atan2      2.5 ulp (radians), 2.3 ulp (cycles)
//    exp        0.99 ulp (AVX512, AVX), 1.14 ulp (AVX2); past |x| = 87.3 x
//               is clamped to [-104, 89], where e^x rounds to 0 or INF anyway
//    log        2.4 ulp (AVX512), 2.8 ulp (AVX2), 0.83 ulp (AVX)
//  Zeros (with their signs), infinities and NaN give the same results as the
//  C library's atan2f, sinf, cosf, expf and logf, and so do exp's denormal
//  results and log's denormal inputs (with FTZ/DAZ on, a denormal result is
//  0 and a denormal input acts as 0, as in the C library).
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

  // sin(2*pi*t) and cos(2*pi*t) for |t| <= 1/8
  TA5 inline void sincosCycPoly512(const __m512 t, __m512 &s, __m512 &c)
  {
    const __m512 z = _mm512_mul_ps(t, t);
    s = _mm512_fmadd_ps(z, _mm512_set1_ps(-75.4456558f),
        _mm512_set1_ps(81.5937653f));
    s = _mm512_fmadd_ps(z, s, _mm512_set1_ps(-41.341671f));
    s = _mm512_fmadd_ps(z, s, _mm512_set1_ps(6.28318548f));
    s = _mm512_mul_ps(t, s);                         // sin(2*pi*t)
    c = _mm512_fmadd_ps(z, _mm512_set1_ps(59.3497009f),
        _mm512_set1_ps(-85.4471436f));
    c = _mm512_fmadd_ps(z, c, _mm512_set1_ps(64.9393616f));
    c = _mm512_fmadd_ps(z, c, _mm512_set1_ps(-19.7392082f));
    c = _mm512_fmadd_ps(z, c, _mm512_set1_ps(1.0f)); // cos(2*pi*t)
  }

  TA5 inline void sincosCyc512(const __m512 ph, __m512 &sn, __m512 &cs)
  {
    // q4 = ph rounded to the nearest quarter cycle (roundscale keeping 2
    // fraction bits, so nothing can overflow), ph - q4 is exact (|t| <= 1/8),
    // so the reduction adds no error however large ph is
    const __m512 q4 = _mm512_roundscale_ps(ph, (2<<4) |
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m512 s, c;
    sincosCycPoly512(_mm512_sub_ps(ph, q4), s, c);
    // quadrant from 4*q4 (for |ph| >= 2^29 the conversion overflows to
    // 0x80000000, quadrant 0, which is right for a whole number of cycles)
    quadrant512(_mm512_cvtps_epi32(_mm512_mul_ps(q4, _mm512_set1_ps(4.0f))),
        s, c, sn, cs);
  }

  // sincosCyc of 16 phases given as doubles, lanes 0-7 from lo and 8-15 from
  // hi.  The nearest quarter cycle is removed in double, so the float
  // argument of the polynomials keeps the accuracy of the double phase (for
  // |ph| < 2^49): y = 4*ph + 1.5*2^52 rounds 4*ph to an integer q held in the
  // low bits of y (the quadrant), q/4 = y/4 - 1.5*2^50 exactly, and
  // t = ph - q/4 is exact (|t| <= 1/8)
  TA5 inline void sincosCycDbl512(const __m512d lo, const __m512d hi,
      __m512 &sn, __m512 &cs)
  {
    const __m512d four = _mm512_set1_pd(4.0);
    const __m512d quarter = _mm512_set1_pd(0.25);
    const __m512d magic = _mm512_set1_pd(6755399441055744.0);  // 1.5*2^52
    const __m512d magic4 = _mm512_set1_pd(1688849860263936.0); // 1.5*2^50
    const __m512d ylo = _mm512_fmadd_pd(lo, four, magic);
    const __m512d yhi = _mm512_fmadd_pd(hi, four, magic);
    const __m256 tlo = _mm512_cvtpd_ps(_mm512_sub_pd(lo,
        _mm512_fmsub_pd(ylo, quarter, magic4)));
    const __m256 thi = _mm512_cvtpd_ps(_mm512_sub_pd(hi,
        _mm512_fmsub_pd(yhi, quarter, magic4)));
    const __m512 t = _mm512_castpd_ps(_mm512_insertf64x4(
        _mm512_castpd256_pd512(_mm256_castps_pd(tlo)), _mm256_castps_pd(thi),
        1));
    // the low 32 bits of each y, in the same lane order as t
    const __m512i q = _mm512_permutex2var_epi32(_mm512_castpd_si512(ylo),
        _mm512_setr_epi32(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26,
        28, 30), _mm512_castpd_si512(yhi));
    __m512 s, c;
    sincosCycPoly512(t, s, c);
    quadrant512(q, s, c, sn, cs);
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
    // both infinite: 0/mx = 0/INF, so the angle is the pi/4 base; both zero:
    // 0/1, so the angle is 0.  NaN needs no special case: one of mn and mx
    // is the NaN and reaches u through n or d (a NaN x with an infinite y
    // gives 0/NaN, since d stays mx there rather than becoming 1)
    const __mmask16 inf = _mm512_cmp_ps_mask(mn,
        _mm512_set1_ps(__builtin_inff()), _CMP_EQ_OQ);
    const __mmask16 zero = _mm512_cmp_ps_mask(mx, _mm512_setzero_ps(),
        _CMP_EQ_OQ);
    __m512 n = _mm512_mask_sub_ps(mn, big, mn, mx);
    __m512 d = _mm512_mask_add_ps(mx, big, mn, mx);
    n = _mm512_mask_mov_ps(n, inf, _mm512_setzero_ps());
    d = _mm512_mask_mov_ps(d, zero, _mm512_set1_ps(1.0f));
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
    return a;
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

  // sin(2*pi*t) and cos(2*pi*t) for |t| <= 1/8
  TA2 inline void sincosCycPolyAvx2(const __m256 t, __m256 &s, __m256 &c)
  {
    const __m256 z = _mm256_mul_ps(t, t);
    s = _mm256_fmadd_ps(z, _mm256_set1_ps(-75.4456558f),
        _mm256_set1_ps(81.5937653f));
    s = _mm256_fmadd_ps(z, s, _mm256_set1_ps(-41.341671f));
    s = _mm256_fmadd_ps(z, s, _mm256_set1_ps(6.28318548f));
    s = _mm256_mul_ps(t, s);
    c = _mm256_fmadd_ps(z, _mm256_set1_ps(59.3497009f),
        _mm256_set1_ps(-85.4471436f));
    c = _mm256_fmadd_ps(z, c, _mm256_set1_ps(64.9393616f));
    c = _mm256_fmadd_ps(z, c, _mm256_set1_ps(-19.7392082f));
    c = _mm256_fmadd_ps(z, c, _mm256_set1_ps(1.0f));
  }

  TA2 inline void sincosCycAvx2(const __m256 ph, __m256 &sn, __m256 &cs)
  {
    // q = nearest quarter cycle, ph - q/4 is exact (|t| <= 1/8), so the
    // reduction adds no error however large ph is, up to |ph| = 8.5e37 where
    // 4*ph overflows and the result is NaN (removing whole cycles first would
    // avoid that, but costs 13-20%)
    const __m256 q = _mm256_round_ps(_mm256_mul_ps(ph, _mm256_set1_ps(4.0f)),
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m256 s, c;
    sincosCycPolyAvx2(_mm256_fnmadd_ps(q, _mm256_set1_ps(0.25f), ph), s, c);
    quadrantAvx2(_mm256_cvtps_epi32(q), s, c, sn, cs);
  }

  // sincosCyc of 8 phases given as doubles, lanes 0-3 from lo and 4-7 from
  // hi, reduced in double as in sincosCycDbl512
  TA2 inline void sincosCycDblAvx2(const __m256d lo, const __m256d hi,
      __m256 &sn, __m256 &cs)
  {
    const __m256d four = _mm256_set1_pd(4.0);
    const __m256d quarter = _mm256_set1_pd(0.25);
    const __m256d magic = _mm256_set1_pd(6755399441055744.0);  // 1.5*2^52
    const __m256d magic4 = _mm256_set1_pd(1688849860263936.0); // 1.5*2^50
    const __m256d ylo = _mm256_fmadd_pd(lo, four, magic);
    const __m256d yhi = _mm256_fmadd_pd(hi, four, magic);
    const __m128 tlo = _mm256_cvtpd_ps(_mm256_sub_pd(lo,
        _mm256_fmsub_pd(ylo, quarter, magic4)));
    const __m128 thi = _mm256_cvtpd_ps(_mm256_sub_pd(hi,
        _mm256_fmsub_pd(yhi, quarter, magic4)));
    const __m256 t = _mm256_insertf128_ps(_mm256_castps128_ps256(tlo), thi, 1);
    // the low 32 bits of each y: shuffle_ps gives the lanes in the order
    // 0,1,4,5,2,3,6,7, which permute4x64 puts back in t's order
    const __m256i q = _mm256_permute4x64_epi64(_mm256_castps_si256(
        _mm256_shuffle_ps(_mm256_castpd_ps(ylo), _mm256_castpd_ps(yhi),
        0x88)), 0xd8);
    __m256 s, c;
    sincosCycPolyAvx2(t, s, c);
    quadrantAvx2(q, s, c, sn, cs);
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
    d = sel256(zero, d, _mm256_set1_ps(1.0f)); // see atan2_512
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
    return a;
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

  // sin(2*pi*t) and cos(2*pi*t) for |t| <= 1/8
  TA1 inline void sincosCycPolyAvx(const __m256 t, __m256 &s, __m256 &c)
  {
    const __m256 z = _mm256_mul_ps(t, t);
    s = _mm256_add_ps(_mm256_mul_ps(z, _mm256_set1_ps(-75.4456558f)),
        _mm256_set1_ps(81.5937653f));
    s = _mm256_add_ps(_mm256_mul_ps(z, s), _mm256_set1_ps(-41.341671f));
    s = _mm256_add_ps(_mm256_mul_ps(z, s), _mm256_set1_ps(6.28318548f));
    s = _mm256_mul_ps(t, s);
    c = _mm256_add_ps(_mm256_mul_ps(z, _mm256_set1_ps(59.3497009f)),
        _mm256_set1_ps(-85.4471436f));
    c = _mm256_add_ps(_mm256_mul_ps(z, c), _mm256_set1_ps(64.9393616f));
    c = _mm256_add_ps(_mm256_mul_ps(z, c), _mm256_set1_ps(-19.7392082f));
    c = _mm256_add_ps(_mm256_mul_ps(z, c), _mm256_set1_ps(1.0f));
  }

  TA1 inline void sincosCycAvx(const __m256 ph, __m256 &sn, __m256 &cs)
  {
    // as in the AVX2 version (q*0.25 is exact, so is ph - q/4)
    const __m256 q = _mm256_round_ps(_mm256_mul_ps(ph, _mm256_set1_ps(4.0f)),
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m256 s, c;
    sincosCycPolyAvx(_mm256_sub_ps(ph, _mm256_mul_ps(q,
        _mm256_set1_ps(0.25f))), s, c);
    quadrantAvx(q, s, c, sn, cs);
  }

  // sincosCyc of 8 phases given as doubles, lanes 0-3 from lo and 4-7 from
  // hi, reduced in double: f = ph - round(ph) and t = f - q/4 with
  // q = round(4*f) are exact, and q (-2 to 2) converts to float exactly
  TA1 inline void sincosCycDblAvx(const __m256d lo, const __m256d hi,
      __m256 &sn, __m256 &cs)
  {
    const int rnd = _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC;
    const __m256d four = _mm256_set1_pd(4.0);
    const __m256d quarter = _mm256_set1_pd(0.25);
    const __m256d flo = _mm256_sub_pd(lo, _mm256_round_pd(lo, rnd));
    const __m256d fhi = _mm256_sub_pd(hi, _mm256_round_pd(hi, rnd));
    const __m256d qlo = _mm256_round_pd(_mm256_mul_pd(flo, four), rnd);
    const __m256d qhi = _mm256_round_pd(_mm256_mul_pd(fhi, four), rnd);
    const __m256 t = _mm256_insertf128_ps(_mm256_castps128_ps256(
        _mm256_cvtpd_ps(_mm256_sub_pd(flo, _mm256_mul_pd(qlo, quarter)))),
        _mm256_cvtpd_ps(_mm256_sub_pd(fhi, _mm256_mul_pd(qhi, quarter))), 1);
    const __m256 q = _mm256_insertf128_ps(_mm256_castps128_ps256(
        _mm256_cvtpd_ps(qlo)), _mm256_cvtpd_ps(qhi), 1);
    __m256 s, c;
    sincosCycPolyAvx(t, s, c);
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
    d = sel256(zero, d, one); // see atan2_512
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
    return a;
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

  //======================================================== exp and log ===//
  // AVX512 and AVX2 are table-driven:
  //   e^x = 2^(n/N) * e^r: t = x*N/ln2 + 1.5*2^23 rounds n into the low bits
  //   of t (the table index), n/N = t/N - 1.5*2^23/N exactly, r = x -
  //   (n/N)*ln2 with ln2 as the float nearest it plus the rest (|r| <=
  //   ln2/(2N)) and e^r - 1 = r + r^2*(P0 + P1*r).  The table holds 2^(j/N).
  //   log(x) = k*ln2 + log(c) + log1p(r) with x = 2^k*m, c the center of the
  //   table entry m falls in (1 for the entries next to 1.0, so that results
  //   near 0 keep their accuracy), r = m/c - 1 and log1p(r) = r +
  //   r^2*P(r).  The tables hold 1/c and log(c) = -log(1/c), so r and log(c)
  //   agree whatever the rounding of 1/c.
  // The polynomials are near-minimax for their r ranges.  AVX uses the Cephes
  // expf and logf polynomials: e^x = 2^n * e^r with r = x - n*ln2 (ln2 as
  // Cephes' 9-bit 0.693359375, which n times is exact, plus the rest),
  // e^r = 1 + r + r^2*P(r), and log(x) = e*ln2 + log(1+f) with 1+f in
  // [sqrt(.5), sqrt(2)), log(1+f) = f - f^2/2 + f^3*P(f).
  //
  // Registers with a lane past |x| = 87.3 (or NaN), where e^x is not a
  // normal float, take a branch to a version that clamps x to [-104, 89]
  // (x second, so a NaN passes through min/max), past which e^x rounds to 0
  // or INF anyway, and applies 2^n so that a denormal result or an overflow
  // rounds right.  Clamping every register cost 10% (AVX512) to 25% (AVX2).
  constexpr float EXP_FAST = 87.3f, EXP_LO = -104.0f, EXP_HI = 89.0f;
  // 2^(j/16) and 2^(j/8), e^r - 1 coefficients for |r| <= ln2/32 and ln2/16
  alignas(64) inline constexpr float EXP_T16[16] = {
    1.000000000e+00f, 1.044273734e+00f, 1.090507746e+00f, 1.138788581e+00f,
    1.189207077e+00f, 1.241857767e+00f, 1.296839595e+00f, 1.354255557e+00f,
    1.414213538e+00f, 1.476826191e+00f, 1.542210817e+00f, 1.610490322e+00f,
    1.681792855e+00f, 1.756252170e+00f, 1.834008098e+00f, 1.915206552e+00f};
  alignas(32) inline constexpr float EXP_T8[8] = {
    1.000000000e+00f, 1.090507746e+00f, 1.189207077e+00f, 1.296839595e+00f,
    1.414213538e+00f, 1.542210817e+00f, 1.681792855e+00f, 1.834008098e+00f};
  constexpr float EXP16_P0 = 5.000162125e-01f, EXP16_P1 = 1.666704565e-01f;
  constexpr float EXP8_P0 = 5.000647306e-01f, EXP8_P1 = 1.666822731e-01f;
  // log, AVX512: m = getmant(x) in [0.75, 1.5), entry j = top 5 mantissa
  // bits of m: [1+j/32, 1+(j+1)/32) for j < 16, [0.5+j/64, 0.5+(j+1)/64)
  // for j >= 16, |r| <= 1/32
  alignas(64) inline constexpr float LOG_INVC32[32] = {
    1.000000000e+00f, 9.552238584e-01f, 9.275362492e-01f, 9.014084339e-01f,
    8.767123222e-01f, 8.533333540e-01f, 8.311688304e-01f, 8.101266026e-01f,
    7.901234627e-01f, 7.710843086e-01f, 7.529411912e-01f, 7.356321812e-01f,
    7.191011310e-01f, 7.032967210e-01f, 6.881720424e-01f, 6.736842394e-01f,
    1.319587588e+00f, 1.292929292e+00f, 1.267326713e+00f, 1.242718458e+00f,
    1.219047666e+00f, 1.196261644e+00f, 1.174311876e+00f, 1.153153181e+00f,
    1.132743359e+00f, 1.113043427e+00f, 1.094017148e+00f, 1.075630307e+00f,
    1.057851195e+00f, 1.040650368e+00f, 1.024000049e+00f, 1.000000000e+00f};
  alignas(64) inline constexpr float LOG_LOGC32[32] = {
    0.000000000e+00f, 4.580955952e-02f, 7.522340119e-02f, 1.037968099e-01f,
    1.315763593e-01f, 1.586050093e-01f, 1.849223375e-01f, 2.105647475e-01f,
    2.355660647e-01f, 2.599575520e-01f, 2.837681472e-01f, 3.070250452e-01f,
    3.297532797e-01f, 3.519763947e-01f, 3.737164140e-01f, 3.949937522e-01f,
    -2.773192525e-01f, -2.569104135e-01f, -2.369097322e-01f, -2.173012793e-01f,
    -1.980699450e-01f, -1.792013943e-01f, -1.606823355e-01f, -1.425000876e-01f,
    -1.246424392e-01f, -1.070980877e-01f, -8.985637873e-02f, -7.290682197e-02f,
    -5.623967573e-02f, -3.984587267e-02f, -2.371657453e-02f, 0.000000000e+00f};
  constexpr float LOG32_P0 = -5.000011921e-01f, LOG32_P1 = 3.333663344e-01f,
    LOG32_P2 = -2.439804077e-01f;
  // log, AVX2: z = x/2^k in [0.71875, 1.4375) and entry j from the bits of x
  // minus those of 0.71875, 2^20 apart, so that 1.0 is in the middle of
  // entry 4 (in bits), |r| <= 1/16
  constexpr int LOG8_OFF = 0x3f380000;
  alignas(32) inline constexpr float LOG_INVC8[8] = {
    1.333333373e+00f, 1.230769277e+00f, 1.142857194e+00f, 1.066666722e+00f,
    1.000000000e+00f, 8.888888955e-01f, 8.000000119e-01f, 7.272727489e-01f};
  alignas(32) inline constexpr float LOG_LOGC8[8] = {
    -2.876821160e-01f, -2.076393962e-01f, -1.335314363e-01f, -6.453857571e-02f,
    0.000000000e+00f, 1.177830249e-01f, 2.231435329e-01f, 3.184536994e-01f};
  constexpr float LOG8_P0 = -4.999994338e-01f, LOG8_P1 = 3.333380520e-01f,
    LOG8_P2 = -2.506930828e-01f, LOG8_P3 = 1.981910318e-01f;
  // ln2 as the float nearest it plus the rest (FMA), and Cephes' (AVX)
  constexpr float LN2_HI = 0.693147182f, LN2_LO = -1.90465421e-09f;
  constexpr float LN2_C1 = 0.693359375f, LN2_C2 = -2.12194442e-4f;
  // Cephes expf and logf, AVX
  constexpr float EXP_P5 = 1.9875691500e-4f, EXP_P4 = 1.3981999507e-3f,
    EXP_P3 = 8.3334519073e-3f, EXP_P2 = 4.1665795894e-2f,
    EXP_P1 = 1.6666665459e-1f, EXP_P0 = 5.0000001201e-1f;
  constexpr float LOG_P8 = 7.0376836292e-2f, LOG_P7 = -1.1514610310e-1f,
    LOG_P6 = 1.1676998740e-1f, LOG_P5 = -1.2420140846e-1f,
    LOG_P4 = 1.4249322787e-1f, LOG_P3 = -1.6668057665e-1f,
    LOG_P2 = 2.0000714765e-1f, LOG_P1 = -2.4999993993e-1f,
    LOG_P0 = 3.3333331174e-1f;

  // 2^(n/16) * e^r of an x already in range (see above)
  TA5 inline __m512 expTable512(const __m512 x)
  {
    const __m512 magic = _mm512_set1_ps(12582912.0f);   // 1.5*2^23
    const __m512 t = _mm512_fmadd_ps(x, _mm512_set1_ps(2.308312035e+01f),
        magic);                                          // n + magic
    const __m512 n16 = _mm512_fmsub_ps(t, _mm512_set1_ps(0.0625f),
        _mm512_set1_ps(786432.0f));                      // n/16
    __m512 r = _mm512_fnmadd_ps(n16, _mm512_set1_ps(LN2_HI), x);
    r = _mm512_fnmadd_ps(n16, _mm512_set1_ps(LN2_LO), r);
    const __m512 tb = _mm512_permutexvar_ps(_mm512_castps_si512(t),
        _mm512_load_ps(EXP_T16));                        // 2^((n mod 16)/16)
    const __m512 em1 = _mm512_fmadd_ps(_mm512_fmadd_ps(_mm512_set1_ps(
        EXP16_P1), r, _mm512_set1_ps(EXP16_P0)), _mm512_mul_ps(r, r), r);
    // times 2^floor(n/16), rounded once, so a denormal result is right and
    // an overflow INF
    return _mm512_scalef_ps(_mm512_fmadd_ps(tb, em1, tb), n16);
  }

  // exp for registers with a lane past |x| = 87.3 or a NaN
  __attribute__((__target__("avx512f"), __noinline__))
  inline __m512 expFar512(const __m512 x)
  {
    return expTable512(_mm512_max_ps(_mm512_set1_ps(EXP_LO),
        _mm512_min_ps(_mm512_set1_ps(EXP_HI), x)));
  }

  TA5 inline __m512 exp512(const __m512 x)
  {
    if(_mm512_cmp_ps_mask(_mm512_abs_ps(x), _mm512_set1_ps(EXP_FAST),
        _CMP_NLE_UQ)) return expFar512(x);
    return expTable512(x);
  }

  TA5 inline __m512 log512(const __m512 x)
  {
    // x = 2^k*m with m in [0.75, 1.5) (getexp/getmant also handle
    // denormals): getmant halves a mantissa of 1.5 or more, and k is then
    // one more than getexp
    const __m512 one = _mm512_set1_ps(1.0f);
    const __m512 m = _mm512_getmant_ps(x, _MM_MANT_NORM_p75_1p5,
        _MM_MANT_SIGN_zero);
    const __m512 ex = _mm512_getexp_ps(x);
    const __m512 k = _mm512_mask_add_ps(ex, _mm512_cmp_ps_mask(m, one,
        _CMP_LT_OQ), ex, one);
    // entry from the top 5 mantissa bits of m (permutex2var uses the low 5
    // bits of the index)
    const __m512i j = _mm512_srli_epi32(_mm512_castps_si512(m), 18);
    const __m512 ic = _mm512_permutex2var_ps(_mm512_load_ps(LOG_INVC32), j,
        _mm512_load_ps(LOG_INVC32 + 16));
    const __m512 lc = _mm512_permutex2var_ps(_mm512_load_ps(LOG_LOGC32), j,
        _mm512_load_ps(LOG_LOGC32 + 16));
    const __m512 r = _mm512_fmsub_ps(m, ic, one);
    __m512 h = _mm512_fmadd_ps(_mm512_set1_ps(LOG32_P2), r,
        _mm512_set1_ps(LOG32_P1));
    h = _mm512_fmadd_ps(h, r, _mm512_set1_ps(LOG32_P0));
    const __m512 y = _mm512_fmadd_ps(_mm512_mul_ps(r, r), h, _mm512_add_ps(
        _mm512_fmadd_ps(k, _mm512_set1_ps(LN2_HI), lc), r));
    // special inputs by class, in one instruction (a denormal counts as zero
    // under DAZ): +-0 gives -INF, 1 gives +0, negative numbers and -INF give
    // NaN, +INF gives +INF and a NaN stays NaN
    return _mm512_fixupimm_ps(y, x, _mm512_set1_epi32(0x03538422), 0);
  }

  // 2^(n/8) * e^r of an x in range: the table entry and n
  TA2 inline __m256 expTableAvx2(const __m256 x, __m256i &n)
  {
    const __m256 magic = _mm256_set1_ps(12582912.0f);
    const __m256 t = _mm256_fmadd_ps(x, _mm256_set1_ps(1.154156017e+01f),
        magic);
    const __m256 n8 = _mm256_fmsub_ps(t, _mm256_set1_ps(0.125f),
        _mm256_set1_ps(1572864.0f));
    __m256 r = _mm256_fnmadd_ps(n8, _mm256_set1_ps(LN2_HI), x);
    r = _mm256_fnmadd_ps(n8, _mm256_set1_ps(LN2_LO), r);
    const __m256 tb = _mm256_permutevar8x32_ps(_mm256_load_ps(EXP_T8),
        _mm256_castps_si256(t));
    const __m256 em1 = _mm256_fmadd_ps(_mm256_fmadd_ps(_mm256_set1_ps(
        EXP8_P1), r, _mm256_set1_ps(EXP8_P0)), _mm256_mul_ps(r, r), r);
    n = _mm256_sub_epi32(_mm256_castps_si256(t), _mm256_castps_si256(magic));
    return _mm256_fmadd_ps(tb, em1, tb);
  }

  // 2^n from the bits (n+127)<<23, for n in -126..127
  TA2 inline __m256 pow2Avx2(const __m256i n)
  {
    return _mm256_castsi256_ps(_mm256_slli_epi32(_mm256_add_epi32(n,
        _mm256_set1_epi32(127)), 23));
  }

  // exp for registers with a lane past |x| = 87.3 or a NaN: 2^floor(n/8)
  // (-150..128) as two factors, which round once, in the last multiply
  __attribute__((__target__("avx2,fma"), __noinline__))
  inline __m256 expFarAvx2(const __m256 x)
  {
    __m256i n;
    const __m256 m = expTableAvx2(_mm256_max_ps(_mm256_set1_ps(EXP_LO),
        _mm256_min_ps(_mm256_set1_ps(EXP_HI), x)), n);
    const __m256i k = _mm256_srai_epi32(n, 3);
    const __m256i k1 = _mm256_srai_epi32(k, 1);
    return _mm256_mul_ps(_mm256_mul_ps(m, pow2Avx2(k1)),
        pow2Avx2(_mm256_sub_epi32(k, k1)));
  }

  TA2 inline __m256 expAvx2(const __m256 x)
  {
    if(_mm256_movemask_ps(_mm256_cmp_ps(_mm256_andnot_ps(
        _mm256_set1_ps(-0.0f), x), _mm256_set1_ps(EXP_FAST), _CMP_NLE_UQ)))
      return expFarAvx2(x);
    // every result is a normal float here: 2^floor(n/8) is added to the
    // exponent bits
    __m256i n;
    const __m256 m = expTableAvx2(x, n);
    return _mm256_castsi256_ps(_mm256_add_epi32(_mm256_castps_si256(m),
        _mm256_slli_epi32(_mm256_srai_epi32(n, 3), 23)));
  }

  // r, k and log(c) of a positive normal x: subtracting the bits of the
  // table's start puts k in the top bits and the entry in the next 3, and
  // removing k*2^23 from the bits of x gives z = x/2^k
  TA2 inline __m256 logReduceAvx2(const __m256 x, __m256 &k, __m256 &lc)
  {
    const __m256i ix = _mm256_castps_si256(x);
    const __m256i tmp = _mm256_sub_epi32(ix, _mm256_set1_epi32(LOG8_OFF));
    const __m256i j = _mm256_srli_epi32(tmp, 20); // vpermps uses 3 bits
    k = _mm256_cvtepi32_ps(_mm256_srai_epi32(tmp, 23));
    const __m256 z = _mm256_castsi256_ps(_mm256_sub_epi32(ix,
        _mm256_and_si256(tmp, _mm256_set1_epi32(-0x00800000))));
    lc = _mm256_permutevar8x32_ps(_mm256_load_ps(LOG_LOGC8), j);
    return _mm256_fmsub_ps(z, _mm256_permutevar8x32_ps(_mm256_load_ps(
        LOG_INVC8), j), _mm256_set1_ps(1.0f));
  }

  // k*ln2 + log(c) + log1p(r)
  TA2 inline __m256 logPolyAvx2(const __m256 r, const __m256 k,
      const __m256 lc)
  {
    __m256 h = _mm256_fmadd_ps(_mm256_set1_ps(LOG8_P3), r,
        _mm256_set1_ps(LOG8_P2));
    h = _mm256_fmadd_ps(h, r, _mm256_set1_ps(LOG8_P1));
    h = _mm256_fmadd_ps(h, r, _mm256_set1_ps(LOG8_P0));
    return _mm256_fmadd_ps(_mm256_mul_ps(r, r), h, _mm256_add_ps(
        _mm256_fmadd_ps(k, _mm256_set1_ps(LN2_HI), lc), r));
  }

  // log for registers with an input that isn't a positive normal float:
  // denormals are scaled by 2^23 first (under DAZ they read as zero), +INF
  // gives +INF, +-0 gives -INF, negative numbers and NaN give NaN
  __attribute__((__target__("avx2,fma"), __noinline__))
  inline __m256 logSpecialAvx2(const __m256 x)
  {
    const __m256 zero = _mm256_setzero_ps();
    const __m256 den = _mm256_cmp_ps(x, _mm256_set1_ps(1.17549435e-38f),
        _CMP_LT_OQ);
    __m256 k, lc;
    const __m256 r = logReduceAvx2(_mm256_mul_ps(x, sel256(den,
        _mm256_set1_ps(1.0f), _mm256_set1_ps(8388608.0f))), k, lc);
    __m256 y = logPolyAvx2(r, _mm256_sub_ps(k, _mm256_and_ps(den,
        _mm256_set1_ps(23.0f))), lc);
    y = sel256(_mm256_cmp_ps(x, _mm256_set1_ps(__builtin_inff()),
        _CMP_EQ_OQ), y, x);
    y = sel256(_mm256_cmp_ps(x, zero, _CMP_EQ_OQ), y,
        _mm256_set1_ps(-__builtin_inff()));
    return _mm256_or_ps(y, _mm256_cmp_ps(x, zero, _CMP_NGE_UQ)); // all ones
  }

  TA2 inline __m256 logAvx2(const __m256 x)
  {
    // x not a positive normal float: its bits minus 0x00800000 are not below
    // 0x7f000000 (unsigned), i.e. bits + 0x7f800000 > 0xfeffffff (signed)
    if(_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(
        _mm256_add_epi32(_mm256_castps_si256(x),
        _mm256_set1_epi32(0x7f800000)), _mm256_set1_epi32(-0x01000001)))))
      return logSpecialAvx2(x);
    __m256 k, lc;
    const __m256 r = logReduceAvx2(x, k, lc);
    return logPolyAvx2(r, k, lc);
  }

  // e^r polynomial of x and n (Cephes, no FMA)
  TA1 inline __m256 expPolyAvx(const __m256 x, __m256 &n)
  {
    n = _mm256_round_ps(_mm256_mul_ps(x, _mm256_set1_ps(1.44269502f)),
        _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m256 r = _mm256_sub_ps(x, _mm256_mul_ps(n, _mm256_set1_ps(LN2_C1)));
    r = _mm256_sub_ps(r, _mm256_mul_ps(n, _mm256_set1_ps(LN2_C2)));
    __m256 p = _mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(EXP_P5), r),
        _mm256_set1_ps(EXP_P4));
    p = _mm256_add_ps(_mm256_mul_ps(p, r), _mm256_set1_ps(EXP_P3));
    p = _mm256_add_ps(_mm256_mul_ps(p, r), _mm256_set1_ps(EXP_P2));
    p = _mm256_add_ps(_mm256_mul_ps(p, r), _mm256_set1_ps(EXP_P1));
    p = _mm256_add_ps(_mm256_mul_ps(p, r), _mm256_set1_ps(EXP_P0));
    return _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(p, _mm256_mul_ps(r, r)),
        r), _mm256_set1_ps(1.0f));
  }

  // 2^n for an integer valued n in -126..127: converting (n+127)*2^23
  // (exact) to an integer gives the bits (n+127)<<23
  TA1 inline __m256 pow2Avx(const __m256 n)
  {
    return _mm256_castsi256_ps(_mm256_cvtps_epi32(_mm256_mul_ps(
        _mm256_add_ps(n, _mm256_set1_ps(127.0f)),
        _mm256_set1_ps(8388608.0f))));
  }

  // as expFarAvx2: n is -150..129, 2^n as two factors
  __attribute__((__target__("avx"), __noinline__))
  inline __m256 expFarAvx(const __m256 x)
  {
    __m256 n;
    const __m256 p = expPolyAvx(_mm256_max_ps(_mm256_set1_ps(EXP_LO),
        _mm256_min_ps(_mm256_set1_ps(EXP_HI), x)), n);
    const __m256 n1 = _mm256_round_ps(_mm256_mul_ps(n, _mm256_set1_ps(0.5f)),
        _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
    return _mm256_mul_ps(_mm256_mul_ps(p, pow2Avx(n1)),
        pow2Avx(_mm256_sub_ps(n, n1)));
  }

  TA1 inline __m256 expAvx(const __m256 x)
  {
    if(_mm256_movemask_ps(_mm256_cmp_ps(_mm256_andnot_ps(
        _mm256_set1_ps(-0.0f), x), _mm256_set1_ps(EXP_FAST), _CMP_NLE_UQ)))
      return expFarAvx(x);
    __m256 n;
    const __m256 p = expPolyAvx(x, n);
    return _mm256_mul_ps(p, pow2Avx(n));
  }

  // log(1+f) + e*ln2 (Cephes, no FMA)
  TA1 inline __m256 logPolyAvx(const __m256 f, const __m256 e)
  {
    const __m256 z = _mm256_mul_ps(f, f);
    __m256 y = _mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(LOG_P8), f),
        _mm256_set1_ps(LOG_P7));
    y = _mm256_add_ps(_mm256_mul_ps(y, f), _mm256_set1_ps(LOG_P6));
    y = _mm256_add_ps(_mm256_mul_ps(y, f), _mm256_set1_ps(LOG_P5));
    y = _mm256_add_ps(_mm256_mul_ps(y, f), _mm256_set1_ps(LOG_P4));
    y = _mm256_add_ps(_mm256_mul_ps(y, f), _mm256_set1_ps(LOG_P3));
    y = _mm256_add_ps(_mm256_mul_ps(y, f), _mm256_set1_ps(LOG_P2));
    y = _mm256_add_ps(_mm256_mul_ps(y, f), _mm256_set1_ps(LOG_P1));
    y = _mm256_add_ps(_mm256_mul_ps(y, f), _mm256_set1_ps(LOG_P0));
    y = _mm256_mul_ps(_mm256_mul_ps(y, f), z);
    y = _mm256_add_ps(y, _mm256_mul_ps(e, _mm256_set1_ps(LN2_C2)));
    y = _mm256_sub_ps(y, _mm256_mul_ps(z, _mm256_set1_ps(0.5f)));
    return _mm256_add_ps(_mm256_add_ps(f, y), _mm256_mul_ps(e,
        _mm256_set1_ps(LN2_C1)));
  }

  // f and e of a positive normal x with float instructions: 1+f is the
  // mantissa in [1,2), halved (and e raised by 1) above sqrt(2), and e comes
  // from converting the exponent bits as an integer (exactly)
  TA1 inline __m256 logReduceAvx(const __m256 x, __m256 &e)
  {
    __m256 m = _mm256_or_ps(_mm256_and_ps(x, _mm256_castsi256_ps(
        _mm256_set1_epi32(0x007fffff))), _mm256_set1_ps(1.0f));
    const __m256 big = _mm256_cmp_ps(m, _mm256_set1_ps(1.41421354f),
        _CMP_GT_OQ);
    e = _mm256_cvtepi32_ps(_mm256_castps_si256(_mm256_and_ps(x,
        _mm256_castsi256_ps(_mm256_set1_epi32(0x7f800000)))));
    e = _mm256_add_ps(_mm256_sub_ps(_mm256_mul_ps(e,
        _mm256_set1_ps(1.0f/8388608.0f)), _mm256_set1_ps(127.0f)),
        _mm256_and_ps(big, _mm256_set1_ps(1.0f)));
    m = _mm256_sub_ps(m, _mm256_and_ps(big, _mm256_mul_ps(m,
        _mm256_set1_ps(0.5f))));
    return _mm256_sub_ps(m, _mm256_set1_ps(1.0f));
  }

  // as logSpecialAvx2
  __attribute__((__target__("avx"), __noinline__))
  inline __m256 logSpecialAvx(const __m256 x)
  {
    const __m256 zero = _mm256_setzero_ps();
    const __m256 den = _mm256_cmp_ps(x, _mm256_set1_ps(1.17549435e-38f),
        _CMP_LT_OQ);
    __m256 e;
    const __m256 f = logReduceAvx(_mm256_mul_ps(x, sel256(den,
        _mm256_set1_ps(1.0f), _mm256_set1_ps(8388608.0f))), e);
    __m256 y = logPolyAvx(f, _mm256_sub_ps(e, _mm256_and_ps(den,
        _mm256_set1_ps(23.0f))));
    y = sel256(_mm256_cmp_ps(x, _mm256_set1_ps(__builtin_inff()),
        _CMP_EQ_OQ), y, x);
    y = sel256(_mm256_cmp_ps(x, zero, _CMP_EQ_OQ), y,
        _mm256_set1_ps(-__builtin_inff()));
    return _mm256_or_ps(y, _mm256_cmp_ps(x, zero, _CMP_NGE_UQ));
  }

  TA1 inline __m256 logAvx(const __m256 x)
  {
    // x below FLT_MIN (or NaN) or +INF
    if(_mm256_movemask_ps(_mm256_or_ps(_mm256_cmp_ps(x,
        _mm256_set1_ps(1.17549435e-38f), _CMP_NGE_UQ), _mm256_cmp_ps(x,
        _mm256_set1_ps(__builtin_inff()), _CMP_EQ_OQ))))
      return logSpecialAvx(x);
    __m256 e;
    const __m256 f = logReduceAvx(x, e);
    return logPolyAvx(f, e);
  }

  #undef TA5
  #undef TA2
  #undef TA1
} // namespace TRIG
} // namespace EVM
#endif // EVM_TRIG_H
