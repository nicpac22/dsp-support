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
 * 
 * Unlike most of my other software, this was vibe-coded.  Use with caution.
 */
//==============================================================================
// Name: evm_unit_test.cc
//
// Purpose: Correctness unit test for evm.h
//
// Created: 2026/09/27
//
// Description:
//  Randomly generates input vectors and checks every EVM function against a
//  double precision scalar baseline computed with std::complex<double>.  Each
//  function reports PASS or FAIL; on failure the length, element index, the
//  EVM result, the baseline result and both the absolute and relative
//  difference are printed.
//
//  Every function is exercised over a sweep of lengths rather than a single
//  length.  This matters because the SIMD kernels handle the final partial
//  register with separate remainder code, so a bug can be invisible at one
//  length and fatal at another (a masked-store defect that only appeared when
//  len%8==4 is exactly how this sweep earns its keep).  The sweep covers every
//  remainder for 4, 8 and 16 element registers, plus negative lengths, which
//  must be treated as a no-op.
//
//  Each function is also checked for writes past the end of the requested
//  output, and for correct operation in-place (the header documents that the
//  output buffer may alias an input buffer).
//
//  The random seed is printed on every run and may be supplied as argv[1] so
//  that a failing run can be reproduced exactly.
//
//  Tolerances are per algorithm class and are deliberately tight enough to
//  tell the classes apart: a "true divide" routine that silently lost its
//  divide and fell back to a bare reciprocal approximation would move from
//  ~1e-07 to ~2e-04 and fail TOL_EXACT by two orders of magnitude.
//
//  Exit status is 0 if every test passed, 1 otherwise.
//
//  Compile:
//    g++ -std=c++17 -march=native -O3 -I../include/ -o evm_unit_test \
//        evm_unit_test.cc
//    icpx -std=c++17 -xHost -O3 -I../include/ -o evm_unit_test evm_unit_test.cc
//
//  The runtime dispatch branch can be forced down with the usual macros, which
//  is worth doing since the branches are independent implementations:
//    g++ ... -DDISABLE_AVX512 ...   // exercise the AVX2 branch
//    g++ ... -DDISABLE_AVX2 ...     // exercise the AVX branch
//    g++ ... -DDISABLE_SIMD ...     // exercise the scalar branch
//==============================================================================
#include <algorithm>
#include <complex>
#include <vector>
#include <random>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <type_traits>
#if defined(__unix__)
#include <csetjmp>
#include <csignal>
#include <sys/mman.h>
#include <unistd.h>
#endif

#include "evm.h"

using std::complex;
typedef complex<double> cdbl;

//==============================================================================
// Configuration
//==============================================================================
// Longest vector exercised.  Input buffers are this long; output buffers carry
// an extra guard region so that writes past the end can be detected.
static const int MAXLEN = 1024;
static const int GUARD  = 32;
// Input length for the divide-gather functions, which read every stride-th
// element: enough for MAXLEN outputs at strides up to 16.
static const int GLEN   = 16*MAXLEN;
// Divisor floor passed to the eps overloads.  The divisors are in [-2,-0.25]
// and [0.25,2] (some exactly -0.5 or 0.5), so it raises every negative one and
// some of the positive ones, and ties with some.
static const float GEPS = 0.5f;

// Relative error tolerances by algorithm class.  The measured worst case for
// each class on AVX2/AVX/scalar is given for reference.
static const double TOL_EXACT = 1.0e-6;   // true multiply/divide   (~1.7e-07)
static const double TOL_NR    = 5.0e-6;   // reciprocal + 1x Newton (~2.5e-07)
static const double TOL_RCP   = 2.0e-3;   // bare reciprocal        (~2.9e-04)
static const double TOL_CMP   = 0.0;      // comparison/selection   (exact)
static const double TOL_TRIG  = 2.0e-6;   // sin/cos/atan2          (~8.1e-07)
static const double TOL_EXPLOG = 6.0e-7;  // exp/log, 5 ulp         (~2.5e-07)

// Number of input buffers for signedAbsMin/signedAbsMax.
static const int NBUF = 5;

// Most failing elements to print per function before summarising the rest.
static const int MAX_REPORT = 8;

//==============================================================================
// Random input data
//==============================================================================
static std::mt19937                 gRng;
static std::vector<float>           gAf, gBf;   // real input vectors
static std::vector<complex<float> > gAc, gBc;   // complex input vectors
static std::vector<complex<float> > gZc;        // gAc with zeros, see below
static float                        gSf;        // real scale factor
static complex<float>               gSc;        // complex scale factor
static std::vector<std::vector<float> > gNf;    // comparison inputs, see below
static float                        gTh;        // comparison threshold
static std::vector<float>           gCf;        // real third input (multSum)
static std::vector<complex<float> > gCc;        // complex third input
static std::vector<int8_t>          gI8;        // integer conversion inputs
static std::vector<int16_t>         gI16;
static std::vector<int32_t>         gI32;
static std::vector<complex<float> > gLc;        // gAc with large phases
static std::vector<complex<float> > gLr;
static std::vector<float>           gPf;        // |gBf|, positive divisors
static std::vector<float>           gGa, gGb;   // divide-gather inputs
static std::vector<float>           gEf, gLf;   // exp and log inputs

static float randReal()
{
  return std::uniform_real_distribution<float>(-2.0f, 2.0f)(gRng);
}

// Divisors are kept away from zero.  Division by (or very near) zero is a
// documented limitation of the reciprocal-approximation variants rather than a
// defect, so feeding it here would only produce noise.
static float randDivisorReal()
{
  float v = randReal();
  if(std::fabs(v) < 0.25f) v = (v < 0.0f) ? -0.5f : 0.5f;
  return v;
}

static complex<float> randComplex()
{
  return complex<float>(randReal(), randReal());
}

static complex<float> randDivisorComplex()
{
  complex<float> v = randComplex();
  const float m = std::abs(v);
  if(m < 0.25f) v = (m > 0.0f) ? v*(0.5f/m) : complex<float>(0.5f, 0.0f);
  return v;
}

// Comparison inputs are a multiple of 0.5 in [-2,2] half of the time, so ties
// are common: elements equal to the threshold, and equal magnitudes of
// opposite sign, which is where a > in place of >= or a tie going to the wrong
// buffer shows up.
static float randTie()
{
  if(std::uniform_int_distribution<int>(0, 1)(gRng))
    return 0.5f*std::uniform_int_distribution<int>(-4, 4)(gRng);
  return randReal();
}

static void generateData(unsigned int seed)
{
  gRng.seed(seed);
  gAf.resize(MAXLEN);  gBf.resize(MAXLEN);
  gAc.resize(MAXLEN);  gBc.resize(MAXLEN);
  gZc.resize(MAXLEN);
  for(int i=0; i<MAXLEN; ++i)
  {
    gAf[i] = randReal();
    gBf[i] = randDivisorReal();
    gAc[i] = randComplex();
    gBc[i] = randDivisorComplex();
    // input for the complex magnitudes with every third element zero (3 is
    // coprime to the register widths, so the zeros visit every lane): a
    // magnitude computed as x*rsqrt(x) gives 0*inf = NaN at zero unless the
    // kernel guards against it.  The zeros cycle through the four sign
    // combinations, since atan2(+-0, -0) is +-pi rather than +-0 (cartToPol)
    gZc[i] = (i%3 == 0) ? complex<float>((i&2) ? -0.0f : 0.0f,
                                         (i&4) ? -0.0f : 0.0f) : gAc[i];
  }
  gSf = randDivisorReal();
  gSc = randDivisorComplex();
  // drawn after the inputs above so that a seed still gives the same values
  // for them
  gNf.assign(NBUF, std::vector<float>(MAXLEN));
  for(int j=0; j<NBUF; ++j)
    for(int i=0; i<MAXLEN; ++i) gNf[j][i] = randTie();
  gTh = 0.5f*std::uniform_int_distribution<int>(-3, 3)(gRng);
  // also drawn last, for the same reason
  gCf.resize(MAXLEN);  gCc.resize(MAXLEN);
  for(int i=0; i<MAXLEN; ++i)
  {
    gCf[i] = randReal();
    gCc[i] = randComplex();
  }
  // also drawn last, for the same reason.  The integers cover each type's full
  // range; half of the 32-bit ones are within +/-2^24, where the conversion
  // is exact, and the rest are mostly larger, where it rounds
  gI8.resize(MAXLEN);  gI16.resize(MAXLEN);  gI32.resize(MAXLEN);
  for(int i=0; i<MAXLEN; ++i)
  {
    gI8[i]  = static_cast<int8_t>(std::uniform_int_distribution<int>(
        INT8_MIN, INT8_MAX)(gRng));
    gI16[i] = static_cast<int16_t>(std::uniform_int_distribution<int>(
        INT16_MIN, INT16_MAX)(gRng));
    gI32[i] = (i&1) ?
        std::uniform_int_distribution<int32_t>(INT32_MIN, INT32_MAX)(gRng) :
        std::uniform_int_distribution<int32_t>(-(1<<24), 1<<24)(gRng);
  }
  // polar inputs with large phases, derived from gAc (no draws): up to 1e4
  // cycles for polToCart, and up to 1e7 radians for polToCartRad, well past
  // the 65536 where evm_trig.h's reduction hands lanes over to scalar sincosf
  // (and far enough that the lanes would be wrong without it)
  gLc.resize(MAXLEN);  gLr.resize(MAXLEN);
  for(int i=0; i<MAXLEN; ++i)
  {
    gLc[i] = complex<float>(gAc[i].real(), 5000.0f*gAc[i].imag());
    gLr[i] = complex<float>(gAc[i].real(), 5.0e6f*gAc[i].imag());
  }
  // positive divisors for the square roots of divSqrt with eps = 0, derived
  // from gBf (no draws)
  gPf.resize(MAXLEN);
  for(int i=0; i<MAXLEN; ++i) gPf[i] = std::fabs(gBf[i]);
  // also drawn last, for the same reason as above
  gGa.resize(GLEN);  gGb.resize(GLEN);
  for(int i=0; i<GLEN; ++i)
  {
    gGa[i] = randReal();
    gGb[i] = randDivisorReal();
  }
  // exp and log inputs, derived from gAf and gBf (no draws).  exp's cover
  // [-88, 88], whose results go from denormals to near FLT_MAX.  log's are
  // |gBf| times 2^-124..2^124 (FLT_MIN to 2^125), and every fourth is within
  // 1/8 of 1, where the result is near 0
  gEf.resize(MAXLEN);  gLf.resize(MAXLEN);
  for(int i=0; i<MAXLEN; ++i)
  {
    gEf[i] = 44.0f*gAf[i];
    gLf[i] = (i%4 == 0) ? 1.0f + 0.0625f*gAf[i] :
        std::ldexp(std::fabs(gBf[i]), static_cast<int>(62.0f*gAf[i]));
  }
}

//==============================================================================
// Small helpers
//==============================================================================
template<class T> struct IsComplex            { static const bool value = false; };
template<> struct IsComplex<complex<float> >  { static const bool value = true;  };

static cdbl toCdbl(const float v)          { return cdbl(v, 0.0); }
static cdbl toCdbl(const complex<float> v) { return cdbl(v.real(), v.imag()); }

// The sentinel is a quiet NaN with a distinctive payload and is matched by its
// bits, so a NaN a kernel computes (0*inf, 0/0) is reported as a wrong value
// or a stray write rather than taken for an element that was never touched.
static const uint32_t SENTINEL_BITS = 0x7fc0deadu;

template<class T> static T sentinel();
template<> float sentinel<float>()
{
  float f;
  std::memcpy(&f, &SENTINEL_BITS, sizeof(f));
  return f;
}
template<> complex<float> sentinel<complex<float> >()
{
  const float n = sentinel<float>();
  return complex<float>(n, n);
}

// z raised to the power 2^squarings, the baseline for pow2/pow4/pow8
static cdbl powSq(cdbl z, const int squarings)
{
  for(int k=0; k<squarings; ++k) z *= z;
  return z;
}

// element k of the first nbufs comparison inputs with the smallest (or largest)
// magnitude, the earliest buffer winning a tie, the baseline for signedAbsMin
// and signedAbsMax
static cdbl signedAbsRef(const int k, const int nbufs, const bool wantMax)
{
  double val = gNf[0][k];
  for(int j=1; j<nbufs; ++j)
  {
    const double v = gNf[j][k];
    if(wantMax ? (std::fabs(v) > std::fabs(val))
               : (std::fabs(v) < std::fabs(val))) val = v;
  }
  return cdbl(val);
}

// (magnitude, phase*scale) to cartesian, the baseline for polToCart and
// polToCartRad (std::polar is undefined for a negative magnitude)
static cdbl polarRef(const complex<float> v, const double scale)
{
  const double m  = v.real();
  const double ph = scale*v.imag();
  return cdbl(m*std::cos(ph), m*std::sin(ph));
}

// (magnitude, phase in cycles) to cartesian with whole cycles removed exactly
// first, the baseline for polToCart at large phases
static cdbl polarRefCyc(const complex<float> v)
{
  const double m  = v.real();
  const double r  = v.imag() - std::nearbyint(static_cast<double>(v.imag()));
  return cdbl(m*std::cos(2.0*M_PI*r), m*std::sin(2.0*M_PI*r));
}

// zneg and zpos work in place.  They are run on a copy of the input whose
// elements past the end hold a value the function would change (-1 for zneg,
// 1 for zpos), since a stray write over the sentinel would put back the same
// NaN; the results, and any element past the end that changed, are passed on
// in out.
static void runZeroInPlace(const bool neg, const float* const in, const int len,
                           float* const out)
{
  const float guard = neg ? -1.0f : 1.0f;
  const int   valid = (len > 0) ? len : 0;
  std::vector<float> buf(MAXLEN + GUARD, guard);
  for(int k=0; k<valid; ++k) buf[k] = in[k];
  if(neg) EVM::zneg(buf.data(), len);
  else    EVM::zpos(buf.data(), len);
  for(int k=0; k<valid; ++k) out[k] = buf[k];
  for(int k=valid; k<MAXLEN+GUARD; ++k) if(buf[k] != guard) out[k] = buf[k];
}

// The integer to float conversions work in place through their void*
// wrappers, where the output is wider than the input.  The integers are copied
// to the front of a buffer that is otherwise the sentinel, and the whole
// buffer is passed on in out so that writes past the end show up.
template<class T>
static void runCnvrtInPlace(const T* const in, const int len, float* const out,
                            void (*cnvrt)(void*, int))
{
  std::vector<float> buf(MAXLEN + GUARD, sentinel<float>());
  if(len > 0) std::memcpy(buf.data(), in, len*sizeof(T));
  cnvrt(buf.data(), len);
  std::memcpy(out, buf.data(), buf.size()*sizeof(float));
}

static bool isSentinel(const float v)
{
  uint32_t b;
  std::memcpy(&b, &v, sizeof(b));
  return b == SENTINEL_BITS;
}
static bool isSentinel(const complex<float> v)
{
  return isSentinel(v.real()) && isSentinel(v.imag());
}

// Render a value the way its type reads: scalars plain, complex as (re,im).
static std::string format(const cdbl v, const bool cplx)
{
  char buf[96];
  if(cplx) snprintf(buf, sizeof(buf), "(%.9g,%.9g)", v.real(), v.imag());
  else     snprintf(buf, sizeof(buf), "%.9g", v.real());
  return std::string(buf);
}

// Lengths exercised.  Negative lengths must be a no-op; 0..64 covers every
// remainder for 4, 8 and 16 element SIMD registers several times over; the
// larger values check that nothing depends on the vector being short.
static std::vector<int> buildLengths()
{
  std::vector<int> v;
  const int neg[] = {-33, -17, -8, -3, -1};
  for(size_t i=0; i<sizeof(neg)/sizeof(neg[0]); ++i) v.push_back(neg[i]);
  for(int i=0; i<=64; ++i) v.push_back(i);
  const int big[] = {100, 127, 128, 129, 255, 256, 257, 1000, 1024};
  for(size_t i=0; i<sizeof(big)/sizeof(big[0]); ++i) v.push_back(big[i]);
  return v;
}
static const std::vector<int> gLengths = buildLengths();

static int gTests = 0, gFailed = 0;

//==============================================================================
// Test driver
//
//   runEvm(a, b, len, out)  calls the EVM function under test
//   baseline(k)             returns the expected value of element k as cdbl
//
// bBuf may be null for the single-input functions (scale), in which case the
// second pointer handed to runEvm is unused.
//
// Relative error is measured against max(|baseline|, FLT_MIN): float results
// below FLT_MIN are subnormal and carry only absolute precision (pow8 of a
// small input lands there).  With flush to zero on (icpx turns it on by
// default) such a result comes out as 0, which is accepted as well.
//
// magnitude(k), if given, takes the place of |baseline|.  multSum passes
// |in1*in2| + |in3|: its result can cancel to near zero, where the rounding
// error of the product is large next to the result but not next to the terms
// (the dot-products are measured the same way).
//
// overComplex = false leaves out mode 3 (below) for the functions that add
// to the output buffer, which may not be on top of the complex input.
//==============================================================================
template<class TA, class TB, class TO, class FEVM, class FREF,
         class FMAG = std::nullptr_t>
static void check(const char* name, const double tol,
                  const TA* const aBuf, const TB* const bBuf,
                  FEVM runEvm, FREF baseline, FMAG magnitude = nullptr,
                  const bool overComplex = true)
{
  ++gTests;
  const bool cplx = IsComplex<TO>::value;

  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0;
  double worstRel = 0.0;

  std::vector<TO> out(MAXLEN + GUARD);
  // flush to zero, under which a result below FLT_MIN may come out as 0
  const bool ftz = (_mm_getcsr() & _MM_FLUSH_ZERO_MASK) == _MM_FLUSH_ZERO_ON;
  // complex input buffer whose front receives the real results in mode 3
  const bool cplxToReal = IsComplex<TA>::value && !cplx;
  std::vector<TA> io(MAXLEN + GUARD);

  // mode 0: separate buffers.  mode 1: out aliases in1.  mode 2: out aliases
  // in2.  mode 3: a real output overwrites the front of a complex in1 (the
  // magnitude functions).  The header states the output buffer may be the
  // same as an input buffer, so the aliasing modes are part of the contract,
  // not a bonus.
  for(int mode=0; mode<4; ++mode)
  {
    if(mode == 1 && !std::is_same<TA, TO>::value) continue;
    if(mode == 2 && (!std::is_same<TB, TO>::value || bBuf == NULL)) continue;
    if(mode == 3 && (!cplxToReal || !overComplex)) continue;
    const char* modeName = (mode==0) ? "" :
                           (mode==1) ? " [in-place on in1]" :
                           (mode==2) ? " [in-place on in2]" :
                                       " [in-place over complex in1]";

    for(size_t li=0; li<gLengths.size(); ++li)
    {
      const int len   = gLengths[li];
      const int valid = (len > 0) ? len : 0;

      // results and how many of them to check for stray writes
      TO* res = out.data();
      int resLen = MAXLEN + GUARD;
      for(int k=0; k<MAXLEN+GUARD; ++k) out[k] = sentinel<TO>();

      if(mode == 3)
      {
        if constexpr (IsComplex<TA>::value && !IsComplex<TO>::value)
        {
          for(int k=0; k<MAXLEN+GUARD; ++k)
            io[k] = (k < valid) ? aBuf[k] : sentinel<TA>();
          res = reinterpret_cast<TO*>(io.data());
          resLen = 2*(MAXLEN + GUARD);
          runEvm(io.data(), bBuf, len, res);
        }
      }
      else if(mode == 0)
      {
        runEvm(aBuf, bBuf, len, out.data());
      }
      else if(mode == 1)
      {
        // seed the output with in1, then ask the function to overwrite it
        if constexpr (std::is_same<TA, TO>::value)
        {
          for(int k=0; k<valid; ++k) out[k] = aBuf[k];
          runEvm(reinterpret_cast<const TA*>(out.data()), bBuf, len, out.data());
        }
      }
      else
      {
        if constexpr (std::is_same<TB, TO>::value)
        {
          for(int k=0; k<valid; ++k) out[k] = bBuf[k];
          runEvm(aBuf, reinterpret_cast<const TB*>(out.data()), len, out.data());
        }
      }

      // --- values ---
      for(int k=0; k<valid; ++k)
      {
        const cdbl   got   = toCdbl(res[k]);
        const cdbl   want  = baseline(k);
        const double adiff = std::abs(got - want);
        double size = std::abs(want);
        if constexpr (!std::is_same<FMAG, std::nullptr_t>::value)
          size = magnitude(k);
        const double mag   = std::max(size,
            static_cast<double>(std::numeric_limits<float>::min()));
        const double rdiff = adiff/mag;
        ++compared;

        if(rdiff <= tol)
        {
          if(rdiff > worstRel) worstRel = rdiff;
          continue;
        }
        if(ftz && got == cdbl(0.0, 0.0) &&
           std::abs(want) < std::numeric_limits<float>::min()) continue;
        ++nbad;
        if(nshown < MAX_REPORT)
        {
          ++nshown;
          if(isSentinel(res[k]))
          {
            printf("      len=%-5d element %-5d NOT WRITTEN   baseline=%s%s\n",
                   len, k, format(want, cplx).c_str(), modeName);
          }
          else
          {
            printf("      len=%-5d element %-5d evm=%s baseline=%s "
                   "absdiff=%.6g reldiff=%.6g%s\n",
                   len, k, format(got, cplx).c_str(),
                   format(want, cplx).c_str(), adiff, rdiff, modeName);
          }
        }
      }

      // --- writes past the requested length (covers negative len too) ---
      // in mode 3 the floats after the results still hold the unused part of
      // the complex input, which must be left as it was
      const float* inFloats = cplxToReal ?
          reinterpret_cast<const float*>(aBuf) : NULL;
      for(int k=valid; k<resLen; ++k)
      {
        if(mode == 3 && k < 2*valid)
        {
          if(std::memcmp(&res[k], &inFloats[k], sizeof(float)) == 0) continue;
        }
        else if(isSentinel(res[k])) continue;
        ++nstray;
        if(nshown < MAX_REPORT)
        {
          ++nshown;
          printf("      len=%-5d element %-5d WROTE PAST END evm=%s%s\n",
                 len, k, format(toCdbl(res[k]), cplx).c_str(), modeName);
        }
      }
    }
  }

  const bool pass = (nbad == 0 && nstray == 0);
  if(!pass) ++gFailed;

  printf("  %-34s %s  [%ld values, tol %.0e, worst rel %.3e",
         name, pass ? "PASS" : "FAIL", compared, tol, worstRel);
  if(nbad)   printf(", %d wrong", nbad);
  if(nstray) printf(", %d past end", nstray);
  printf("]\n");
  if(nbad + nstray > nshown)
    printf("      ... %d further problems not shown\n",
           nbad + nstray - nshown);
}

//==============================================================================
// Reduction test driver (dot-products)
//
//   runEvm(a, b, len, res)  calls the EVM reduction under test and puts its
//                           result in res[0]; res[1] is a guard that must be
//                           left untouched
//   conjB                   true for the conjugate dot-products, whose
//                           baseline uses conj(b[k])
//   revB                    true for the reversed dot-products (dotpRev,
//                           dotpcRev), whose baseline pairs a[k] with
//                           b[len-1-k]
//
// The baseline is the double precision sum of a[k]*b[k].  Summed in float in
// any order, each component of the result is within (m+1)*2^-24 times the sum
// of the magnitudes of its m product terms (the standard forward error bound;
// m = 2*len when both inputs are complex, whose real and imaginary parts each
// add two products per element, otherwise m = len).  The error is measured
// against that magnitude sum rather than the result, which cancellation can
// bring near zero.  The inputs are copied into NaN-filled buffers, in2 one
// element off in1's alignment, so a read past the end shows up as a NaN
// result (for the reversed ones, a read before the start of in2 too, since
// the element before it is NaN as well).  Negative and zero lengths must give
// exactly 0.
//==============================================================================
template<class TA, class TB, class FEVM>
static void checkDot(const char* name, const TA* const aBuf,
                     const TB* const bBuf, FEVM runEvm,
                     const bool conjB = false, const bool revB = false)
{
  typedef decltype(TA()*TB()) T;       // result type
  ++gTests;
  const bool   cplx  = IsComplex<T>::value;
  // products per element per component
  const int    terms = (IsComplex<TA>::value && IsComplex<TB>::value) ? 2 : 1;
  const double u     = std::ldexp(1.0, -24);
  std::vector<TA> a(MAXLEN + GUARD);
  std::vector<TB> b(MAXLEN + GUARD);

  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0;
  double worstRel = 0.0;

  for(size_t li=0; li<gLengths.size(); ++li)
  {
    const int len   = gLengths[li];
    const int valid = (len > 0) ? len : 0;

    for(int k=0; k<MAXLEN+GUARD; ++k)
    {
      a[k] = sentinel<TA>();
      b[k] = sentinel<TB>();
    }
    for(int k=0; k<valid; ++k) { a[k] = aBuf[k]; b[k+1] = bBuf[k]; }
    T res[2] = {sentinel<T>(), sentinel<T>()};
    runEvm(a.data(), b.data()+1, len, res);
    const cdbl got = toCdbl(res[0]);

    // baseline, and the magnitude sums that scale each component's error
    cdbl   want  = 0.0;
    double absRe = 0.0, absIm = 0.0;
    for(int k=0; k<valid; ++k)
    {
      const cdbl x = toCdbl(aBuf[k]);
      const cdbl z = toCdbl(bBuf[revB ? valid-1-k : k]);
      const cdbl y = conjB ? std::conj(z) : z;
      want  += x*y;
      absRe += std::fabs(x.real()*y.real()) + std::fabs(x.imag()*y.imag());
      absIm += std::fabs(x.real()*y.imag()) + std::fabs(x.imag()*y.real());
    }
    const double tol   = (terms*valid + 1)*u;
    const double dRe   = std::fabs(got.real() - want.real());
    const double dIm   = std::fabs(got.imag() - want.imag());
    const double rdiff = std::max((absRe > 0.0) ? dRe/absRe : dRe,
                                  (absIm > 0.0) ? dIm/absIm : dIm);
    ++compared;

    // --- writes past the single result ---
    if(!isSentinel(res[1]))
    {
      ++nstray;
      if(nshown < MAX_REPORT)
      {
        ++nshown;
        printf("      len=%-5d WROTE PAST THE RESULT evm=%s\n", len,
               format(toCdbl(res[1]), cplx).c_str());
      }
    }

    // --- value (NaN fails both tests) ---
    if((valid > 0) ? (rdiff <= tol) : (got == cdbl(0.0, 0.0)))
    {
      if(rdiff > worstRel) worstRel = rdiff;
      continue;
    }
    ++nbad;
    if(nshown < MAX_REPORT)
    {
      ++nshown;
      printf("      len=%-5d evm=%s baseline=%s reldiff=%.6g "
             "(of sum|a*b|, tol %.3g)\n", len, format(got, cplx).c_str(),
             format(want, cplx).c_str(), rdiff, tol);
    }
  }

  const bool pass = (nbad == 0 && nstray == 0);
  if(!pass) ++gFailed;

  printf("  %-34s %s  [%ld values, tol (%slen+1)*2^-24, worst rel %.3e",
         name, pass ? "PASS" : "FAIL", compared, (terms == 2) ? "2*" : "",
         worstRel);
  if(nbad)   printf(", %d wrong", nbad);
  if(nstray) printf(", %d past end", nstray);
  printf("]\n");
  if(nbad + nstray > nshown)
    printf("      ... %d further problems not shown\n",
           nbad + nstray - nshown);
}

//==============================================================================
// Normalized dot-product test driver (dotpn, dotpnr, dotpcn, dotpcnr)
//
//   runEvm(a, b, len, res)  as for checkDot()
//   conjB                   true for dotpcn/dotpcnr, whose products use
//                           conj(b[k])
//   tol                     error allowed in each normalized product
//                           (TOL_EXACT, or TOL_RCP for the approximations)
//
// The baseline is the double precision sum of p/|p| over the products
// p = a[k]*b[k], skipping those whose squared magnitude is below FLT_MIN (a
// zero product adds nothing).  Each term has magnitude 1, so the error of each
// component is measured against the number of terms n: the bound is
// tol + (len+1)*2^-24 times n, the error of each term plus that of summing
// them in float.  The inputs are copied into NaN-filled buffers as for
// checkDot(), so a read past the end shows up as a NaN result.  Negative and
// zero lengths must give exactly 0.
//==============================================================================
template<class FEVM>
static void checkDotN(const char* name, const complex<float>* const aBuf,
                      const complex<float>* const bBuf, FEVM runEvm,
                      const bool conjB, const double tol)
{
  ++gTests;
  const double u = std::ldexp(1.0, -24);
  std::vector<complex<float> > a(MAXLEN + GUARD);
  std::vector<complex<float> > b(MAXLEN + GUARD);

  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0;
  double worstRel = 0.0;

  for(size_t li=0; li<gLengths.size(); ++li)
  {
    const int len   = gLengths[li];
    const int valid = (len > 0) ? len : 0;

    for(int k=0; k<MAXLEN+GUARD; ++k)
    {
      a[k] = sentinel<complex<float> >();
      b[k] = sentinel<complex<float> >();
    }
    for(int k=0; k<valid; ++k) { a[k] = aBuf[k]; b[k+1] = bBuf[k]; }
    complex<float> res[2] = {sentinel<complex<float> >(),
                             sentinel<complex<float> >()};
    runEvm(a.data(), b.data()+1, len, res);
    const cdbl got = toCdbl(res[0]);

    // baseline and the number of terms
    cdbl want = 0.0;
    int  terms = 0;
    for(int k=0; k<valid; ++k)
    {
      const cdbl y = conjB ? std::conj(toCdbl(bBuf[k])) : toCdbl(bBuf[k]);
      const cdbl p = toCdbl(aBuf[k])*y;
      const double m2 = std::norm(p);
      if(m2 >= FLT_MIN) { want += p/std::sqrt(m2); ++terms; }
    }
    const double scale = (terms > 0) ? terms : 1.0;
    const double bound = tol + (valid + 1)*u;
    const double rdiff = std::max(std::fabs(got.real() - want.real()),
                                  std::fabs(got.imag() - want.imag()))/scale;
    ++compared;

    // --- writes past the single result ---
    if(!isSentinel(res[1]))
    {
      ++nstray;
      if(nshown < MAX_REPORT)
      {
        ++nshown;
        printf("      len=%-5d WROTE PAST THE RESULT evm=%s\n", len,
               format(toCdbl(res[1]), true).c_str());
      }
    }

    // --- value (NaN fails both tests) ---
    if((valid > 0) ? (rdiff <= bound) : (got == cdbl(0.0, 0.0)))
    {
      if(rdiff > worstRel) worstRel = rdiff;
      continue;
    }
    ++nbad;
    if(nshown < MAX_REPORT)
    {
      ++nshown;
      printf("      len=%-5d evm=%s baseline=%s reldiff=%.6g "
             "(of %d terms, tol %.3g)\n", len, format(got, true).c_str(),
             format(want, true).c_str(), rdiff, terms, bound);
    }
  }

  const bool pass = (nbad == 0 && nstray == 0);
  if(!pass) ++gFailed;

  printf("  %-34s %s  [%ld values, tol %.0e+(len+1)*2^-24, worst rel %.3e",
         name, pass ? "PASS" : "FAIL", compared, tol, worstRel);
  if(nbad)   printf(", %d wrong", nbad);
  if(nstray) printf(", %d past end", nstray);
  printf("]\n");
  if(nbad + nstray > nshown)
    printf("      ... %d further problems not shown\n",
           nbad + nstray - nshown);
}

//==============================================================================
// Expand test driver (multe, dive, divre, divnre)
//
//   runEvm(a, b, len, out)  calls the EVM function under test
//   baseline(k)             returns the expected real part of element k
//
// These put in1*in2 (or in1/in2) in the real part of each complex output
// element and leave the imaginary part as it was.  The output buffer starts as
// the sentinel, so every imaginary part must still be the sentinel afterwards,
// as must every float past the end.  In place, the complex output is laid over
// in1 or in2, which it is twice as wide as (the kernels run from the end to the
// start for this), and each imaginary part must still hold the input float
// that was there.  The real parts must be within tol of the double baseline.
//==============================================================================
template<class FEVM, class FREF>
static void checkExpand(const char* name, const double tol,
                        const float* const aBuf, const float* const bBuf,
                        FEVM runEvm, FREF baseline)
{
  ++gTests;
  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0;
  double worstRel = 0.0;

  const int nfl = 2*(MAXLEN + GUARD);
  std::vector<float> buf(nfl), before(nfl);

  // mode 0: separate buffers.  mode 1: out over in1.  mode 2: out over in2.
  for(int mode=0; mode<3; ++mode)
  {
    const char* modeName = (mode==0) ? "" :
                           (mode==1) ? " [in-place on in1]" :
                                       " [in-place on in2]";
    for(size_t li=0; li<gLengths.size(); ++li)
    {
      const int len   = gLengths[li];
      const int valid = (len > 0) ? len : 0;

      for(int k=0; k<nfl; ++k) buf[k] = sentinel<float>();
      const float* a = aBuf;
      const float* b = bBuf;
      if(mode == 1)
      {
        for(int k=0; k<valid; ++k) buf[k] = aBuf[k];
        a = buf.data();
      }
      if(mode == 2)
      {
        for(int k=0; k<valid; ++k) buf[k] = bBuf[k];
        b = buf.data();
      }
      before = buf;
      runEvm(a, b, len, reinterpret_cast<complex<float>*>(buf.data()));

      for(int k=0; k<nfl; ++k)
      {
        const int e = k>>1; // complex element
        if(e < valid && !(k&1)) // real part, must hold the result
        {
          const double want  = baseline(e);
          const double adiff = std::abs(static_cast<double>(buf[k]) - want);
          const double mag   = std::max(std::abs(want),
              static_cast<double>(std::numeric_limits<float>::min()));
          const double rdiff = adiff/mag;
          ++compared;
          if(rdiff <= tol)
          {
            if(rdiff > worstRel) worstRel = rdiff;
            continue;
          }
          ++nbad;
          if(nshown < MAX_REPORT)
          {
            ++nshown;
            if(std::memcmp(&buf[k], &before[k], sizeof(float)) == 0)
              printf("      len=%-5d element %-5d NOT WRITTEN   "
                     "baseline=%.9g%s\n", len, e, want, modeName);
            else
              printf("      len=%-5d element %-5d evm=%.9g baseline=%.9g "
                     "absdiff=%.6g reldiff=%.6g%s\n",
                     len, e, buf[k], want, adiff, rdiff, modeName);
          }
        }
        // imaginary part or past the end, must be left as it was
        else if(std::memcmp(&buf[k], &before[k], sizeof(float)) != 0)
        {
          const bool imag = (e < valid);
          if(imag) ++nbad;
          else     ++nstray;
          if(nshown < MAX_REPORT)
          {
            ++nshown;
            printf("      len=%-5d element %-5d %s evm=%.9g was=%.9g%s\n",
                   len, e, imag ? "IMAG CHANGED" : "WROTE PAST END",
                   buf[k], before[k], modeName);
          }
        }
      }
    }
  }

  const bool pass = (nbad == 0 && nstray == 0);
  if(!pass) ++gFailed;

  printf("  %-34s %s  [%ld values, tol %.0e, worst rel %.3e",
         name, pass ? "PASS" : "FAIL", compared, tol, worstRel);
  if(nbad)   printf(", %d wrong", nbad);
  if(nstray) printf(", %d past end", nstray);
  printf("]\n");
  if(nbad + nstray > nshown)
    printf("      ... %d further problems not shown\n",
           nbad + nstray - nshown);
}

//==============================================================================
// Divide-gather test driver (divg, divrg, divnrg)
//
//   runEvm(a, b, len, stride, out)  calls the EVM function under test
//   ref(a, b)                       returns the expected quotient of the input
//                                   elements a and b, as a double
//
// These compute out[k] = in1[k*stride]/in2[k*stride] for the
// n = ceil(len/stride) outputs, on gGa and gGb.  Each stride runs the output
// counts of the length sweep, each with the shortest and the longest input
// that gives that many outputs (the last output comes from the last input
// element, or from stride-1 elements before it).  Negative and zero lengths
// and strides must be a no-op.  The output runs in its own buffer and on top
// of in1 or in2, and must not change anything past out[n-1]: in its own
// buffer that is the sentinel, in place it is the rest of the input.  The two
// largest strides only give one output, and put the offsets of the unused
// lanes of the register past 2^31, where they wrap around.  The kernels load
// small strides as whole registers, so on unix-like systems a fourth mode
// puts the end of each input right before an inaccessible page: a read past
// the end faults, and the fault is caught and reported.
//==============================================================================
#if defined(__unix__)
// input buffers that end at an inaccessible page
static sigjmp_buf gFaultJmp;
static void onFault(int) { siglongjmp(gFaultJmp, 1); }
struct GuardedBuf
{
  char*  base;
  size_t bytes;
  float* end;  // first float of the inaccessible page
  explicit GuardedBuf(const size_t nfloats)
  {
    const size_t page = static_cast<size_t>(sysconf(_SC_PAGESIZE));
    const size_t data = (nfloats*sizeof(float) + page - 1)/page*page;
    bytes = data + page;
    base  = static_cast<char*>(mmap(NULL, bytes, PROT_READ | PROT_WRITE,
                                    MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
    if(base == MAP_FAILED || mprotect(base + data, page, PROT_NONE) != 0)
    {
      printf("mmap/mprotect failed\n");
      exit(1);
    }
    end = reinterpret_cast<float*>(base + data);
  }
  ~GuardedBuf() { munmap(base, bytes); }
};
static const int GATHER_MODES = 4;
#else
static const int GATHER_MODES = 3;
#endif

template<class FEVM, class FREF>
static void checkGather(const char* name, const double tol, FEVM runEvm,
                        FREF ref)
{
  ++gTests;
  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0, nfault = 0;
  double worstRel = 0.0;

  // (len, stride) cases
  std::vector<std::pair<int, int> > cases;
  const int strides[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 16, 17, 33,
                         64, 1000, 1<<28, std::numeric_limits<int>::max()};
  for(const int s : strides)
    for(size_t li=0; li<gLengths.size(); ++li)
    {
      const long long n = gLengths[li];
      if(n <= 0) { cases.push_back(std::make_pair(gLengths[li], s)); continue; }
      const long long shortest = (n-1)*s + 1, longest = n*s;
      if(shortest <= GLEN) cases.push_back(std::make_pair(int(shortest), s));
      if(longest <= GLEN && longest != shortest)
        cases.push_back(std::make_pair(int(longest), s));
    }
  const int noop[][2] = {{100, 0}, {100, -1}, {100, -16}, {-5, 0}};
  for(size_t c=0; c<sizeof(noop)/sizeof(noop[0]); ++c)
    cases.push_back(std::make_pair(noop[c][0], noop[c][1]));
  cases.push_back(std::make_pair(GLEN, 1<<28));
  cases.push_back(std::make_pair(GLEN, std::numeric_limits<int>::max()));

  std::vector<float> buf(GLEN + GUARD), before(GLEN + GUARD);
  for(int k=0; k<GLEN+GUARD; ++k) buf[k] = sentinel<float>();
#if defined(__unix__)
  GuardedBuf ga(GLEN), gb(GLEN);
#endif

  // mode 0: separate buffers.  mode 1: out over in1.  mode 2: out over in2.
  // mode 3: separate buffers, the inputs end at an inaccessible page.
  for(int mode=0; mode<GATHER_MODES; ++mode)
  {
    const char* modeName = (mode==0) ? "" :
                           (mode==1) ? " [in-place on in1]" :
                           (mode==2) ? " [in-place on in2]" :
                                       " [inputs end at a guard page]";
    for(size_t c=0; c<cases.size(); ++c)
    {
      const int len    = cases[c].first;
      const int stride = cases[c].second;
      const int valid  = (len > 0) ? len : 0;
      const int n      = (len > 0 && stride > 0) ? (len-1)/stride + 1 : 0;
      const int span   = valid + GUARD; // floats checked

      for(int k=0; k<span; ++k) buf[k] = sentinel<float>();
      const float* a = gGa.data();
      const float* b = gGb.data();
      if(mode == 1)
      {
        for(int k=0; k<valid; ++k) buf[k] = gGa[k];
        a = buf.data();
      }
      if(mode == 2)
      {
        for(int k=0; k<valid; ++k) buf[k] = gGb[k];
        b = buf.data();
      }
      std::copy(buf.begin(), buf.begin()+span, before.begin());
#if defined(__unix__)
      if(mode == 3)
      {
        float* const ca = ga.end - valid;
        float* const cb = gb.end - valid;
        std::copy(gGa.begin(), gGa.begin()+valid, ca);
        std::copy(gGb.begin(), gGb.begin()+valid, cb);
        struct sigaction sa, old;
        std::memset(&sa, 0, sizeof(sa));
        sa.sa_handler = onFault;
        sigaction(SIGSEGV, &sa, &old);
        const bool fault = (sigsetjmp(gFaultJmp, 1) != 0);
        if(!fault) runEvm(ca, cb, len, stride, buf.data());
        sigaction(SIGSEGV, &old, NULL);
        if(fault)
        {
          ++nfault;
          if(nshown < MAX_REPORT)
          {
            ++nshown;
            printf("      len=%-5d stride=%-10d READ PAST THE END OF THE "
                   "INPUT%s\n", len, stride, modeName);
          }
          continue;
        }
      }
      else
#endif
      runEvm(a, b, len, stride, buf.data());

      for(int k=0; k<span; ++k)
      {
        if(k < n) // must hold the quotient
        {
          const size_t j     = static_cast<size_t>(k)*stride;
          const double want  = ref(gGa[j], gGb[j]);
          const double adiff = std::abs(static_cast<double>(buf[k]) - want);
          const double mag   = std::max(std::abs(want),
              static_cast<double>(std::numeric_limits<float>::min()));
          const double rdiff = adiff/mag;
          ++compared;
          if(rdiff <= tol)
          {
            if(rdiff > worstRel) worstRel = rdiff;
            continue;
          }
          ++nbad;
          if(nshown < MAX_REPORT)
          {
            ++nshown;
            if(std::memcmp(&buf[k], &before[k], sizeof(float)) == 0)
              printf("      len=%-5d stride=%-10d element %-5d NOT WRITTEN   "
                     "baseline=%.9g%s\n", len, stride, k, want, modeName);
            else
              printf("      len=%-5d stride=%-10d element %-5d evm=%.9g "
                     "baseline=%.9g reldiff=%.6g%s\n",
                     len, stride, k, buf[k], want, rdiff, modeName);
          }
        }
        // past the last output, must be left as it was
        else if(std::memcmp(&buf[k], &before[k], sizeof(float)) != 0)
        {
          ++nstray;
          if(nshown < MAX_REPORT)
          {
            ++nshown;
            printf("      len=%-5d stride=%-10d element %-5d WROTE PAST END "
                   "evm=%.9g was=%.9g%s\n",
                   len, stride, k, buf[k], before[k], modeName);
          }
        }
      }
    }
  }

  const bool pass = (nbad == 0 && nstray == 0 && nfault == 0);
  if(!pass) ++gFailed;

  printf("  %-34s %s  [%ld values, tol %.0e, worst rel %.3e",
         name, pass ? "PASS" : "FAIL", compared, tol, worstRel);
  if(nbad)   printf(", %d wrong", nbad);
  if(nstray) printf(", %d past end", nstray);
  if(nfault) printf(", %d read past the input", nfault);
  printf("]\n");
  if(nbad + nstray + nfault > nshown)
    printf("      ... %d further problems not shown\n",
           nbad + nstray + nfault - nshown);
}

//==============================================================================
// Two-output test driver (minMax of two buffers)
//
//   runEvm(a, b, len, o1, o2)  calls the EVM function under test
//   base1(k), base2(k)         return the expected values of o1[k] and o2[k]
//
// Both outputs must match the baseline exactly.  Each output is run in its own
// buffer and on top of in1 or in2 (the two outputs never share a buffer): a
// kernel that stores the first output before it has read everything the
// second one needs only fails when they alias.  All four buffers are checked
// for writes past the requested length.  Past the end in1 holds 1e30 and in2
// holds -1e30 rather than the sentinel, since the min or max of two sentinels
// is the sentinel again and a stray write of it would go unseen.
//==============================================================================
template<class FEVM, class F1, class F2>
static void checkPair(const char* name, const float* const aBuf,
                      const float* const bBuf, FEVM runEvm, F1 base1, F2 base2)
{
  ++gTests;
  long compared = 0;
  int  nbad = 0, nshown = 0, nstray = 0;

  // buffers 0 and 1 hold in1 and in2, 2 and 3 are the outputs' own buffers
  std::vector<std::vector<float> > buf(4, std::vector<float>(MAXLEN + GUARD));
  const float fill[4] = {1.0e30f, -1.0e30f, sentinel<float>(),
                         sentinel<float>()};
  const char* const where[]   = {"in1", "in2", "own buffer"};
  const char* const bufName[] = {"in1", "in2", "out1", "out2"};

  // o1 and o2 each go on top of in1 (0) or in2 (1), or in their own buffer (2)
  for(int m1=0; m1<3; ++m1)
  for(int m2=0; m2<3; ++m2)
  {
    if(m1 == m2 && m1 != 2) continue;
    float* const o1 = buf[(m1 == 2) ? 2 : m1].data();
    float* const o2 = buf[(m2 == 2) ? 3 : m2].data();

    for(size_t li=0; li<gLengths.size(); ++li)
    {
      const int len   = gLengths[li];
      const int valid = (len > 0) ? len : 0;
      for(int k=0; k<MAXLEN+GUARD; ++k)
      {
        buf[0][k] = (k < valid) ? aBuf[k] : fill[0];
        buf[1][k] = (k < valid) ? bBuf[k] : fill[1];
        buf[2][k] = fill[2];
        buf[3][k] = fill[3];
      }
      runEvm(buf[0].data(), buf[1].data(), len, o1, o2);

      // --- values ---
      for(int k=0; k<valid; ++k)
      {
        const double got[2]  = {o1[k], o2[k]};
        const double want[2] = {base1(k).real(), base2(k).real()};
        for(int n=0; n<2; ++n)
        {
          ++compared;
          if(got[n] == want[n]) continue;
          ++nbad;
          if(nshown < MAX_REPORT)
          {
            ++nshown;
            printf("      len=%-5d element %-5d out%d=%.9g baseline=%.9g "
                   "[out1 on %s, out2 on %s]\n", len, k, n+1, got[n],
                   want[n], where[m1], where[m2]);
          }
        }
      }

      // --- writes past the requested length (covers negative len too) ---
      for(int b=0; b<4; ++b)
      {
        for(int k=valid; k<MAXLEN+GUARD; ++k)
        {
          if(std::memcmp(&buf[b][k], &fill[b], sizeof(float)) == 0) continue;
          ++nstray;
          if(nshown < MAX_REPORT)
          {
            ++nshown;
            printf("      len=%-5d element %-5d WROTE PAST END of %s "
                   "evm=%.9g [out1 on %s, out2 on %s]\n", len, k,
                   bufName[b], buf[b][k], where[m1], where[m2]);
          }
        }
      }
    }
  }

  const bool pass = (nbad == 0 && nstray == 0);
  if(!pass) ++gFailed;

  printf("  %-34s %s  [%ld values, exact", name, pass ? "PASS" : "FAIL",
         compared);
  if(nbad)   printf(", %d wrong", nbad);
  if(nstray) printf(", %d past end", nstray);
  printf("]\n");
  if(nbad + nstray > nshown)
    printf("      ... %d further problems not shown\n",
           nbad + nstray - nshown);
}

//==============================================================================
// minMax reduction test driver
//
// The result must equal the baseline exactly.  Every length is run with noInit
// false, and with noInit true from starting values inside, outside and the
// wrong way round relative to the data range, on the data as given, made all
// positive and made all negative: a kernel that let the zeros a masked load
// leaves in the unused lanes take part would return 0 as the minimum of an
// all-positive buffer or the maximum of an all-negative one.  Then a single
// minimum, and then a single maximum, is put at every position of every
// length, so that each lane of every accumulator and of the partial register
// has to be the one that comes out.  The input is copied into a buffer whose
// elements past the end alternate between -1e30 and 1e30, so a read past the
// end changes the result.  Negative and zero lengths must leave minOut and
// maxOut untouched.
//==============================================================================
static void checkMinMax(const char* name, const float* const aBuf)
{
  ++gTests;
  long compared = 0;
  int  nbad = 0, nshown = 0;
  std::vector<float> a(MAXLEN + GUARD);

  // starting values when noInit is true, case 0 runs with noInit false
  const float initMin[] = {0.0f, -0.5f, -3.0f,  3.0f};
  const float initMax[] = {0.0f,  0.5f,  3.0f, -3.0f};
  const char* const dataName[] = {"as given", "all positive", "all negative",
                                  "single min/max"};

  // runs minMax over a with noInit from case init and checks the result
  auto runCase = [&](const int len, const int init, const int data,
                     const int at)
  {
    const int  valid  = (len > 0) ? len : 0;
    const bool noInit = (init != 0);
    // without noInit the outputs start as the sentinel, which a negative or
    // zero length must leave in place
    float mn = noInit ? initMin[init] : sentinel<float>();
    float mx = noInit ? initMax[init] : sentinel<float>();
    float wantMin = mn;
    float wantMax = mx;
    if(valid > 0)
    {
      wantMin = noInit ? std::min(initMin[init], a[0]) : a[0];
      wantMax = noInit ? std::max(initMax[init], a[0]) : a[0];
      for(int k=1; k<valid; ++k)
      {
        wantMin = std::min(wantMin, a[k]);
        wantMax = std::max(wantMax, a[k]);
      }
    }
    EVM::minMax(a.data(), len, mn, mx, noInit);
    ++compared;

    const bool ok = (valid == 0 && !noInit) ?
        (isSentinel(mn) && isSentinel(mx)) :
        (mn == wantMin && mx == wantMax);
    if(ok) return;
    ++nbad;
    if(nshown < MAX_REPORT)
    {
      ++nshown;
      printf("      len=%-5d evm min,max=%.9g,%.9g baseline=%.9g,%.9g "
             "[%s", len, mn, mx, wantMin, wantMax, dataName[data]);
      if(at >= 0) printf(" at element %d", at);
      printf(", noInit=%d]\n", noInit ? 1 : 0);
    }
  };

  for(size_t li=0; li<gLengths.size(); ++li)
  {
    const int len   = gLengths[li];
    const int valid = (len > 0) ? len : 0;
    for(int data=0; data<4; ++data)
    {
      for(int k=0; k<MAXLEN+GUARD; ++k)
      {
        if(k >= valid) a[k] = (k&1) ? 1.0e30f : -1.0e30f;
        else if(data == 1) a[k] = std::fabs(aBuf[k]) + 0.5f;
        else if(data == 2) a[k] = -std::fabs(aBuf[k]) - 0.5f;
        else               a[k] = aBuf[k];
      }
      if(data < 3)
      {
        for(int init=0; init<4; ++init) runCase(len, init, data, -1);
        continue;
      }
      // a minimum then a maximum outside the data range at each position
      for(int at=0; at<valid; ++at)
      {
        for(int s=0; s<2; ++s)
        {
          a[at] = s ? 5.0f : -5.0f;
          runCase(len, 0, data, at);
        }
        a[at] = aBuf[at];
      }
    }
  }

  const bool pass = (nbad == 0);
  if(!pass) ++gFailed;

  printf("  %-34s %s  [%ld values, exact", name, pass ? "PASS" : "FAIL",
         compared);
  if(nbad) printf(", %d wrong", nbad);
  printf("]\n");
  if(nbad > nshown)
    printf("      ... %d further problems not shown\n", nbad - nshown);
}

//==============================================================================
// Integer to float conversion with min/max test driver
//
//   runEvm(in, len, out, minOut, maxOut, noInit)  calls the EVM function
//
// The converted values, minOut and maxOut must all equal the baseline exactly
// (static_cast<float> of the input).  As in checkMinMax(), every length is run
// with noInit false, and with noInit true from starting values inside, outside
// and the wrong way round relative to the data range, on the data as given,
// made all positive and made all negative, and then with a single minimum,
// and then a single maximum, at every position (on the data halved, so that
// the type's lowest and highest values are unique).  Past the end the input
// alternates between the type's lowest and highest values, so a read past the
// end that took part would change the result, and the output is checked for
// writes past the end.  Each length is also run once in place, with the
// integers at the front of the float output buffer.
//==============================================================================
template<class T, class FEVM>
static void checkCnvrtMinMax(const char* name, const T* const aBuf,
                             FEVM runEvm)
{
  ++gTests;
  long compared = 0;
  int  nbad = 0, nshown = 0, nstray = 0;
  const T lo = std::numeric_limits<T>::lowest();
  const T hi = std::numeric_limits<T>::max();
  std::vector<T>     a(MAXLEN + GUARD);
  std::vector<float> out(MAXLEN + GUARD);

  // starting values when noInit is true, case 0 runs with noInit false
  const float initMin[] = {0.0f, -0.5f, -3.0e9f,  3.0e9f};
  const float initMax[] = {0.0f,  0.5f,  3.0e9f, -3.0e9f};
  const char* const dataName[] = {"as given", "all positive", "all negative",
                                  "single min/max"};

  // runs the conversion of a with noInit from case init and checks the result
  auto runCase = [&](const int len, const int init, const int data,
                     const int at, const bool inPlace)
  {
    const int  valid  = (len > 0) ? len : 0;
    const bool noInit = (init != 0);
    // without noInit the outputs start as the sentinel, which a negative or
    // zero length must leave in place
    float mn = noInit ? initMin[init] : sentinel<float>();
    float mx = noInit ? initMax[init] : sentinel<float>();
    float wantMin = mn;
    float wantMax = mx;
    if(valid > 0)
    {
      if(!noInit)
      {
        wantMin = static_cast<float>(a[0]);
        wantMax = static_cast<float>(a[0]);
      }
      for(int k=0; k<valid; ++k)
      {
        wantMin = std::min(wantMin, static_cast<float>(a[k]));
        wantMax = std::max(wantMax, static_cast<float>(a[k]));
      }
    }
    for(int k=0; k<MAXLEN+GUARD; ++k) out[k] = sentinel<float>();
    if(inPlace)
    {
      if(valid > 0) std::memcpy(out.data(), a.data(), valid*sizeof(T));
      runEvm(reinterpret_cast<const T*>(out.data()), len, out.data(), mn, mx,
             noInit);
    }
    else
    {
      runEvm(a.data(), len, out.data(), mn, mx, noInit);
    }
    ++compared;

    const char* where = inPlace ? ", in place" : "";
    for(int k=0; k<MAXLEN+GUARD; ++k)
    {
      if((k < valid) ? (out[k] == static_cast<float>(a[k]))
                     : isSentinel(out[k])) continue;
      if(k < valid) ++nbad;
      else          ++nstray;
      if(nshown < MAX_REPORT)
      {
        ++nshown;
        if(k < valid)
          printf("      len=%-5d element %-5d evm=%.9g baseline=%.9g "
                 "[%s, noInit=%d%s]\n", len, k, out[k],
                 static_cast<float>(a[k]), dataName[data], noInit ? 1 : 0,
                 where);
        else
          printf("      len=%-5d element %-5d WROTE PAST END evm=%.9g "
                 "[%s, noInit=%d%s]\n", len, k, out[k], dataName[data],
                 noInit ? 1 : 0, where);
      }
    }

    const bool ok = (valid == 0 && !noInit) ?
        (isSentinel(mn) && isSentinel(mx)) :
        (mn == wantMin && mx == wantMax);
    if(ok) return;
    ++nbad;
    if(nshown < MAX_REPORT)
    {
      ++nshown;
      printf("      len=%-5d evm min,max=%.9g,%.9g baseline=%.9g,%.9g "
             "[%s", len, mn, mx, wantMin, wantMax, dataName[data]);
      if(at >= 0) printf(" at element %d", at);
      printf(", noInit=%d%s]\n", noInit ? 1 : 0, where);
    }
  };

  for(size_t li=0; li<gLengths.size(); ++li)
  {
    const int len   = gLengths[li];
    const int valid = (len > 0) ? len : 0;
    for(int data=0; data<4; ++data)
    {
      for(int k=0; k<MAXLEN+GUARD; ++k)
      {
        // |x|/2 + 1 can't overflow, even for the lowest value
        const T half = static_cast<T>(std::abs(static_cast<int64_t>(
            aBuf[(k < valid) ? k : 0]))/2 + 1);
        if(k >= valid)     a[k] = (k&1) ? hi : lo;
        else if(data == 1) a[k] = half;
        else if(data == 2) a[k] = static_cast<T>(-half);
        else if(data == 3) a[k] = static_cast<T>(aBuf[k]/2);
        else               a[k] = aBuf[k];
      }
      if(data < 3)
      {
        for(int init=0; init<4; ++init) runCase(len, init, data, -1, false);
        if(data == 0) runCase(len, 0, data, -1, true);
        continue;
      }
      // a minimum then a maximum outside the data range at each position
      for(int at=0; at<valid; ++at)
      {
        for(int s=0; s<2; ++s)
        {
          a[at] = s ? hi : lo;
          runCase(len, 0, data, at, false);
        }
        a[at] = static_cast<T>(aBuf[at]/2);
      }
    }
  }

  const bool pass = (nbad == 0 && nstray == 0);
  if(!pass) ++gFailed;

  printf("  %-34s %s  [%ld values, exact", name, pass ? "PASS" : "FAIL",
         compared);
  if(nbad)   printf(", %d wrong", nbad);
  if(nstray) printf(", %d past end", nstray);
  printf("]\n");
  if(nbad + nstray > nshown)
    printf("      ... %d further problems not shown\n",
           nbad + nstray - nshown);
}

//==============================================================================
// Special values test driver (polar/cartesian conversions)
//
//   runEvm(in, len, out)  calls the EVM function under test
//   want(v, part)         the expected value of part 0 (real) or 1 (imag) of
//                         the output for the input v, in double
//   size(v, part)         what the error of that part is measured against
//
// Every combination of a real part from re[] and an imaginary part from im[]
// (zeros of both signs, +-1, infinities, NaN, denormals, large values) goes
// through the function, out of place and then in place; the inputs fill
// several registers and a partial one.  A NaN must give NaN and an infinity must match exactly.  If zeroSign
// is set, a zero angle (part 1) must also keep its sign, since atan2(-0, 1)
// is -0.  Otherwise the error must be within TOL_TRIG of size, plus FLT_MIN
// since a denormal result may be flushed to zero (icpx turns on FTZ/DAZ by
// default; under DAZ the double references see denormal inputs as zero too,
// as the kernels do).  The output buffer is also checked for writes past the end.
//==============================================================================
template<class FEVM, class FWANT, class FSIZE>
static void checkPolarSpecial(const char* name, const float* const re,
                              const int nre, const float* const im,
                              const int nim, FEVM runEvm, FWANT want,
                              FSIZE size, const bool zeroSign)
{
  ++gTests;
  long compared = 0;
  int  nbad = 0, nshown = 0;
  const int n = nre*nim;
  std::vector<complex<float> > in(n), out(n + GUARD);
  for(int a=0; a<nre; ++a)
    for(int b=0; b<nim; ++b) in[a*nim + b] = complex<float>(re[a], im[b]);

  for(int mode=0; mode<2; ++mode)
  {
    for(int k=0; k<n+GUARD; ++k)
      out[k] = (mode == 1 && k < n) ? in[k] : sentinel<complex<float> >();
    if(mode == 0) runEvm(in.data(), n, out.data());
    else          runEvm(out.data(), n, out.data());
    for(int k=0; k<n+GUARD; ++k)
    {
      bool ok[2] = {true, true};
      if(k >= n) ok[0] = isSentinel(out[k]);
      else
      {
        const double got[2] = {out[k].real(), out[k].imag()};
        for(int part=0; part<2; ++part)
        {
          const double w = want(in[k], part);
          const double g = got[part];
          ++compared;
          if(std::isnan(w))      ok[part] = std::isnan(g);
          else if(std::isinf(w)) ok[part] = (g == w);
          else if(zeroSign && part == 1 && w == 0)
            ok[part] = (g == 0 && std::signbit(g) == std::signbit(w));
          else ok[part] = (std::fabs(g - w) <= TOL_TRIG*size(in[k], part) +
                                       std::numeric_limits<float>::min());
        }
      }
      if(ok[0] && ok[1]) continue;
      ++nbad;
      if(nshown < MAX_REPORT)
      {
        ++nshown;
        if(k >= n)
          printf("      element %-5d WROTE PAST END%s\n", k,
                 mode ? " [in place]" : "");
        else
          printf("      in=(%g,%g) evm=(%.9g,%.9g) want=(%.9g,%.9g)%s\n",
                 in[k].real(), in[k].imag(), out[k].real(), out[k].imag(),
                 want(in[k], 0), want(in[k], 1), mode ? " [in place]" : "");
      }
    }
  }

  const bool pass = (nbad == 0);
  if(!pass) ++gFailed;
  printf("  %-34s %s  [%ld values, NaN/inf exact", name, pass ? "PASS" : "FAIL",
         compared);
  if(nbad) printf(", %d wrong", nbad);
  printf("]\n");
  if(nbad > nshown)
    printf("      ... %d further problems not shown\n", nbad - nshown);
}


//==============================================================================
// Test driver for the polar to cartesian conversions with the magnitude and
// the phase in separate real buffers (complex output)
//
//   runEvm(mag, ph, len, out)  calls the EVM function under test
//   baseline(k)                returns the expected value of element k
//
// Mode 0 uses separate buffers; modes 1 and 2 put the complex output on top
// of mag or of the phases (the real input in the front of the complex
// output, which the header allows).  The error of each element is measured
// against |mag[k]|, since one part can be near zero, at every length, and
// nothing past out[len-1] may be written.  With flush to zero on, a result
// below FLT_MIN may come out as 0.
//==============================================================================
template<class FEVM, class FREF>
static void checkPolarSplit(const char* name, const float* const magBuf,
                            const float* const phBuf, FEVM runEvm,
                            FREF baseline)
{
  ++gTests;
  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0;
  double worstRel = 0.0;
  const bool ftz = (_mm_getcsr() & _MM_FLUSH_ZERO_MASK) == _MM_FLUSH_ZERO_ON;
  std::vector<complex<float> > out(MAXLEN + GUARD);
  float* const of = reinterpret_cast<float*>(out.data());
  for(int mode=0; mode<3; ++mode)
  {
    const char* modeName = (mode==0) ? "" :
                           (mode==1) ? " [out on mag]" : " [out on phase]";
    for(size_t li=0; li<gLengths.size(); ++li)
    {
      const int len   = gLengths[li];
      const int valid = (len > 0) ? len : 0;
      for(int k=0; k<MAXLEN+GUARD; ++k) out[k] = sentinel<complex<float> >();
      if(mode == 1)
      {
        for(int k=0; k<valid; ++k) of[k] = magBuf[k];
        runEvm(of, phBuf, len, out.data());
      }
      else if(mode == 2)
      {
        for(int k=0; k<valid; ++k) of[k] = phBuf[k];
        runEvm(magBuf, of, len, out.data());
      }
      else runEvm(magBuf, phBuf, len, out.data());
      for(int k=0; k<valid; ++k)
      {
        const cdbl   got   = toCdbl(out[k]);
        const cdbl   want  = baseline(k);
        const double adiff = std::abs(got - want);
        const double size  = std::max(
            std::fabs(static_cast<double>(magBuf[k])),
            static_cast<double>(std::numeric_limits<float>::min()));
        const double rdiff = adiff/size;
        ++compared;
        if(rdiff <= TOL_TRIG)
        {
          if(rdiff > worstRel) worstRel = rdiff;
          continue;
        }
        if(ftz && got == cdbl(0.0, 0.0) &&
           std::abs(want) < std::numeric_limits<float>::min()) continue;
        ++nbad;
        if(nshown < MAX_REPORT)
        {
          ++nshown;
          printf("      len=%-5d element %-5d evm=%s baseline=%s "
                 "reldiff=%.6g%s\n", len, k, format(got, true).c_str(),
                 format(want, true).c_str(), rdiff, modeName);
        }
      }
      for(int k=valid; k<MAXLEN+GUARD; ++k)
      {
        if(isSentinel(out[k])) continue;
        ++nstray;
        if(nshown < MAX_REPORT)
        {
          ++nshown;
          printf("      len=%-5d element %-5d WROTE PAST END evm=%s%s\n", len,
                 k, format(toCdbl(out[k]), true).c_str(), modeName);
        }
      }
    }
  }
  const bool pass = (nbad == 0 && nstray == 0);
  if(!pass) ++gFailed;
  printf("  %-34s %s  [%ld values, tol %.0e, worst rel %.3e", name,
         pass ? "PASS" : "FAIL", compared, TOL_TRIG, worstRel);
  if(nbad)   printf(", %d wrong", nbad);
  if(nstray) printf(", %d past end", nstray);
  printf("]\n");
  if(nbad + nstray > nshown)
    printf("      ... %d further problems not shown\n",
           nbad + nstray - nshown);
}

//==============================================================================
// Special values for the inversion functions
//
//   runEvm(in, len, out, useEps, scale)  calls the function under test, with
//                                        eps = GEPS if useEps (the square
//                                        roots without it use eps = 0) and
//                                        scale for the scaled ones
//
// Zeros of both signs, infinities, NaN, +-1, the largest and smallest normal
// floats and a few large and small values go through every lane position,
// including the partial last register, at lengths 1 to 40, with and without
// eps, and for the scaled functions with scale = gSf and -0.  The expected
// results follow the rules of 1/x and 1/sqrt(x) rather than arithmetic on
// infinities: NaN or a negative square root gives NaN, a zero gives an
// infinity and an infinity a zero, of the right signs, and a zero scale times
// an infinity gives NaN.  These must match exactly (bit tests, since icpx's
// fast math may drop NaN checks).  A result below FLT_MIN may come out as a
// zero of the same sign (rcp_ps always flushes it, and so does FTZ);
// otherwise the error must be within tol.  Denormal inputs are included for
// the true versions, and for the approximations only under DAZ, where every
// path sees them as zero: without it, AVX/AVX2 rcp_ps and rsqrt_ps still
// treat them as zero (giving INF) while AVX-512 and the defaults do not.
//==============================================================================
static uint32_t floatBits(const float v)
{
  uint32_t b;
  std::memcpy(&b, &v, sizeof(b));
  return b;
}
static bool bitsNan(const float v)
{
  return (floatBits(v) & 0x7fffffffu) > 0x7f800000u;
}
static bool bitsInf(const float v)
{
  return (floatBits(v) & 0x7fffffffu) == 0x7f800000u;
}
static bool bitsZero(const float v)
{
  return (floatBits(v) & 0x7fffffffu) == 0;
}
static bool bitsNeg(const float v)
{
  return (floatBits(v) >> 31) != 0;
}
static bool bitsDenormal(const float v)
{
  const uint32_t m = floatBits(v) & 0x7fffffffu;
  return m != 0 && m < 0x00800000u;
}

// expected result: kind 0 = finite value, 1 = NaN, 2 = infinity, 3 = zero
struct InvWant { int kind; bool neg; double val; };

static InvWant invWant(float x, const bool sqrtForm, const bool useEps,
                       const bool scaled, const float scale, const bool daz)
{
  // std::max(x, eps), which keeps a NaN x
  if(useEps || sqrtForm)
  {
    const float e = useEps ? GEPS : 0.0f;
    if(!bitsNan(x) && x < e) x = e;
  }
  if(daz && bitsDenormal(x)) x = bitsNeg(x) ? -0.0f : 0.0f;
  InvWant w = {0, false, 0.0};
  if(bitsNan(x) || (sqrtForm && bitsNeg(x) && !bitsZero(x)))
  {
    w.kind = 1;
    return w;
  }
  if(bitsZero(x))     { w.kind = 2; w.neg = bitsNeg(x); }
  else if(bitsInf(x)) { w.kind = 3; w.neg = bitsNeg(x); }
  else
  {
    w.val = sqrtForm ? 1.0/std::sqrt(static_cast<double>(x))
                     : 1.0/static_cast<double>(x);
  }
  if(scaled)
  {
    if(w.kind == 2 && bitsZero(scale)) w.kind = 1;
    else if(w.kind != 0) w.neg = (w.neg != bitsNeg(scale));
    else if(bitsZero(scale))
    {
      w.kind = 3;
      w.neg = ((w.val < 0) != bitsNeg(scale));
    }
    else w.val *= scale;
  }
  // a finite result too large for a float (1/x of a small denormal)
  if(w.kind == 0 && std::fabs(w.val) > std::numeric_limits<float>::max())
  {
    w.kind = 2;
    w.neg = w.val < 0;
  }
  return w;
}

template<class FEVM>
static void checkInvSpecial(const char* name, const double tol,
                            const bool approx, const bool sqrtForm,
                            const bool scaled, FEVM runEvm)
{
  ++gTests;
  const bool daz = (_mm_getcsr() & _MM_DENORMALS_ZERO_MASK) ==
                   _MM_DENORMALS_ZERO_ON;
  const float fmax = std::numeric_limits<float>::max();
  const float fmin = std::numeric_limits<float>::min();
  std::vector<float> vals = {0.0f, -0.0f,
      std::numeric_limits<float>::infinity(),
      -std::numeric_limits<float>::infinity(),
      std::numeric_limits<float>::quiet_NaN(), 1.0f, -1.0f, 0.5f, -2.0f,
      fmax, -fmax, fmin, -fmin, 1.0e30f, -3.0e-30f, 3.0f};
  if(!approx || daz)
  {
    vals.push_back(std::numeric_limits<float>::denorm_min());
    vals.push_back(-3.0e-39f);
    vals.push_back(1.0e-38f);
  }
  const int nv = static_cast<int>(vals.size());
  const float scales[2] = {gSf, -0.0f};
  long compared = 0;
  int  nbad = 0, nshown = 0;
  std::vector<float> in(64), out(64);
  for(int useEps=0; useEps<2; ++useEps)
  for(int si=0; si<(scaled ? 2 : 1); ++si)
  for(int len=1; len<=40; ++len)
  {
    for(int i=0; i<len; ++i) in[i] = vals[(i + len)%nv];
    runEvm(in.data(), len, out.data(), useEps != 0, scales[si]);
    for(int i=0; i<len; ++i)
    {
      const InvWant w = invWant(in[i], sqrtForm, useEps != 0, scaled,
                                scales[si], daz);
      const float g = out[i];
      bool ok;
      ++compared;
      if(w.kind == 1)      ok = bitsNan(g);
      else if(w.kind == 2) ok = bitsInf(g) && bitsNeg(g) == w.neg;
      else if(w.kind == 3) ok = bitsZero(g) && bitsNeg(g) == w.neg;
      else if(std::fabs(w.val) < fmin && bitsZero(g))
        ok = (bitsNeg(g) == (w.val < 0));
      else
        ok = !bitsNan(g) && !bitsInf(g) &&
             std::fabs(g - w.val)/std::max(std::fabs(w.val),
                                           static_cast<double>(fmin)) <= tol;
      if(ok) continue;
      ++nbad;
      if(nshown < MAX_REPORT)
      {
        ++nshown;
        char want[32];
        if(w.kind == 0) snprintf(want, sizeof(want), "%.9g", w.val);
        else snprintf(want, sizeof(want), "%s%s", w.kind == 1 ? "" :
                      (w.neg ? "-" : "+"), w.kind == 1 ? "NaN" :
                      (w.kind == 2 ? "INF" : "0"));
        printf("      len=%-3d in=%-13g%s%s evm=%-13g want=%s\n", len, in[i],
               useEps ? " eps" : "", (scaled && si) ? " scale=-0" : "", g,
               want);
      }
    }
  }
  const bool pass = (nbad == 0);
  if(!pass) ++gFailed;
  printf("  %-34s %s  [%ld values, NaN/INF/0 exact", name,
         pass ? "PASS" : "FAIL", compared);
  if(nbad) printf(", %d wrong", nbad);
  printf("]\n");
  if(nbad > nshown)
    printf("      ... %d further problems not shown\n", nbad - nshown);
}

//==============================================================================
// Special values for exp and log
//
//   runEvm(in, len, out)  calls the function under test
//
// Zeros of both signs, infinities, NaN, +-1, the largest and smallest normal
// floats, denormals, large values and values around exp's overflow (88.72)
// and underflow (denormal results from -87.34, 0 below -103.97) go through
// every lane position, including the partial last register, at lengths 1 to
// 40, out of place and in place.  NaN, infinities and zeros must match
// exactly (bit tests, since icpx's fast math may drop NaN checks): exp gives
// +INF for +INF and on overflow and +0 for -INF and below -103.97; log gives
// -INF for +-0, NaN for a negative number or -INF, +INF for +INF and +0 for
// 1.  Under DAZ a denormal input acts as a zero of its sign (exp gives 1, log
// -INF), and under FTZ a result below FLT_MIN may come out as +0; otherwise
// the error must be within tol of max(|result|, FLT_MIN).  negDenNaN also
// accepts NaN from log for a negative denormal under DAZ (icpx's default
// version, see the call).
//==============================================================================
template<class FEVM>
static void checkExpLogSpecial(const char* name, const bool isLog,
                               const bool negDenNaN, FEVM runEvm)
{
  ++gTests;
  const unsigned int csr = _mm_getcsr();
  const bool daz = (csr & _MM_DENORMALS_ZERO_MASK) == _MM_DENORMALS_ZERO_ON;
  const bool ftz = (csr & _MM_FLUSH_ZERO_MASK) == _MM_FLUSH_ZERO_ON;
  const float inf  = std::numeric_limits<float>::infinity();
  const float fmax = std::numeric_limits<float>::max();
  const float fmin = std::numeric_limits<float>::min();
  // e^x rounds to INF from FLT_MAX plus half an ulp, to 0 below 2^-150
  const double expInf = fmax + std::ldexp(1.0, 103);
  const double expZero = std::ldexp(1.0, -150);
  const std::vector<float> vals = {0.0f, -0.0f, inf, -inf,
      std::numeric_limits<float>::quiet_NaN(), 1.0f, -1.0f, 0.5f, 2.0f,
      fmax, -fmax, fmin, -fmin, std::numeric_limits<float>::denorm_min(),
      -3.0e-39f, 1.0e-38f, 1.0e30f, -1.0e30f, 88.72283f, 88.7229f, 89.0f,
      100.0f, -87.33654f, -87.5f, -103.9f, -103.98f, -104.0f, -200.0f,
      1.0e-30f, 1.0000001f, 0.99999994f};
  const int nv = static_cast<int>(vals.size());
  long compared = 0;
  int  nbad = 0, nshown = 0;
  std::vector<float> in(64), out(64);
  for(int mode=0; mode<2; ++mode)
  for(int len=1; len<=40; ++len)
  {
    for(int i=0; i<len; ++i) in[i] = vals[(i + len)%nv];
    if(mode == 0) runEvm(in.data(), len, out.data());
    else
    {
      for(int i=0; i<len; ++i) out[i] = in[i];
      runEvm(out.data(), len, out.data());
    }
    for(int i=0; i<len; ++i)
    {
      float x = in[i];
      if(daz && bitsDenormal(x)) x = bitsNeg(x) ? -0.0f : 0.0f;
      // kind 0 = finite value, 1 = NaN, 2 = infinity, 3 = zero (as invWant)
      InvWant w = {0, false, 0.0};
      if(bitsNan(x)) w.kind = 1;
      else if(isLog)
      {
        if(bitsNeg(x) && !bitsZero(x)) w.kind = 1;
        else if(bitsZero(x))           { w.kind = 2; w.neg = true; }
        else if(bitsInf(x))            w.kind = 2;
        else if(x == 1.0f)             w.kind = 3;
        else w.val = std::log(static_cast<double>(x));
      }
      else
      {
        w.val = bitsInf(x) ? (bitsNeg(x) ? 0.0 : HUGE_VAL) :
                std::exp(static_cast<double>(x));
        if(w.val >= expInf)      w.kind = 2;
        else if(w.val < expZero) w.kind = 3;
      }
      const float g = out[i];
      bool ok;
      ++compared;
      if(w.kind == 1 || (negDenNaN && isLog && daz && bitsDenormal(in[i]) &&
                         bitsNeg(in[i]) && bitsNan(g))) ok = bitsNan(g);
      else if(w.kind == 2) ok = bitsInf(g) && bitsNeg(g) == w.neg;
      else if(w.kind == 3) ok = bitsZero(g) && !bitsNeg(g);
      else if(ftz && std::fabs(w.val) < fmin && bitsZero(g)) ok = true;
      else
        ok = !bitsNan(g) && !bitsInf(g) &&
             std::fabs(g - w.val)/std::max(std::fabs(w.val),
                 static_cast<double>(fmin)) <= TOL_EXPLOG;
      if(ok) continue;
      ++nbad;
      if(nshown < MAX_REPORT)
      {
        ++nshown;
        char want[32];
        if(w.kind == 0) snprintf(want, sizeof(want), "%.9g", w.val);
        else snprintf(want, sizeof(want), "%s%s", w.kind == 1 ? "" :
                      (w.neg ? "-" : "+"), w.kind == 1 ? "NaN" :
                      (w.kind == 2 ? "INF" : "0"));
        printf("      len=%-3d in=%-13g evm=%-13g want=%s%s\n", len, in[i], g,
               want, mode ? " [in place]" : "");
      }
    }
  }
  const bool pass = (nbad == 0);
  if(!pass) ++gFailed;
  printf("  %-34s %s  [%ld values, NaN/INF/0 exact", name,
         pass ? "PASS" : "FAIL", compared);
  if(nbad) printf(", %d wrong", nbad);
  printf("]\n");
  if(nbad > nshown)
    printf("      ... %d further problems not shown\n", nbad - nshown);
}

//==============================================================================
// NCO tuning (tune, tunecx, tunesp)
//
//   runEvm(in, len, phInit, phIncr, out, phOut)  calls the function under test
//       (the overloads without phOut ignore it) and returns its result
//
// in is gAc or gAf, and element k must come out as
// in[k]*e^(j*scale*(phInit + k*phIncr)), scale = 1 for radians and 2*pi for
// cycles, within TOL_TRIG + m*drift of |in[k]|: tune and tunecx advance their
// phasor with complex multiplies, so their error may ramp up along the buffer
// (drift), m elements from where it starts: element 0 for complex input, the
// last element for real input, which runs from the end to the start.
// tunesp's error may not ramp up (drift = 0).  With ph, phOut[k] must be the
// phasor within the same bound.  Each pair of phases runs every length with
// separate buffers and in place: the complex input under the output, or the
// real input in the front of the (wider) complex output buffer, which the
// real versions handle by running from the end to the start.  Nothing may be
// written past out[len-1] or phOut[len-1] (len <= 0 writes nothing), and the
// result must be the phase for the next call, phInit + len*phIncr with whole
// cycles removed (fmod(., 2*pi) for tune), or phInit itself for len <= 0.
// A NaN or infinite phInit must give NaN in every element and phasor, and so
// must a NaN or infinite phIncr, except in element 0 with k0Finite (the
// complex tune and tunecx, whose first phasor comes from phInit alone).
//==============================================================================
template<class TI, class FEVM>
static void checkTune(const char* name, const bool cyc, const bool ph,
                      const bool k0Finite, const double drift,
                      const std::vector<double>& phases, FEVM runEvm)
{
  ++gTests;
  const double scale = cyc ? 2.0*M_PI : 1.0;
  const TI* src;
  if constexpr (std::is_same<TI, float>::value) src = gAf.data();
  else src = gAc.data();
  long compared = 0;
  int  nbad = 0, nshown = 0;
  double worst = 0.0;  // worst error as a fraction of its bound
  std::vector<complex<float> > out(MAXLEN + GUARD), pho(MAXLEN + GUARD);
  const float* of = reinterpret_cast<const float*>(out.data());
  const float* pf = reinterpret_cast<const float*>(pho.data());
  auto report = [&](const int len, const int k, const char* what,
                    const cdbl got, const cdbl want, const char* mode)
  {
    ++nbad;
    if(nshown >= MAX_REPORT) return;
    ++nshown;
    printf("      len=%-5d element %-5d %s evm=%s want=%s%s\n", len, k, what,
           format(got, true).c_str(), format(want, true).c_str(), mode);
  };
  for(size_t pi=0; pi+1<phases.size(); pi+=2)
  {
    const double p0 = phases[pi], dp = phases[pi+1];
    for(int mode=0; mode<2; ++mode)
    for(size_t li=0; li<gLengths.size(); ++li)
    {
      const int len = gLengths[li];
      const int valid = (len > 0) ? len : 0;
      const char* mname = mode ? " [in place]" : "";
      for(int k=0; k<MAXLEN+GUARD; ++k)
        out[k] = pho[k] = sentinel<complex<float> >();
      double ret;
      if(mode == 0) ret = runEvm(src, len, p0, dp, out.data(), pho.data());
      else
      {
        // the input in the front of the output buffer
        std::memcpy(out.data(), src, valid*sizeof(TI));
        ret = runEvm(reinterpret_cast<const TI*>(out.data()), len, p0, dp,
                     out.data(), pho.data());
      }
      for(int k=0; k<valid; ++k)
      {
        const cdbl w = std::polar(1.0, scale*(p0 + k*dp));
        const cdbl x = toCdbl(src[k]);
        const int m = std::is_same<TI, float>::value ? len-1-k : k;
        const double bound = TOL_TRIG + m*drift;
        const double e = std::abs(toCdbl(out[k]) - x*w)/std::max(std::abs(x),
            static_cast<double>(std::numeric_limits<float>::min()));
        ++compared;
        if(e <= bound) worst = std::max(worst, e/bound);
        else report(len, k, isSentinel(out[k]) ? "NOT WRITTEN" : "",
                    toCdbl(out[k]), x*w, mname);
        if(!ph) continue;
        const double ep = std::abs(toCdbl(pho[k]) - w);
        if(ep <= bound) worst = std::max(worst, ep/bound);
        else report(len, k, "phasor", toCdbl(pho[k]), w, mname);
      }
      for(int k=2*valid; k<2*(MAXLEN+GUARD); ++k)
      {
        if(!isSentinel(of[k]))
          report(len, k/2, "WROTE PAST END", toCdbl(out[k/2]), cdbl(), mname);
        if(ph && !isSentinel(pf[k]))
          report(len, k/2, "phasor WROTE PAST END", toCdbl(pho[k/2]), cdbl(),
                 mname);
      }
      if(len <= 0)
      {
        if(!(ret == p0))
          report(len, 0, "phase changed", cdbl(ret), cdbl(p0), mname);
      }
      else
      {
        double want = p0 + len*dp;
        want = cyc ? want - std::trunc(want) : std::fmod(want, 2.0*M_PI);
        if(!(std::fabs(ret - want) <= 1.0e-9))
          report(len, len, "next phase", cdbl(ret), cdbl(want), mname);
      }
    }
    // NaN and infinite phases
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double bad[4][2] = {{nan, dp}, {HUGE_VAL, dp}, {p0, nan},
                              {p0, -HUGE_VAL}};
    for(int b=0; b<4; ++b)
    for(int len=1; len<=40; ++len)
    {
      for(int k=0; k<MAXLEN+GUARD; ++k)
        out[k] = pho[k] = sentinel<complex<float> >();
      runEvm(src, len, bad[b][0], bad[b][1], out.data(), pho.data());
      for(int k=0; k<len; ++k)
      {
        ++compared;
        if(k == 0 && k0Finite && b >= 2)
        {
          const cdbl w = std::polar(1.0, scale*p0);
          const double e = std::abs(toCdbl(out[0]) - toCdbl(src[0])*w)/
              std::max(std::abs(toCdbl(src[0])), 1.0e-30);
          if(!(e <= TOL_TRIG))
            report(len, 0, "phIncr NaN/INF", toCdbl(out[0]),
                   toCdbl(src[0])*w, "");
          continue;
        }
        if(!bitsNan(out[k].real()) || !bitsNan(out[k].imag()))
          report(len, k, b < 2 ? "phInit NaN/INF" : "phIncr NaN/INF",
                 toCdbl(out[k]), cdbl(nan, nan), "");
        if(ph && (!bitsNan(pho[k].real()) || !bitsNan(pho[k].imag())))
          report(len, k, b < 2 ? "phasor, phInit NaN/INF" :
                 "phasor, phIncr NaN/INF", toCdbl(pho[k]), cdbl(nan, nan), "");
      }
    }
  }
  const bool pass = (nbad == 0);
  if(!pass) ++gFailed;
  printf("  %-34s %s  [%ld values, worst %.2f of bound", name,
         pass ? "PASS" : "FAIL", compared, worst);
  if(nbad) printf(", %d wrong", nbad);
  printf("]\n");
  if(nbad > nshown)
    printf("      ... %d further problems not shown\n", nbad - nshown);
}

//==============================================================================
// Complex exponentials (cxexp, cxexpsi, cxexpso, cxexpss)
//
//   runEvm(in, len, out)  calls the function under test
//
// Element k must come out as os*e^(j*(is*in[k] + eps)) for the function's
// input scale is, output scale os and eps.  The phase is computed in float,
// so its rounding (about 2^-24 of |is*in[k]| for the product and of the
// phase for the sum, or once for an FMA) is part of the bound, which is
// TOL_TRIG plus that, of |os|.  Every length runs with separate buffers and
// in place (the real input in the front of the wider complex output buffer,
// which the functions handle by running from the end to the start), and
// nothing may be written past out[len-1] (len <= 0 writes nothing).
//
// The special values (zeros of both signs, +-1, infinities, NaN, denormals,
// large values, in every lane at lengths 1 to 40, out of place and in place)
// expect the phase in float, ph = is*in[k] + eps.  An infinite or NaN ph
// must give NaN in both parts, and ph = +-0 exactly (os, os*(+-0)) with the
// sign of the zero; the overloads without eps pass eps = -0, which keeps the
// sign of a zero phase.  evm_trig.h's sincosRad gives sin(-0) = +0 (keeping
// the sign would cost it 2.5-8%), so on the SIMD paths without SVML, and on
// the AVX path, which always uses it, the zero may have either sign.  The
// huge values only go through the functions whose phase is exact (no scale,
// or eps = -0), since a rounding of the phase of 1 ulp at 1e30 changes the
// result completely.
//==============================================================================
template<class FEVM>
static void checkCxexp(const char* name, const std::vector<float>& src,
                       const double is, const double os, const float eps,
                       const bool huge, FEVM runEvm)
{
  ++gTests;
  long compared = 0;
  int  nbad = 0, nshown = 0;
  double worst = 0.0;
  const int n0 = static_cast<int>(src.size());
  std::vector<complex<float> > out(n0 + GUARD);
  const std::string lvl = EVM::getOptLevel();
#if defined(EVM_SVML)
  const bool anyZeroSign = (lvl == "AVX"); // AVX always uses evm_trig.h
#else
  const bool anyZeroSign = (lvl == "AVX512" || lvl == "AVX2" || lvl == "AVX");
#endif
  const float* of = reinterpret_cast<const float*>(out.data());
  const float* cur = src.data(); // the input of the current run
  auto report = [&](const int len, const int k, const char* what,
                    const cdbl got, const cdbl want, const char* mode)
  {
    ++nbad;
    if(nshown >= MAX_REPORT) return;
    ++nshown;
    printf("      len=%-5d element %-5d %s in=%.9g evm=%s want=%s%s\n", len, k,
           what, k < len ? cur[k] : 0.0f, format(got, true).c_str(),
           format(want, true).c_str(), mode);
  };
  auto run = [&](const float* in, const int len, const int mode,
                 std::vector<complex<float> >& o)
  {
    for(size_t k=0; k<o.size(); ++k) o[k] = sentinel<complex<float> >();
    const int valid = (len > 0) ? len : 0;
    if(mode == 0) runEvm(in, len, o.data());
    else
    {
      std::memcpy(o.data(), in, valid*sizeof(float));
      runEvm(reinterpret_cast<const float*>(o.data()), len, o.data());
    }
  };
  // ordinary values, every length
  for(int mode=0; mode<2; ++mode)
  for(size_t li=0; li<gLengths.size(); ++li)
  {
    const int len = std::min(gLengths[li], n0);
    const int valid = (len > 0) ? len : 0;
    const char* mname = mode ? " [in place]" : "";
    run(src.data(), len, mode, out);
    for(int k=0; k<valid; ++k)
    {
      const double x = src[k], p = is*x + eps;
      const double bound = TOL_TRIG + (std::fabs(is*x) + std::fabs(p))*
          std::ldexp(1.0, -24);
      const cdbl w = os*std::polar(1.0, p);
      const double e = std::abs(toCdbl(out[k]) - w)/std::fabs(os);
      ++compared;
      if(e <= bound) worst = std::max(worst, e/bound);
      else report(len, k, isSentinel(out[k]) ? "NOT WRITTEN" : "",
                  toCdbl(out[k]), w, mname);
    }
    for(int k=2*valid; k<2*(n0+GUARD); ++k)
      if(!isSentinel(of[k]))
        report(len, k/2, "WROTE PAST END", toCdbl(out[k/2]), cdbl(), mname);
  }
  // special values
  {
    const float inf = std::numeric_limits<float>::infinity();
    const float nan = std::numeric_limits<float>::quiet_NaN();
    std::vector<float> sp = {0.0f, -0.0f, 1.0f, -1.0f, 0.5f, inf, -inf, nan,
                             1.0e-40f, -1.0e-40f, 1.0e-30f};
    if(huge)
      for(float v : {3.4e38f, -3.4e38f, 1.0e30f, -1.0e7f, 65537.0f})
        sp.push_back(v);
    const int nsp = static_cast<int>(sp.size());
    std::vector<complex<float> > o(64 + GUARD);
    for(int mode=0; mode<2; ++mode)
    for(int len=1; len<=40; ++len)
    for(int start=0; start<nsp; ++start)
    {
      std::vector<float> in(len);
      for(int k=0; k<len; ++k) in[k] = sp[(start + k) % nsp];
      const char* mname = mode ? " [in place]" : "";
      cur = in.data();
      run(in.data(), len, mode, o);
      for(int k=0; k<len; ++k)
      {
        const float x = in[k];
        // the float phase.  Zero signs are worked out on the bits (bitsNeg),
        // since icpx's fast fp-model folds x + 0 to x and drops the sign of
        // zero constants here: the sum of two zeros is -0 only if both are -0
        const float px = static_cast<float>(is*x);
        const float pf = px + eps;
        const bool pzero = bitsZero(px) && bitsZero(eps);
        const bool pneg = bitsNeg(px) && bitsNeg(eps);
        const cdbl got = toCdbl(o[k]);
        ++compared;
        if(!std::isfinite(pf))
        {
          if(!bitsNan(o[k].real()) || !bitsNan(o[k].imag()))
            report(len, k, "NaN/INF phase", got, cdbl(nan, nan), mname);
          continue;
        }
        if(pzero)
        {
          // exactly (os, os*(+-0)): the zero's sign is the phase's, flipped
          // for a negative os
          const bool wneg = pneg != (os < 0);
          if(!(o[k].real() == static_cast<float>(os) && bitsZero(o[k].imag()) &&
               (bitsNeg(o[k].imag()) == wneg || anyZeroSign)))
            report(len, k, "zero phase", got, cdbl(os, wneg ? -0.0 : 0.0),
                   mname);
          continue;
        }
        const double p = static_cast<double>(pf);
        const double bound = TOL_TRIG + (std::fabs(is*x) + std::fabs(p))*
            std::ldexp(1.0, -24);
        const cdbl w = os*std::polar(1.0, p);
        const double e = std::abs(got - w)/std::fabs(os);
        if(!(e <= bound + std::numeric_limits<float>::min()))
          report(len, k, "special value", got, w, mname);
      }
      for(int k=2*len; k<2*(64+GUARD); ++k)
        if(!isSentinel(reinterpret_cast<const float*>(o.data())[k]))
          report(len, k/2, "WROTE PAST END", toCdbl(o[k/2]), cdbl(), mname);
    }
  }
  const bool pass = (nbad == 0);
  if(!pass) ++gFailed;
  printf("  %-34s %s  [%ld values, worst %.2f of bound", name,
         pass ? "PASS" : "FAIL", compared, worst);
  if(nbad) printf(", %d wrong", nbad);
  printf("]\n");
  if(nbad > nshown)
    printf("      ... %d further problems not shown\n", nbad - nshown);
}

//==============================================================================
// Special values for the arctangents (atan2, atan2s, multcatan2, multcatan2s)
//
//   runEvm(in1, in2, len, out)  calls the function under test (atan2 and
//                               atan2s ignore in2); TO is float for a real
//                               output and complex<float> for a complex one
//
// atan2 and atan2s run on every combination of a real and an imaginary part
// from zeros of both signs, +-1, infinities, NaN, denormals, the largest
// floats and 0.5.  multcatan2 and multcatan2s run on every pair of inputs
// whose parts are zeros of both signs, +-1, 0.5, -3, infinities and NaN:
// their products are exact, and finite ones neither overflow nor underflow,
// so the double reference sees the same product as the kernels.  The list
// runs from each of 17 starts to its end, so every value goes through every
// lane and the last register has every length, out of place and in place
// (on in1, and for multcatan2 on in2).  An angle of NaN must give NaN and a
// zero angle a zero of the same sign (times the scale's), checked on the
// bits since icpx's fast math may drop NaN checks; otherwise the error must
// be within TOL_TRIG of the angle, plus FLT_MIN since a denormal angle may
// be flushed to zero (under DAZ the double references see denormal inputs
// as zero too, as the kernels do).  A complex output's imaginary part must
// be +0.  With zeroAngle (multcatan2 and multcatan2s with the flag set), a
// product whose parts are both zero must give +0 instead.  Nothing may be
// written past out[len-1], and a real output on top of a complex input must
// leave the rest of the input as it was.
//==============================================================================
static bool dblNan(const double v)
{
  uint64_t b;
  std::memcpy(&b, &v, sizeof(b));
  return (b & 0x7fffffffffffffffull) > 0x7ff0000000000000ull;
}

template<class TO, class FEVM>
static void checkAtan2Special(const char* name, const bool two,
                              const double scale, const bool zeroAngle,
                              FEVM runEvm)
{
  ++gTests;
  long compared = 0;
  int  nbad = 0, nshown = 0;
  const bool cplx = IsComplex<TO>::value;
  const float inf = std::numeric_limits<float>::infinity();
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const std::vector<float> sp = two ?
      std::vector<float>{0.0f, -0.0f, 1.0f, -1.0f, 0.5f, -3.0f, inf, -inf,
                         nan} :
      std::vector<float>{0.0f, -0.0f, 1.0f, -1.0f, inf, -inf, nan, 1.0e-40f,
                         -1.0e-40f, 3.4e38f, -3.4e38f, 0.5f};
  std::vector<complex<float> > vals, a1, a2;
  for(float re : sp)
    for(float im : sp) vals.push_back(complex<float>(re, im));
  for(const complex<float> x : vals)
  {
    if(!two) { a1.push_back(x); a2.push_back(x); continue; }
    for(const complex<float> y : vals) { a1.push_back(x); a2.push_back(y); }
  }
  const int n = static_cast<int>(a1.size());
  const int nf = 2*(n + GUARD); // floats in each buffer
  // the angle in double, times the scale
  auto want = [&](const complex<float> x, const complex<float> y)
  {
    double re = x.real(), im = x.imag();
    if(two)
    {
      const double xr = x.real(), xi = x.imag();
      const double yr = y.real(), yi = y.imag();
      re = xr*yr + xi*yi;
      im = xi*yr - xr*yi;
    }
    return std::atan2(im, re)*scale;
  };
  // a zero product, which gives +0 with zeroAngle
  auto zeroProd = [&](const complex<float> x, const complex<float> y)
  {
    const double xr = x.real(), xi = x.imag();
    const double yr = y.real(), yi = y.imag();
    return xr*yr + xi*yi == 0.0 && xi*yr - xr*yi == 0.0;
  };
  std::vector<complex<float> > b1(n + GUARD), b2(n + GUARD), ob(n + GUARD);
  std::vector<complex<float> > orig(n + GUARD);
  for(int mode=0; mode<(two ? 3 : 2); ++mode)
  for(int s=0; s<=16; ++s)
  {
    const int len = n - s;
    const char* mname = (mode == 0) ? "" : (mode == 1) ? " [in place on in1]"
                                                       : " [in place on in2]";
    for(int k=0; k<n+GUARD; ++k)
    {
      b1[k] = (k < len) ? a1[s+k] : sentinel<complex<float> >();
      b2[k] = (k < len) ? a2[s+k] : sentinel<complex<float> >();
      ob[k] = sentinel<complex<float> >();
    }
    std::vector<complex<float> >& res = (mode == 0) ? ob : (mode == 1) ? b1
                                                                       : b2;
    orig = res;
    runEvm(b1.data(), b2.data(), len, reinterpret_cast<TO*>(res.data()));
    const float* rf = reinterpret_cast<const float*>(res.data());
    const float* of = reinterpret_cast<const float*>(orig.data());
    for(int k=0; k<len; ++k)
    {
      const double w = want(a1[s+k], a2[s+k]);
      const float g = cplx ? rf[2*k] : rf[k];
      bool ok;
      ++compared;
      if(zeroAngle && zeroProd(a1[s+k], a2[s+k])) ok = floatBits(g) == 0;
      else if(dblNan(w)) ok = bitsNan(g);
      else if(w == 0.0)
        ok = bitsZero(g) && (bitsNeg(g) == std::signbit(w));
      else ok = !bitsNan(g) && std::fabs(g - w) <= TOL_TRIG*std::fabs(w) +
                std::numeric_limits<float>::min();
      if(cplx && floatBits(rf[2*k+1]) != 0) ok = false;
      if(ok) continue;
      ++nbad;
      if(nshown >= MAX_REPORT) continue;
      ++nshown;
      printf("      len=%-5d element %-5d in1=(%g,%g)", len, k,
             a1[s+k].real(), a1[s+k].imag());
      if(two) printf(" in2=(%g,%g)", a2[s+k].real(), a2[s+k].imag());
      if(cplx) printf(" evm=(%.9g,%.9g)", rf[2*k], rf[2*k+1]);
      else     printf(" evm=%.9g", g);
      printf(" want=%.9g%s\n", w, mname);
    }
    // past the results: the rest of a complex input under a real output,
    // otherwise the sentinel
    for(int k=(cplx ? 2*len : len); k<nf; ++k)
    {
      const bool kept = (mode > 0 && !cplx && k < 2*len) ?
          std::memcmp(&rf[k], &of[k], sizeof(float)) == 0 :
          floatBits(rf[k]) == SENTINEL_BITS;
      if(kept) continue;
      ++nbad;
      if(nshown >= MAX_REPORT) continue;
      ++nshown;
      printf("      len=%-5d float %-5d WROTE PAST END evm=%.9g%s\n", len, k,
             rf[k], mname);
    }
  }
  const bool pass = (nbad == 0);
  if(!pass) ++gFailed;
  printf("  %-34s %s  [%ld values, NaN/zeros exact", name,
         pass ? "PASS" : "FAIL", compared);
  if(nbad) printf(", %d wrong", nbad);
  printf("]\n");
  if(nbad > nshown)
    printf("      ... %d further problems not shown\n", nbad - nshown);
}

//==============================================================================
// FM modulation (fmmod, and fmmods with an output scale)
//
//   runEvm(in, len, phscale, phstart, phincr, out)  calls the function under
//                                                   test (with the output
//                                                   scale of the run) and
//                                                   returns its result
//
// Element k must come out as scale*e^(j*2*pi*ph) with the phase in cycles
// ph = phscale*in[k] + phstart + k*phincr, computed in double with whole
// cycles removed before the float sin/cos (from phstart first), so a large
// phstart or k*phincr costs no accuracy.  The error must be within TOL_TRIG
// of |scale|.  The return value must be the phase of element len, phstart +
// len*phincr, with whole cycles removed (|r| < 1, within 1e-9 modulo whole
// cycles), and exactly phstart for len <= 0.  Every length runs out of place
// and in place (the real input in the front of the wider complex output,
// which the functions handle by running from the end to the start), and
// nothing may be written past out[len-1].
//
// The special values (zeros of both signs, quarter and half cycles, +-1,
// infinities, NaN, denormals and large values up to FLT_MAX, in every lane
// at lengths 1 to 40, out of place and in place) run with phscale = 1 and
// phstart = phincr = 0, so the phase is the input itself: an infinite or NaN
// one must give NaN in both parts.  A phstart past 2^49 cycles with
// fractional bits (2^50 + 0.25 and others) must come out right too: the
// AVX-512/AVX2 reduction is exact below 2^49 cycles, and whole cycles are
// removed from phstart first.  phscale*in[k] past 2^49 cycles is not tested,
// since there only the AVX and default versions are exact.  An infinite or
// NaN phstart or phincr must give NaN everywhere and as the return value.
//==============================================================================
template<class FEVM>
static void checkFmmod(const char* name, const std::vector<float>& src,
                       const double phscale, const double phstart,
                       const double phincr, const double scale, FEVM runEvm)
{
  ++gTests;
  long compared = 0;
  int  nbad = 0, nshown = 0;
  double worst = 0.0;
  const int n0 = static_cast<int>(src.size());
  std::vector<complex<float> > out(n0 + GUARD);
  const float* cur = src.data(); // the input of the current run
  int curLen = 0;
  auto report = [&](const int len, const int k, const char* what,
                    const cdbl got, const cdbl want, const char* mode)
  {
    ++nbad;
    if(nshown >= MAX_REPORT) return;
    ++nshown;
    printf("      len=%-5d element %-5d %s in=%.9g evm=%s want=%s%s\n", len, k,
           what, (k >= 0 && k < curLen) ? cur[k] : 0.0f,
           format(got, true).c_str(), format(want, true).c_str(), mode);
  };
  auto run = [&](const float* in, const int len, const int mode,
                 std::vector<complex<float> >& o, const double ps,
                 const double p0, const double p1)
  {
    for(size_t k=0; k<o.size(); ++k) o[k] = sentinel<complex<float> >();
    const int valid = (len > 0) ? len : 0;
    cur = in;
    curLen = valid;
    if(mode == 0) return runEvm(in, len, ps, p0, p1, o.data());
    std::memcpy(o.data(), in, valid*sizeof(float));
    return runEvm(reinterpret_cast<const float*>(o.data()), len, ps, p0, p1,
                  o.data());
  };
  // the expected element, whole cycles removed in double
  auto want = [&](const float x, const int k, const double ps,
                  const double p0, const double p1)
  {
    double ph = ps*static_cast<double>(x) + ((p0 - std::rint(p0)) + k*p1);
    ph -= std::rint(ph);
    return scale*std::polar(1.0, 2.0*M_PI*ph);
  };
  auto checkPast = [&](const std::vector<complex<float> >& o, const int len,
                       const char* mname)
  {
    const float* of = reinterpret_cast<const float*>(o.data());
    const int valid = (len > 0) ? len : 0;
    for(size_t k=2*valid; k<2*o.size(); ++k)
      if(!isSentinel(of[k]))
        report(len, static_cast<int>(k/2), "WROTE PAST END", toCdbl(o[k/2]),
               cdbl(), mname);
  };
  // ordinary values, every length, and the return value
  for(int mode=0; mode<2; ++mode)
  for(size_t li=0; li<gLengths.size(); ++li)
  {
    const int len = std::min(gLengths[li], n0);
    const int valid = (len > 0) ? len : 0;
    const char* mname = mode ? " [in place]" : "";
    const double r = run(src.data(), len, mode, out, phscale, phstart, phincr);
    // the phase of element len, compared modulo whole cycles
    const double wr = (len <= 0) ? phstart : phstart + len*phincr;
    const double dr = (r - wr) - std::rint(r - wr);
    ++compared;
    if(len <= 0 ? !(r == wr) : !(std::fabs(dr) <= 1.0e-9 &&
                                 std::fabs(r) < 1.0))
      report(len, -1, "RETURN", cdbl(r, 0.0), cdbl(wr, 0.0), mname);
    cur = src.data();
    for(int k=0; k<valid; ++k)
    {
      const cdbl w = want(src[k], k, phscale, phstart, phincr);
      const double e = std::abs(toCdbl(out[k]) - w)/std::fabs(scale);
      ++compared;
      if(e <= TOL_TRIG) worst = std::max(worst, e);
      else report(len, k, isSentinel(out[k]) ? "NOT WRITTEN" : "",
                  toCdbl(out[k]), w, mname);
    }
    checkPast(out, len, mname);
  }
  // special values, the phase is the input itself
  {
    const float inf = std::numeric_limits<float>::infinity();
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const std::vector<float> sp = {0.0f, -0.0f, 0.25f, -0.25f, 0.5f, -0.5f,
        1.0f, -1.0f, 0.125f, inf, -inf, nan, 1.0e-40f, -1.0e-40f, 65537.3f,
        -1.0e7f, 1.0e15f, -3.4e38f, 3.4e38f};
    const int nsp = static_cast<int>(sp.size());
    std::vector<complex<float> > o(64 + GUARD);
    for(int mode=0; mode<2; ++mode)
    for(int len=1; len<=40; ++len)
    for(int start=0; start<nsp; ++start)
    {
      std::vector<float> in(len);
      for(int k=0; k<len; ++k) in[k] = sp[(start + k) % nsp];
      const char* mname = mode ? " [in place]" : "";
      run(in.data(), len, mode, o, 1.0, 0.0, 0.0);
      cur = in.data();
      for(int k=0; k<len; ++k)
      {
        const cdbl got = toCdbl(o[k]);
        ++compared;
        if(bitsNan(in[k]) || bitsInf(in[k]))
        {
          if(!bitsNan(o[k].real()) || !bitsNan(o[k].imag()))
            report(len, k, "NaN/INF phase", got, cdbl(nan, nan), mname);
          continue;
        }
        const cdbl w = want(in[k], k, 1.0, 0.0, 0.0);
        const double e = std::abs(got - w)/std::fabs(scale);
        if(!(e <= TOL_TRIG + std::numeric_limits<float>::min()))
          report(len, k, "special value", got, w, mname);
      }
      checkPast(o, len, mname);
    }
    // phases from 2^49 cycles up, which still have fractional bits in
    // double (2^50 + 0.25 is a quarter cycle): phscale = 0, so every
    // element's phase is phstart
    const double big[4] = {1125899906842624.25, 2251799813685248.5,
                           4503599627370497.0, -1125899906842624.75};
    for(int c=0; c<4; ++c)
    {
      const int len = std::min(37, n0);
      run(src.data(), len, 0, o, 0.0, big[c], 0.0);
      for(int k=0; k<len; ++k)
      {
        const cdbl w = want(0.0f, k, 0.0, big[c], 0.0);
        const double e = std::abs(toCdbl(o[k]) - w)/std::fabs(scale);
        ++compared;
        if(!(e <= TOL_TRIG))
          report(len, k, "phstart past 2^49", toCdbl(o[k]), w, "");
      }
      checkPast(o, len, "");
    }
    // an infinite or NaN phstart or phincr
    const double dinf = std::numeric_limits<double>::infinity();
    const double dnan = std::numeric_limits<double>::quiet_NaN();
    const double p0s[3] = {dinf, dnan, 0.25}, p1s[3] = {0.1, 0.1, -dinf};
    for(int c=0; c<3; ++c)
    {
      const int len = std::min(37, n0);
      const double r = run(src.data(), len, 0, o, phscale, p0s[c], p1s[c]);
      ++compared;
      if(!dblNan(r))
        report(len, -1, "RETURN of a NaN/INF phase", cdbl(r, 0.0),
               cdbl(dnan, 0.0), "");
      for(int k=0; k<len; ++k)
      {
        ++compared;
        if(!bitsNan(o[k].real()) || !bitsNan(o[k].imag()))
          report(len, k, "NaN/INF phstart or phincr", toCdbl(o[k]),
                 cdbl(dnan, dnan), "");
      }
      checkPast(o, len, "");
    }
  }
  const bool pass = (nbad == 0);
  if(!pass) ++gFailed;
  printf("  %-34s %s  [%ld values, tol %.0e, worst rel %.3e", name,
         pass ? "PASS" : "FAIL", compared, TOL_TRIG, worst);
  if(nbad) printf(", %d wrong", nbad);
  printf("]\n");
  if(nbad > nshown)
    printf("      ... %d further problems not shown\n", nbad - nshown);
}

//==============================================================================
// Summation test driver (sum, sumSqDiff)
//
//   runEvm(a, b, len, res)  calls the EVM sum under test and puts its result
//                           in res[0]; res[1] is a guard that must be left
//                           untouched
//   term(x, y)              the term of an element in double, from a[k] and
//                           b[k] (a[k] for sum, |a[k]-b[k]|^2 for sumSqDiff)
//   perElem                 float terms each component of the result adds per
//                           element (2 for sumSqDiff, the squares of both
//                           parts of the difference)
//   extra                   rounding errors of a term before it is added, in
//                           units of 2^-24 (sumSqDiff: 2 from the difference,
//                           which the square doubles, and 1 from the square
//                           where the default version rounds it)
//   special(pass, want)     a value for a[p] that makes the result NaN or
//                           infinite, and in want the result it gives (its
//                           finite parts are not checked)
//
// As in checkDot(), the error of each component of the result is measured
// against the sum of the magnitudes of its terms, within
// (perElem*len + 1 + extra)*2^-24, and the inputs are copied into NaN-filled
// buffers, in2 one element off in1's alignment, so a read past the end shows
// up as a NaN result.  Negative and zero lengths must give exactly 0.  At
// lengths 1 to 40 each of the two special values goes through every position
// in turn, and the result must be the NaN or infinity it gives (bit tests),
// so every element, in every lane and in the partial last register, has to
// reach the sum.
//==============================================================================
template<class TR, class TA, class TB, class FEVM, class FTERM, class FSPEC>
static void checkSum(const char* name, const TA* const aBuf,
                     const TB* const bBuf, FEVM runEvm, FTERM term,
                     const int perElem, const int extra, FSPEC special)
{
  ++gTests;
  const bool   cplx = IsComplex<TR>::value;
  const double u    = std::ldexp(1.0, -24);
  std::vector<TA> a(MAXLEN + GUARD);
  std::vector<TB> b(MAXLEN + GUARD);

  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0;
  double worstRel = 0.0;

  auto fill = [&](const int valid)
  {
    for(int k=0; k<MAXLEN+GUARD; ++k)
    {
      a[k] = sentinel<TA>();
      b[k] = sentinel<TB>();
    }
    for(int k=0; k<valid; ++k) { a[k] = aBuf[k]; b[k+1] = bBuf[k]; }
  };
  auto stray = [&](const TR& guard, const int len)
  {
    if(isSentinel(guard)) return;
    ++nstray;
    if(nshown < MAX_REPORT)
    {
      ++nshown;
      printf("      len=%-5d WROTE PAST THE RESULT evm=%s\n", len,
             format(toCdbl(guard), cplx).c_str());
    }
  };

  for(size_t li=0; li<gLengths.size(); ++li)
  {
    const int len   = gLengths[li];
    const int valid = (len > 0) ? len : 0;
    fill(valid);
    TR res[2] = {sentinel<TR>(), sentinel<TR>()};
    runEvm(a.data(), b.data()+1, len, res);
    const cdbl got = toCdbl(res[0]);
    stray(res[1], len);

    // baseline, and the magnitude sums that scale each component's error
    cdbl   want  = 0.0;
    double absRe = 0.0, absIm = 0.0;
    for(int k=0; k<valid; ++k)
    {
      const cdbl t = term(aBuf[k], bBuf[k]);
      want  += t;
      absRe += std::fabs(t.real());
      absIm += std::fabs(t.imag());
    }
    const double tol   = (perElem*valid + 1 + extra)*u;
    const double dRe   = std::fabs(got.real() - want.real());
    const double dIm   = std::fabs(got.imag() - want.imag());
    const double rdiff = std::max((absRe > 0.0) ? dRe/absRe : dRe,
                                  (absIm > 0.0) ? dIm/absIm : dIm);
    ++compared;

    // --- value (NaN fails both tests) ---
    if((valid > 0) ? (rdiff <= tol) : (got == cdbl(0.0, 0.0)))
    {
      if(rdiff > worstRel) worstRel = rdiff;
      continue;
    }
    ++nbad;
    if(nshown < MAX_REPORT)
    {
      ++nshown;
      printf("      len=%-5d evm=%s baseline=%s reldiff=%.6g "
             "(of the sum of |terms|, tol %.3g)\n", len,
             format(got, cplx).c_str(), format(want, cplx).c_str(), rdiff,
             tol);
    }
  }

  // --- special values at every position ---
  // a part of want that is NaN or infinite must match exactly, compared as
  // floats with bit tests
  auto specOk = [](const float g, const double w)
  {
    const float wf = static_cast<float>(w);
    if(bitsNan(wf)) return bitsNan(g);
    if(bitsInf(wf)) return bitsInf(g) && bitsNeg(g) == bitsNeg(wf);
    return true;
  };
  for(int len=1; len<=40; ++len)
    for(int pass=0; pass<2; ++pass)
      for(int pos=0; pos<len; ++pos)
      {
        fill(len);
        cdbl want;
        a[pos] = special(pass, want);
        TR res[2] = {sentinel<TR>(), sentinel<TR>()};
        runEvm(a.data(), b.data()+1, len, res);
        stray(res[1], len);
        const cdbl got = toCdbl(res[0]);
        ++compared;
        if(specOk(static_cast<float>(got.real()), want.real()) &&
           specOk(static_cast<float>(got.imag()), want.imag())) continue;
        ++nbad;
        if(nshown < MAX_REPORT)
        {
          ++nshown;
          printf("      len=%-5d special value at %-3d evm=%s want=%s\n", len,
                 pos, format(got, cplx).c_str(), format(want, cplx).c_str());
        }
      }

  const bool pass = (nbad == 0 && nstray == 0);
  if(!pass) ++gFailed;

  printf("  %-34s %s  [%ld values, tol (%slen+%d)*2^-24, worst rel %.3e",
         name, pass ? "PASS" : "FAIL", compared, (perElem == 2) ? "2*" : "",
         1 + extra, worstRel);
  if(nbad)   printf(", %d wrong", nbad);
  if(nstray) printf(", %d past end", nstray);
  printf("]\n");
  if(nbad + nstray > nshown)
    printf("      ... %d further problems not shown\n",
           nbad + nstray - nshown);
}

//==============================================================================
// Special values for the magnitude functions
//
//   runEvm(in, len, out, cout)  calls the function under test, which writes
//                               its real results to out (out holds 1 before
//                               the call, for the sums) or its complex ones
//                               to cout
//
// kind says what the function computes from M = |in|: MS_SQ |in|^2, MS_MAG
// M, MS_NORM in/M, MS_BIAS M + 0.5, MS_SUM 0.75*M + 0.25*out.  Zeros of all
// four signs, ordinary values, values whose squared magnitude underflows,
// is a denormal or overflows, infinities, NaN and denormal parts go through
// every lane position, including the partial last register, at lengths 1 to
// 40.  The expected results follow rules rather than arithmetic on
// infinities: a NaN part gives NaN, an infinite part or an overflowing
// |in|^2 gives M = INF, and |in|^2 = 0 gives M = 0 exactly.  A denormal
// |in|^2 may give M = 0 (FTZ, and AVX/AVX2 rsqrt_ps without DAZ) or its
// square root, to within the rounding of the denormal (about 2^-149 of
// |in|^2).  The normalizations give 0 wherever |in|^2 is below FLT_MIN,
// INF or NaN: each part becomes a zero of its own sign, or NaN where it is
// NaN or INF.  NaN, infinities and those zeros must match exactly (bit tests,
// since icpx's fast math may drop NaN checks); otherwise the error must be
// within tol (for the normalizations, of 1).  Under DAZ, denormal parts act
// as zeros of the same sign.  infMayBeNaN lets a magnitude of INF come out as
// NaN instead: magr's AVX-512 version doesn't guard against an overflowing
// |in|^2 (INF*rsqrt14(INF) = INF*0), which is kept for speed.
//==============================================================================
enum MagSpecialKind { MS_SQ, MS_MAG, MS_NORM, MS_BIAS, MS_SUM };

// expected value of one real result: kind as in InvWant (0 = finite value,
// 1 = NaN, 2 = INF, 3 = zero), tiny = a denormal |in|^2 that may also give 0
struct MagWant { int kind; bool neg; double val; bool tiny; };

template<class FEVM>
static void checkMagSpecial(const char* name, const double tol,
                            const MagSpecialKind kind, FEVM runEvm,
                            const bool infMayBeNaN = false)
{
  ++gTests;
  const bool daz = (_mm_getcsr() & _MM_DENORMALS_ZERO_MASK) ==
                   _MM_DENORMALS_ZERO_ON;
  const float inf = std::numeric_limits<float>::infinity();
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const float fmin = std::numeric_limits<float>::min();
  const float fmax = std::numeric_limits<float>::max();
  const complex<float> vals[] = {
      complex<float>(0.0f, 0.0f), complex<float>(-0.0f, 0.0f),
      complex<float>(0.0f, -0.0f), complex<float>(-0.0f, -0.0f),
      complex<float>(1.0f, -1.0f), complex<float>(-3.0f, 4.0f),
      complex<float>(0.5f, 0.0f), complex<float>(0.0f, -2.0f),
      complex<float>(1.0e-25f, 0.0f), complex<float>(-1.0e-25f, 3.0e-26f),
      complex<float>(1.0e-20f, 1.0e-20f), complex<float>(5.0e-20f, 0.0f),
      complex<float>(2.0e-18f, -1.0e-18f), complex<float>(1.0e19f, 1.0e19f),
      complex<float>(2.0e19f, 0.0f), complex<float>(-3.0e19f, 1.0f),
      complex<float>(inf, 1.0f), complex<float>(1.0f, -inf),
      complex<float>(-inf, -inf), complex<float>(nan, 1.0f),
      complex<float>(1.0f, nan), complex<float>(1.0e-40f, 0.0f),
      complex<float>(-1.0e-40f, -2.0f)};
  const int nv = sizeof(vals)/sizeof(vals[0]);
  long compared = 0;
  int  nbad = 0, nshown = 0;
  std::vector<complex<float> > in(64), cout(64);
  std::vector<float> out(64);
  for(int len=1; len<=40; ++len)
  {
    for(int i=0; i<len; ++i)
    {
      in[i] = vals[(i + len)%nv];
      out[i] = 1.0f;
    }
    runEvm(in.data(), len, out.data(), cout.data());
    for(int i=0; i<len; ++i)
    {
      // the input as the kernels see it
      float re = in[i].real(), im = in[i].imag();
      if(daz && bitsDenormal(re)) re = bitsNeg(re) ? -0.0f : 0.0f;
      if(daz && bitsDenormal(im)) im = bitsNeg(im) ? -0.0f : 0.0f;
      const double m = static_cast<double>(re)*re + static_cast<double>(im)*im;
      // M = |in| as the functions see it
      MagWant w = {0, false, 0.0, false};
      if(bitsNan(re) || bitsNan(im)) w.kind = 1;
      else if(bitsInf(re) || bitsInf(im) || m > fmax) w.kind = 2;
      else if(static_cast<float>(m) == 0.0f) w.kind = 3;
      else
      {
        w.val = (kind == MS_SQ) ? m : std::sqrt(m);
        w.tiny = m < fmin;
      }
      const bool zeroScale = (w.kind != 0) || w.tiny;
      bool ok[2] = {true, true};
      float got[2] = {out[i], 0.0f};
      if(kind == MS_NORM)
      {
        got[0] = cout[i].real();
        got[1] = cout[i].imag();
        const float part[2] = {re, im};
        for(int c=0; c<2; ++c)
        {
          if(zeroScale)
          {
            if(bitsNan(part[c]) || bitsInf(part[c])) ok[c] = bitsNan(got[c]);
            else
              ok[c] = bitsZero(got[c]) && bitsNeg(got[c]) == bitsNeg(part[c]);
          }
          else
            ok[c] = !bitsNan(got[c]) && !bitsInf(got[c]) &&
                    std::fabs(got[c] - part[c]/w.val) <= tol;
        }
      }
      else
      {
        // the bias and the sum on top of M
        double add = 0.0, mul = 1.0;
        if(kind == MS_BIAS) add = 0.5;
        if(kind == MS_SUM) { mul = 0.75; add = 0.25; }
        const float g = got[0];
        if(w.kind == 1) ok[0] = bitsNan(g);
        else if(w.kind == 2)
          ok[0] = (bitsInf(g) && !bitsNeg(g)) || (infMayBeNaN && bitsNan(g));
        else
        {
          const double want = mul*w.val + add;
          // a denormal |in|^2 carries only absolute precision (about 2^-149
          // after its rounding), which its square root scales by 1/(2*M)
          const double extra = (w.tiny && kind != MS_SQ) ?
              mul*3.0e-45/(2.0*w.val) : 0.0;
          if(w.kind == 3 && add == 0.0) ok[0] = bitsZero(g) && !bitsNeg(g);
          else if(w.tiny && add == 0.0 && bitsZero(g) && !bitsNeg(g))
            ok[0] = true;
          else
            ok[0] = !bitsNan(g) && !bitsInf(g) &&
                    std::fabs(g - want) <=
                        tol*std::max(want, static_cast<double>(fmin)) + extra;
        }
      }
      ++compared;
      if(ok[0] && ok[1]) continue;
      ++nbad;
      if(nshown < MAX_REPORT)
      {
        ++nshown;
        if(kind == MS_NORM)
          printf("      len=%-3d in=(%g,%g) evm=(%.9g,%.9g)\n", len,
                 in[i].real(), in[i].imag(), got[0], got[1]);
        else
          printf("      len=%-3d in=(%g,%g) evm=%.9g\n", len, in[i].real(),
                 in[i].imag(), got[0]);
      }
    }
  }
  const bool pass = (nbad == 0);
  if(!pass) ++gFailed;
  printf("  %-34s %s  [%ld values, NaN/INF/0 exact", name,
         pass ? "PASS" : "FAIL", compared);
  if(nbad) printf(", %d wrong", nbad);
  printf("]\n");
  if(nbad > nshown)
    printf("      ... %d further problems not shown\n", nbad - nshown);
}

//==============================================================================
int main(int argc, char** argv)
{
  const unsigned int seed = (argc > 1)
      ? static_cast<unsigned int>(strtoul(argv[1], NULL, 0))
      : 20260927u;
  generateData(seed);

  printf("======================================================================\n");
  printf("EVM unit test\n");
  printf("  dispatch branch : %s\n", EVM::getOptLevel().c_str());
  printf("  random seed     : %u   (pass as argv[1] to reproduce)\n", seed);
  printf("  lengths tested  : %zu, from %d to %d\n",
         gLengths.size(), gLengths.front(), gLengths.back());
  printf("======================================================================\n");

  const float*          Af = gAf.data();
  const float*          Bf = gBf.data();
  const complex<float>* Ac = gAc.data();
  const complex<float>* Bc = gBc.data();
  const complex<float>* Zc = gZc.data();
  const float*          nullF = NULL;
  const float*          N0 = gNf[0].data();
  const float*          N1 = gNf[1].data();
  const float*          N2 = gNf[2].data();

  //--- min / max / minMax ---------------------------------------------------
  printf("\nmin / max / minMax\n");
  check<float, float, float>("min(float,float)", TOL_CMP, N0, N1,
    [](const float* a, const float* b, int n, float* o)
      { EVM::min(a, b, n, o); },
    [](int k){ return cdbl(std::min(gNf[0][k], gNf[1][k])); });

  check<float, float, float>("max(float,float)", TOL_CMP, N0, N1,
    [](const float* a, const float* b, int n, float* o)
      { EVM::max(a, b, n, o); },
    [](int k){ return cdbl(std::max(gNf[0][k], gNf[1][k])); });

  checkPair("minMax(float,float)", N0, N1,
    [](const float* a, const float* b, int n, float* mn, float* mx)
      { EVM::minMax(a, b, n, mn, mx); },
    [](int k){ return cdbl(std::min(gNf[0][k], gNf[1][k])); },
    [](int k){ return cdbl(std::max(gNf[0][k], gNf[1][k])); });

  // the reduction uses Af, the many ties in the comparison inputs would hide
  // a lost accumulator
  checkMinMax("minMax(float) of a buffer", Af);

  //--- zneg / zpos : negative or positive values set to zero ----------------
  printf("\nzneg / zpos\n");
  check<float, float, float>("zneg(float)", TOL_CMP, N0, nullF,
    [](const float* a, const float*, int n, float* o)
      { runZeroInPlace(true, a, n, o); },
    [](int k){ return cdbl(std::max(gNf[0][k], 0.0f)); });

  check<float, float, float>("zpos(float)", TOL_CMP, N0, nullF,
    [](const float* a, const float*, int n, float* o)
      { runZeroInPlace(false, a, n, o); },
    [](int k){ return cdbl(std::min(gNf[0][k], 0.0f)); });

  //--- signedAbsMin / signedAbsMax : signed value of min/max magnitude ------
  // in[0] and in[1] are the a and b buffers, so the in-place modes put out on
  // top of the first and second inputs, the rest come straight from gNf
  printf("\nsignedAbsMin / signedAbsMax\n");
  const int nbufsTested[] = {1, 2, NBUF};
  for(int w=0; w<2; ++w) // 0 = signedAbsMin, 1 = signedAbsMax
  {
    for(int t=0; t<3; ++t)
    {
      const int nb = nbufsTested[t];
      char name[64];
      snprintf(name, sizeof(name), "%s(%d buffer%s)",
               w ? "signedAbsMax" : "signedAbsMin", nb, (nb > 1) ? "s" : "");
      check<float, float, float>(name, TOL_CMP, N0, (nb > 1) ? N1 : nullF,
        [nb, w](const float* a, const float* b, int n, float* o)
          { const float* in[NBUF];
            in[0] = a;
            in[1] = b;
            for(int j=2; j<NBUF; ++j) in[j] = gNf[j].data();
            if(w) EVM::signedAbsMax(in, nb, n, o);
            else  EVM::signedAbsMin(in, nb, n, o); },
        [nb, w](int k){ return signedAbsRef(k, nb, w != 0); });
    }
  }
  // no buffers must be a no-op even for a positive length
  {
    ++gTests;
    std::vector<float> out(MAXLEN + GUARD, sentinel<float>());
    const float* in[NBUF];
    for(int j=0; j<NBUF; ++j) in[j] = gNf[j].data();
    for(int nb=-1; nb<=0; ++nb)
    {
      EVM::signedAbsMin(in, nb, MAXLEN, out.data());
      EVM::signedAbsMax(in, nb, MAXLEN, out.data());
    }
    bool ok = true;
    for(int k=0; k<MAXLEN+GUARD; ++k) if(!isSentinel(out[k])) ok = false;
    if(!ok) ++gFailed;
    printf("  %-34s %s  [nbufs 0 and -1 leave out untouched]\n",
           "signedAbsMin/Max(0 buffers)", ok ? "PASS" : "FAIL");
  }

  //--- threshCmp / blend : choose in1 or in2 by comparing with a threshold --
  printf("\nthreshCmp / blend\n");
  check<float, float, float>("threshCmp(float,float)", TOL_CMP, N0, N1,
    [](const float* a, const float* b, int n, float* o)
      { EVM::threshCmp(a, b, n, gTh, o); },
    [](int k){ return cdbl((gNf[0][k] >= gTh) ? gNf[0][k] : gNf[1][k]); });

  check<float, float, float>("blend(float,float,float)", TOL_CMP, N0, N1,
    [](const float* a, const float* b, int n, float* o)
      { EVM::blend(a, b, gNf[2].data(), n, gTh, o); },
    [](int k){ return cdbl((gNf[2][k] >= gTh) ? gNf[0][k] : gNf[1][k]); });

  // the same again with in3 as the a buffer, so the in-place modes also put
  // out on top of in3
  check<float, float, float>("blend(float,float,float) a=in3", TOL_CMP, N2, N1,
    [](const float* a, const float* b, int n, float* o)
      { EVM::blend(gNf[0].data(), b, a, n, gTh, o); },
    [](int k){ return cdbl((gNf[2][k] >= gTh) ? gNf[0][k] : gNf[1][k]); });

  //--- scale : out = in * scaleFactor ---------------------------------------
  printf("\nscale\n");
  check<float, float, float>("scale(float,float)", TOL_EXACT, Af, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::scale(a, gSf, n, o); },
    [](int k){ return toCdbl(gAf[k]) * cdbl(gSf); });

  check<complex<float>, float, complex<float> >("scale(complex,float)",
    TOL_EXACT, Ac, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::scale(a, gSf, n, o); },
    [](int k){ return toCdbl(gAc[k]) * cdbl(gSf); });

  check<float, float, complex<float> >("scale(float,complex)",
    TOL_EXACT, Af, nullF,
    [](const float* a, const float*, int n, complex<float>* o)
      { EVM::scale(a, gSc, n, o); },
    [](int k){ return toCdbl(gAf[k]) * toCdbl(gSc); });

  check<complex<float>, float, complex<float> >("scale(complex,complex)",
    TOL_EXACT, Ac, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::scale(a, gSc, n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gSc); });

  //--- mult : out = in1 * in2 -----------------------------------------------
  printf("\nmult\n");
  check<float, float, float>("mult(float,float)", TOL_EXACT, Af, Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::mult(a, b, n, o); },
    [](int k){ return toCdbl(gAf[k]) * toCdbl(gBf[k]); });

  check<complex<float>, complex<float>, complex<float> >("mult(complex,complex)",
    TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::mult(a, b, n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gBc[k]); });

  check<float, complex<float>, complex<float> >("mult(float,complex)",
    TOL_EXACT, Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::mult(a, b, n, o); },
    [](int k){ return toCdbl(gAf[k]) * toCdbl(gBc[k]); });

  //--- multc : out = in1 * conj(in2) ----------------------------------------
  printf("\nmultc\n");
  check<complex<float>, complex<float>, complex<float> >("multc(complex,complex)",
    TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::multc(a, b, n, o); },
    [](int k){ return toCdbl(gAc[k]) * std::conj(toCdbl(gBc[k])); });

  check<float, complex<float>, complex<float> >("multc(float,complex)",
    TOL_EXACT, Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::multc(a, b, n, o); },
    [](int k){ return toCdbl(gAf[k]) * std::conj(toCdbl(gBc[k])); });

  //--- mults : out = in1 * in2 * scale --------------------------------------
  printf("\nmults\n");
  check<float, float, float>("mults(float,float)", TOL_EXACT, Af, Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::mults(a, b, n, gSf, o); },
    [](int k){ return toCdbl(gAf[k]) * toCdbl(gBf[k]) * cdbl(gSf); });

  check<complex<float>, complex<float>, complex<float> >("mults(complex,complex)",
    TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::mults(a, b, n, gSf, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gBc[k]) * cdbl(gSf); });

  check<float, complex<float>, complex<float> >("mults(float,complex)",
    TOL_EXACT, Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::mults(a, b, n, gSf, o); },
    [](int k){ return toCdbl(gAf[k]) * toCdbl(gBc[k]) * cdbl(gSf); });

  //--- multcs : out = in1 * conj(in2) * scale -------------------------------
  printf("\nmultcs\n");
  check<complex<float>, complex<float>, complex<float> >("multcs(complex,complex)",
    TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::multcs(a, b, n, gSf, o); },
    [](int k){ return toCdbl(gAc[k]) * std::conj(toCdbl(gBc[k])) * cdbl(gSf); });

  //--- multe : out.real = in1 * in2, out.imag left as it was ---------------
  printf("\nmulte\n");
  checkExpand("multe(float,float)", TOL_EXACT, Af, Bf,
    [](const float* a, const float* b, int n, complex<float>* o)
      { EVM::multe(a, b, n, o); },
    [](int k){ return static_cast<double>(gAf[k])*gBf[k]; });

  //--- remult / immult : out = in1 * real(in2) or in1 * imag(in2) -----------
  printf("\nremult / immult\n");
  check<complex<float>, complex<float>, complex<float> >("remult(complex,complex)",
    TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::remult(a, b, n, o); },
    [](int k){ return toCdbl(gAc[k]) * static_cast<double>(gBc[k].real()); });

  check<complex<float>, complex<float>, complex<float> >("immult(complex,complex)",
    TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::immult(a, b, n, o); },
    [](int k){ return toCdbl(gAc[k]) * static_cast<double>(gBc[k].imag()); });

  //--- multSum : out = in1 * in2 + in3 --------------------------------------
  // error measured against |in1*in2| + |in3|, see the test driver; the a=in3
  // runs pass in3 as the a buffer so the in-place modes also put out on in3
  printf("\nmultSum\n");
  check<float, float, float>("multSum(float,float,float)", TOL_EXACT, Af, Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::multSum(a, b, gCf.data(), n, o); },
    [](int k){ return toCdbl(gAf[k]) * toCdbl(gBf[k]) + toCdbl(gCf[k]); },
    [](int k){ return std::abs(toCdbl(gAf[k]) * toCdbl(gBf[k])) +
                      std::abs(toCdbl(gCf[k])); });

  check<float, float, float>("multSum(float,float,float) a=in3", TOL_EXACT,
    gCf.data(), Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::multSum(gAf.data(), b, a, n, o); },
    [](int k){ return toCdbl(gAf[k]) * toCdbl(gBf[k]) + toCdbl(gCf[k]); },
    [](int k){ return std::abs(toCdbl(gAf[k]) * toCdbl(gBf[k])) +
                      std::abs(toCdbl(gCf[k])); });

  check<complex<float>, float, complex<float> >(
    "multSum(complex,float,complex)", TOL_EXACT, Ac, Bf,
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::multSum(a, b, gCc.data(), n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gBf[k]) + toCdbl(gCc[k]); },
    [](int k){ return std::abs(toCdbl(gAc[k]) * toCdbl(gBf[k])) +
                      std::abs(toCdbl(gCc[k])); });

  check<complex<float>, float, complex<float> >(
    "multSum(complex,float,complex) a=in3", TOL_EXACT, gCc.data(), Bf,
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::multSum(gAc.data(), b, a, n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gBf[k]) + toCdbl(gCc[k]); },
    [](int k){ return std::abs(toCdbl(gAc[k]) * toCdbl(gBf[k])) +
                      std::abs(toCdbl(gCc[k])); });

  check<complex<float>, complex<float>, complex<float> >(
    "multSum(complex,complex,complex)", TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::multSum(a, b, gCc.data(), n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gBc[k]) + toCdbl(gCc[k]); },
    [](int k){ return std::abs(toCdbl(gAc[k]) * toCdbl(gBc[k])) +
                      std::abs(toCdbl(gCc[k])); });

  check<complex<float>, complex<float>, complex<float> >(
    "multSum(complex,complex,complex) a=in3", TOL_EXACT, gCc.data(), Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::multSum(gAc.data(), b, a, n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gBc[k]) + toCdbl(gCc[k]); },
    [](int k){ return std::abs(toCdbl(gAc[k]) * toCdbl(gBc[k])) +
                      std::abs(toCdbl(gCc[k])); });

  //--- multcSum : out = in1 * conj(in2) + in3 -------------------------------
  // error measured against |in1*in2| + |in3| and a=in3 runs as for multSum
  printf("\nmultcSum\n");
  check<complex<float>, complex<float>, complex<float> >(
    "multcSum(complex,complex,complex)", TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::multcSum(a, b, gCc.data(), n, o); },
    [](int k){ return toCdbl(gAc[k]) * std::conj(toCdbl(gBc[k])) +
                      toCdbl(gCc[k]); },
    [](int k){ return std::abs(toCdbl(gAc[k]) * toCdbl(gBc[k])) +
                      std::abs(toCdbl(gCc[k])); });

  check<complex<float>, complex<float>, complex<float> >(
    "multcSum(complex,complex,complex) a=in3", TOL_EXACT, gCc.data(), Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::multcSum(gAc.data(), b, a, n, o); },
    [](int k){ return toCdbl(gAc[k]) * std::conj(toCdbl(gBc[k])) +
                      toCdbl(gCc[k]); },
    [](int k){ return std::abs(toCdbl(gAc[k]) * toCdbl(gBc[k])) +
                      std::abs(toCdbl(gCc[k])); });

  //--- multMult : out = in1 * in2 * in3 -------------------------------------
  printf("\nmultMult\n");
  check<complex<float>, complex<float>, complex<float> >(
    "multMult(complex,complex,float)", TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, complex<float>* o)
      { EVM::multMult(a, b, gCf.data(), n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gBc[k]) * toCdbl(gCf[k]); });

  //--- multDiv / multDivr / multDivnr : out = in1 * in2 / in3 ---------------
  // the divisor in3 is gBf, which is kept away from zero
  printf("\nmultDiv / multDivr / multDivnr\n");
  check<complex<float>, float, complex<float> >(
    "multDiv(complex,float,float)", TOL_EXACT, Ac, Af,
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::multDiv(a, b, gBf.data(), n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gAf[k]) / toCdbl(gBf[k]); });

  check<complex<float>, float, complex<float> >(
    "multDivr(complex,float,float)", TOL_RCP, Ac, Af,
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::multDivr(a, b, gBf.data(), n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gAf[k]) / toCdbl(gBf[k]); });

  check<complex<float>, float, complex<float> >(
    "multDivnr(complex,float,float)", TOL_NR, Ac, Af,
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::multDivnr(a, b, gBf.data(), n, o); },
    [](int k){ return toCdbl(gAc[k]) * toCdbl(gAf[k]) / toCdbl(gBf[k]); });

  //--- div / divnr / divr : out = in1 / in2 ---------------------------------
  // Same maths three ways: true divide, reciprocal + one Newton-Raphson
  // refinement stage, and bare reciprocal.  Only the tolerance differs.
  struct DivVariant { const char* suffix; double tol; int which; };
  const DivVariant variants[] = {
    {"div",   TOL_EXACT, 0},
    {"divnr", TOL_NR,    1},
    {"divr",  TOL_RCP,   2},
  };

  for(int v=0; v<3; ++v)
  {
    const DivVariant& V = variants[v];
    const int w = V.which;
    printf("\n%s\n", V.suffix);

    check<float, float, float>((std::string(V.suffix)+"(float,float)").c_str(),
      V.tol, Af, Bf,
      [w](const float* a, const float* b, int n, float* o)
        { if(w==0) EVM::div(a,b,n,o); else if(w==1) EVM::divnr(a,b,n,o);
          else EVM::divr(a,b,n,o); },
      [](int k){ return toCdbl(gAf[k]) / toCdbl(gBf[k]); });

    check<complex<float>, complex<float>, complex<float> >(
      (std::string(V.suffix)+"(complex,complex)").c_str(), V.tol, Ac, Bc,
      [w](const complex<float>* a, const complex<float>* b, int n,
          complex<float>* o)
        { if(w==0) EVM::div(a,b,n,o); else if(w==1) EVM::divnr(a,b,n,o);
          else EVM::divr(a,b,n,o); },
      [](int k){ return toCdbl(gAc[k]) / toCdbl(gBc[k]); });

    check<float, complex<float>, complex<float> >(
      (std::string(V.suffix)+"(float,complex)").c_str(), V.tol, Af, Bc,
      [w](const float* a, const complex<float>* b, int n, complex<float>* o)
        { if(w==0) EVM::div(a,b,n,o); else if(w==1) EVM::divnr(a,b,n,o);
          else EVM::divr(a,b,n,o); },
      [](int k){ return toCdbl(gAf[k]) / toCdbl(gBc[k]); });

    check<complex<float>, float, complex<float> >(
      (std::string(V.suffix)+"(complex,float)").c_str(), V.tol, Ac, Bf,
      [w](const complex<float>* a, const float* b, int n, complex<float>* o)
        { if(w==0) EVM::div(a,b,n,o); else if(w==1) EVM::divnr(a,b,n,o);
          else EVM::divr(a,b,n,o); },
      [](int k){ return toCdbl(gAc[k]) / toCdbl(gBf[k]); });
  }

  //--- divc1 family : out = conj(in1) / in2 ---------------------------------
  const DivVariant c1variants[] = {
    {"divc1",   TOL_EXACT, 0},
    {"divc1nr", TOL_NR,    1},
    {"divc1r",  TOL_RCP,   2},
  };
  for(int v=0; v<3; ++v)
  {
    const DivVariant& V = c1variants[v];
    const int w = V.which;
    printf("\n%s  (first input conjugated)\n", V.suffix);

    check<complex<float>, complex<float>, complex<float> >(
      (std::string(V.suffix)+"(complex,complex)").c_str(), V.tol, Ac, Bc,
      [w](const complex<float>* a, const complex<float>* b, int n,
          complex<float>* o)
        { if(w==0) EVM::divc1(a,b,n,o); else if(w==1) EVM::divc1nr(a,b,n,o);
          else EVM::divc1r(a,b,n,o); },
      [](int k){ return std::conj(toCdbl(gAc[k])) / toCdbl(gBc[k]); });

    check<complex<float>, float, complex<float> >(
      (std::string(V.suffix)+"(complex,float)").c_str(), V.tol, Ac, Bf,
      [w](const complex<float>* a, const float* b, int n, complex<float>* o)
        { if(w==0) EVM::divc1(a,b,n,o); else if(w==1) EVM::divc1nr(a,b,n,o);
          else EVM::divc1r(a,b,n,o); },
      [](int k){ return std::conj(toCdbl(gAc[k])) / toCdbl(gBf[k]); });
  }

  //--- divc2 family : out = in1 / conj(in2) ---------------------------------
  const DivVariant c2variants[] = {
    {"divc2",   TOL_EXACT, 0},
    {"divc2nr", TOL_NR,    1},
    {"divc2r",  TOL_RCP,   2},
  };
  for(int v=0; v<3; ++v)
  {
    const DivVariant& V = c2variants[v];
    const int w = V.which;
    printf("\n%s  (second input conjugated)\n", V.suffix);

    check<float, complex<float>, complex<float> >(
      (std::string(V.suffix)+"(float,complex)").c_str(), V.tol, Af, Bc,
      [w](const float* a, const complex<float>* b, int n, complex<float>* o)
        { if(w==0) EVM::divc2(a,b,n,o); else if(w==1) EVM::divc2nr(a,b,n,o);
          else EVM::divc2r(a,b,n,o); },
      [](int k){ return toCdbl(gAf[k]) / std::conj(toCdbl(gBc[k])); });

    check<complex<float>, complex<float>, complex<float> >(
      (std::string(V.suffix)+"(complex,complex)").c_str(), V.tol, Ac, Bc,
      [w](const complex<float>* a, const complex<float>* b, int n,
          complex<float>* o)
        { if(w==0) EVM::divc2(a,b,n,o); else if(w==1) EVM::divc2nr(a,b,n,o);
          else EVM::divc2r(a,b,n,o); },
      [](int k){ return toCdbl(gAc[k]) / std::conj(toCdbl(gBc[k])); });
  }

  //--- dive / divre / divnre : out.real = in1 / in2, out.imag as it was ----
  printf("\ndive / divre / divnre\n");
  checkExpand("dive(float,float)", TOL_EXACT, Af, Bf,
    [](const float* a, const float* b, int n, complex<float>* o)
      { EVM::dive(a, b, n, o); },
    [](int k){ return static_cast<double>(gAf[k])/gBf[k]; });

  checkExpand("divre(float,float)", TOL_RCP, Af, Bf,
    [](const float* a, const float* b, int n, complex<float>* o)
      { EVM::divre(a, b, n, o); },
    [](int k){ return static_cast<double>(gAf[k])/gBf[k]; });

  checkExpand("divnre(float,float)", TOL_NR, Af, Bf,
    [](const float* a, const float* b, int n, complex<float>* o)
      { EVM::divnre(a, b, n, o); },
    [](int k){ return static_cast<double>(gAf[k])/gBf[k]; });

  //--- divg / divrg / divnrg : out[k] = in1[k*stride] / in2[k*stride] -------
  // the eps overloads divide by max(in2, eps) instead
  printf("\ndivg / divrg / divnrg\n");
  checkGather("divg(float,float)", TOL_EXACT,
    [](const float* a, const float* b, int n, int s, float* o)
      { EVM::divg(a, b, n, s, o); },
    [](float a, float b){ return static_cast<double>(a)/b; });

  checkGather("divg(float,float) eps", TOL_EXACT,
    [](const float* a, const float* b, int n, int s, float* o)
      { EVM::divg(a, b, n, s, o, GEPS); },
    [](float a, float b){ return static_cast<double>(a)/std::max(b, GEPS); });

  checkGather("divrg(float,float)", TOL_RCP,
    [](const float* a, const float* b, int n, int s, float* o)
      { EVM::divrg(a, b, n, s, o); },
    [](float a, float b){ return static_cast<double>(a)/b; });

  checkGather("divrg(float,float) eps", TOL_RCP,
    [](const float* a, const float* b, int n, int s, float* o)
      { EVM::divrg(a, b, n, s, o, GEPS); },
    [](float a, float b){ return static_cast<double>(a)/std::max(b, GEPS); });

  checkGather("divnrg(float,float)", TOL_NR,
    [](const float* a, const float* b, int n, int s, float* o)
      { EVM::divnrg(a, b, n, s, o); },
    [](float a, float b){ return static_cast<double>(a)/b; });

  checkGather("divnrg(float,float) eps", TOL_NR,
    [](const float* a, const float* b, int n, int s, float* o)
      { EVM::divnrg(a, b, n, s, o, GEPS); },
    [](float a, float b){ return static_cast<double>(a)/std::max(b, GEPS); });

  //--- divSqrt / divrSqrt / divnrSqrt : out = in1 / sqrt(max(in2, eps)) -----
  // eps defaults to 0, those runs use the positive divisors gPf
  printf("\ndivSqrt / divrSqrt / divnrSqrt\n");
  check<complex<float>, float, complex<float> >(
    "divSqrt(complex,float)", TOL_EXACT, Ac, gPf.data(),
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::divSqrt(a, b, n, o); },
    [](int k){ return toCdbl(gAc[k]) / std::sqrt(toCdbl(gPf[k])); });

  check<complex<float>, float, complex<float> >(
    "divSqrt(complex,float) eps", TOL_EXACT, Ac, Bf,
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::divSqrt(a, b, n, o, GEPS); },
    [](int k)
      { return toCdbl(gAc[k]) / std::sqrt(toCdbl(std::max(gBf[k], GEPS))); });

  check<complex<float>, float, complex<float> >(
    "divrSqrt(complex,float)", TOL_RCP, Ac, gPf.data(),
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::divrSqrt(a, b, n, o); },
    [](int k){ return toCdbl(gAc[k]) / std::sqrt(toCdbl(gPf[k])); });

  check<complex<float>, float, complex<float> >(
    "divrSqrt(complex,float) eps", TOL_RCP, Ac, Bf,
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::divrSqrt(a, b, n, o, GEPS); },
    [](int k)
      { return toCdbl(gAc[k]) / std::sqrt(toCdbl(std::max(gBf[k], GEPS))); });

  check<complex<float>, float, complex<float> >(
    "divnrSqrt(complex,float)", TOL_NR, Ac, gPf.data(),
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::divnrSqrt(a, b, n, o); },
    [](int k){ return toCdbl(gAc[k]) / std::sqrt(toCdbl(gPf[k])); });

  check<complex<float>, float, complex<float> >(
    "divnrSqrt(complex,float) eps", TOL_NR, Ac, Bf,
    [](const complex<float>* a, const float* b, int n, complex<float>* o)
      { EVM::divnrSqrt(a, b, n, o, GEPS); },
    [](int k)
      { return toCdbl(gAc[k]) / std::sqrt(toCdbl(std::max(gBf[k], GEPS))); });

  //--- inv / invr / invnr / invs / invrs / invnrs : scale/max(in, eps) ----
  // the inputs are the divisors gBf, scale is gSf
  printf("\ninv / invr / invnr / invs / invrs / invnrs\n");
  check<float, float, float>("inv(float)", TOL_EXACT, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::inv(a, n, o); },
    [](int k){ return cdbl(1.0) / toCdbl(gBf[k]); });

  check<float, float, float>("inv(float) eps", TOL_EXACT, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::inv(a, n, o, GEPS); },
    [](int k){ return cdbl(1.0) / toCdbl(std::max(gBf[k], GEPS)); });

  check<float, float, float>("invr(float)", TOL_RCP, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invr(a, n, o); },
    [](int k){ return cdbl(1.0) / toCdbl(gBf[k]); });

  check<float, float, float>("invr(float) eps", TOL_RCP, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invr(a, n, o, GEPS); },
    [](int k){ return cdbl(1.0) / toCdbl(std::max(gBf[k], GEPS)); });

  check<float, float, float>("invnr(float)", TOL_NR, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invnr(a, n, o); },
    [](int k){ return cdbl(1.0) / toCdbl(gBf[k]); });

  check<float, float, float>("invnr(float) eps", TOL_NR, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invnr(a, n, o, GEPS); },
    [](int k){ return cdbl(1.0) / toCdbl(std::max(gBf[k], GEPS)); });

  check<float, float, float>("invs(float)", TOL_EXACT, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invs(a, n, gSf, o); },
    [](int k){ return toCdbl(gSf) / toCdbl(gBf[k]); });

  check<float, float, float>("invs(float) eps", TOL_EXACT, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invs(a, n, gSf, o, GEPS); },
    [](int k){ return toCdbl(gSf) / toCdbl(std::max(gBf[k], GEPS)); });

  check<float, float, float>("invrs(float)", TOL_RCP, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invrs(a, n, gSf, o); },
    [](int k){ return toCdbl(gSf) / toCdbl(gBf[k]); });

  check<float, float, float>("invrs(float) eps", TOL_RCP, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invrs(a, n, gSf, o, GEPS); },
    [](int k){ return toCdbl(gSf) / toCdbl(std::max(gBf[k], GEPS)); });

  check<float, float, float>("invnrs(float)", TOL_NR, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invnrs(a, n, gSf, o); },
    [](int k){ return toCdbl(gSf) / toCdbl(gBf[k]); });

  check<float, float, float>("invnrs(float) eps", TOL_NR, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invnrs(a, n, gSf, o, GEPS); },
    [](int k){ return toCdbl(gSf) / toCdbl(std::max(gBf[k], GEPS)); });

  //--- invSqrt / invrSqrt / invnrSqrt / invSqrts / invrSqrts / invnrSqrts --
  // out = scale/sqrt(max(in, eps))
  // eps defaults to 0, those runs use the positive divisors gPf
  printf("\ninvSqrt / invrSqrt / invnrSqrt / invSqrts / invrSqrts / "
         "invnrSqrts\n");
  check<float, float, float>("invSqrt(float)", TOL_EXACT, gPf.data(), nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invSqrt(a, n, o); },
    [](int k){ return cdbl(1.0) / std::sqrt(toCdbl(gPf[k])); });

  check<float, float, float>("invSqrt(float) eps", TOL_EXACT, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invSqrt(a, n, o, GEPS); },
    [](int k){ return cdbl(1.0) / std::sqrt(toCdbl(std::max(gBf[k], GEPS))); });

  check<float, float, float>("invrSqrt(float)", TOL_RCP, gPf.data(), nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invrSqrt(a, n, o); },
    [](int k){ return cdbl(1.0) / std::sqrt(toCdbl(gPf[k])); });

  check<float, float, float>("invrSqrt(float) eps", TOL_RCP, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invrSqrt(a, n, o, GEPS); },
    [](int k){ return cdbl(1.0) / std::sqrt(toCdbl(std::max(gBf[k], GEPS))); });

  check<float, float, float>("invnrSqrt(float)", TOL_NR, gPf.data(), nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invnrSqrt(a, n, o); },
    [](int k){ return cdbl(1.0) / std::sqrt(toCdbl(gPf[k])); });

  check<float, float, float>("invnrSqrt(float) eps", TOL_NR, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invnrSqrt(a, n, o, GEPS); },
    [](int k){ return cdbl(1.0) / std::sqrt(toCdbl(std::max(gBf[k], GEPS))); });

  check<float, float, float>("invSqrts(float)", TOL_EXACT, gPf.data(), nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invSqrts(a, n, gSf, o); },
    [](int k){ return toCdbl(gSf) / std::sqrt(toCdbl(gPf[k])); });

  check<float, float, float>("invSqrts(float) eps", TOL_EXACT, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invSqrts(a, n, gSf, o, GEPS); },
    [](int k)
      { return toCdbl(gSf) / std::sqrt(toCdbl(std::max(gBf[k], GEPS))); });

  check<float, float, float>("invrSqrts(float)", TOL_RCP, gPf.data(), nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invrSqrts(a, n, gSf, o); },
    [](int k){ return toCdbl(gSf) / std::sqrt(toCdbl(gPf[k])); });

  check<float, float, float>("invrSqrts(float) eps", TOL_RCP, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invrSqrts(a, n, gSf, o, GEPS); },
    [](int k)
      { return toCdbl(gSf) / std::sqrt(toCdbl(std::max(gBf[k], GEPS))); });

  check<float, float, float>("invnrSqrts(float)", TOL_NR, gPf.data(), nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invnrSqrts(a, n, gSf, o); },
    [](int k){ return toCdbl(gSf) / std::sqrt(toCdbl(gPf[k])); });

  check<float, float, float>("invnrSqrts(float) eps", TOL_NR, Bf, nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::invnrSqrts(a, n, gSf, o, GEPS); },
    [](int k)
      { return toCdbl(gSf) / std::sqrt(toCdbl(std::max(gBf[k], GEPS))); });

  //--- inversion special values -------------------------------------------
  // icpx's fast math turns 1.0f/sqrtf(x) in the default versions into an
  // rsqrt approximation and a Newton-Raphson step, which gives NaN for 0 and
  // INF, so the square roots are only checked on the SIMD paths with icpx
  printf("\ninversion special values\n");
  {
    const std::string lvl = EVM::getOptLevel();
#if defined(__INTEL_LLVM_COMPILER)
    const bool sqrtSpecial = (lvl == "AVX512" || lvl == "AVX2" ||
                              lvl == "AVX");
#else
    const bool sqrtSpecial = true;
#endif
    checkInvSpecial("inv(float) special", TOL_EXACT, false, false, false,
      [](const float* a, int n, float* o, bool e, float)
      {
        if(e) EVM::inv(a, n, o, GEPS);
        else  EVM::inv(a, n, o);
      });
    checkInvSpecial("invr(float) special", TOL_RCP, true, false, false,
      [](const float* a, int n, float* o, bool e, float)
      {
        if(e) EVM::invr(a, n, o, GEPS);
        else  EVM::invr(a, n, o);
      });
    checkInvSpecial("invnr(float) special", TOL_NR, true, false, false,
      [](const float* a, int n, float* o, bool e, float)
      {
        if(e) EVM::invnr(a, n, o, GEPS);
        else  EVM::invnr(a, n, o);
      });
    checkInvSpecial("invs(float) special", TOL_EXACT, false, false, true,
      [](const float* a, int n, float* o, bool e, float s)
      {
        if(e) EVM::invs(a, n, s, o, GEPS);
        else  EVM::invs(a, n, s, o);
      });
    checkInvSpecial("invrs(float) special", TOL_RCP, true, false, true,
      [](const float* a, int n, float* o, bool e, float s)
      {
        if(e) EVM::invrs(a, n, s, o, GEPS);
        else  EVM::invrs(a, n, s, o);
      });
    checkInvSpecial("invnrs(float) special", TOL_NR, true, false, true,
      [](const float* a, int n, float* o, bool e, float s)
      {
        if(e) EVM::invnrs(a, n, s, o, GEPS);
        else  EVM::invnrs(a, n, s, o);
      });
    if(sqrtSpecial)
    {
      checkInvSpecial("invSqrt(float) special", TOL_EXACT, false, true, false,
        [](const float* a, int n, float* o, bool e, float)
        {
          if(e) EVM::invSqrt(a, n, o, GEPS);
          else  EVM::invSqrt(a, n, o);
        });
      checkInvSpecial("invrSqrt(float) special", TOL_RCP, true, true, false,
        [](const float* a, int n, float* o, bool e, float)
        {
          if(e) EVM::invrSqrt(a, n, o, GEPS);
          else  EVM::invrSqrt(a, n, o);
        });
      checkInvSpecial("invnrSqrt(float) special", TOL_NR, true, true, false,
        [](const float* a, int n, float* o, bool e, float)
        {
          if(e) EVM::invnrSqrt(a, n, o, GEPS);
          else  EVM::invnrSqrt(a, n, o);
        });
      checkInvSpecial("invSqrts(float) special", TOL_EXACT, false, true, true,
        [](const float* a, int n, float* o, bool e, float s)
        {
          if(e) EVM::invSqrts(a, n, s, o, GEPS);
          else  EVM::invSqrts(a, n, s, o);
        });
      checkInvSpecial("invrSqrts(float) special", TOL_RCP, true, true, true,
        [](const float* a, int n, float* o, bool e, float s)
        {
          if(e) EVM::invrSqrts(a, n, s, o, GEPS);
          else  EVM::invrSqrts(a, n, s, o);
        });
      checkInvSpecial("invnrSqrts(float) special", TOL_NR, true, true, true,
        [](const float* a, int n, float* o, bool e, float s)
        {
          if(e) EVM::invnrSqrts(a, n, s, o, GEPS);
          else  EVM::invnrSqrts(a, n, s, o);
        });
    }
    else
      printf("  %-34s skipped (icpx default versions)\n",
             "invSqrt ... invnrSqrts special");
  }

  //--- dotp : sum of in1 * in2 ----------------------------------------------
  // the reference versions write the result to out, the by-value wrappers
  // return it
  printf("\ndotp\n");
  checkDot("dotp(float,float)", Af, Bf,
    [](const float* a, const float* b, int n, float* r)
      { EVM::dotp(a, b, n, r[0]); });

  checkDot("dotp(float,float) by value", Af, Bf,
    [](const float* a, const float* b, int n, float* r)
      { r[0] = EVM::dotp(a, b, n); });

  checkDot("dotp(complex,complex)", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotp(a, b, n, r[0]); });

  checkDot("dotp(complex,complex) by value", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { r[0] = EVM::dotp(a, b, n); });

  checkDot("dotp(float,complex)", Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* r)
      { EVM::dotp(a, b, n, r[0]); });

  checkDot("dotp(float,complex) by value", Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* r)
      { r[0] = EVM::dotp(a, b, n); });

  //--- dotpc : sum of in1 * conj(in2) ---------------------------------------
  printf("\ndotpc  (second input conjugated)\n");
  checkDot("dotpc(complex,complex)", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpc(a, b, n, r[0]); }, true);

  checkDot("dotpc(complex,complex) by value", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { r[0] = EVM::dotpc(a, b, n); }, true);

  checkDot("dotpc(float,complex)", Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* r)
      { EVM::dotpc(a, b, n, r[0]); }, true);

  checkDot("dotpc(float,complex) by value", Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* r)
      { r[0] = EVM::dotpc(a, b, n); }, true);

  //--- dotpRev : sum of in1[k] * in2[len-1-k] -------------------------------
  printf("\ndotpRev  (second input reversed)\n");
  checkDot("dotpRev(float,float)", Af, Bf,
    [](const float* a, const float* b, int n, float* r)
      { EVM::dotpRev(a, b, n, r[0]); }, false, true);

  checkDot("dotpRev(float,float) by value", Af, Bf,
    [](const float* a, const float* b, int n, float* r)
      { r[0] = EVM::dotpRev(a, b, n); }, false, true);

  checkDot("dotpRev(complex,complex)", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpRev(a, b, n, r[0]); }, false, true);

  checkDot("dotpRev(complex,complex) by value", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { r[0] = EVM::dotpRev(a, b, n); }, false, true);

  checkDot("dotpRev(float,complex)", Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* r)
      { EVM::dotpRev(a, b, n, r[0]); }, false, true);

  checkDot("dotpRev(float,complex) by value", Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* r)
      { r[0] = EVM::dotpRev(a, b, n); }, false, true);

  //--- dotpcRev : sum of in1[k] * conj(in2[len-1-k]) ------------------------
  printf("\ndotpcRev  (second input conjugated and reversed)\n");
  checkDot("dotpcRev(complex,complex)", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpcRev(a, b, n, r[0]); }, true, true);

  checkDot("dotpcRev(complex,complex) by value", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { r[0] = EVM::dotpcRev(a, b, n); }, true, true);

  checkDot("dotpcRev(float,complex)", Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* r)
      { EVM::dotpcRev(a, b, n, r[0]); }, true, true);

  checkDot("dotpcRev(float,complex) by value", Af, Bc,
    [](const float* a, const complex<float>* b, int n, complex<float>* r)
      { r[0] = EVM::dotpcRev(a, b, n); }, true, true);

  //--- dotpn / dotpnr / dotpcn / dotpcnr : sum of p/|p|, p = in1*in2 -------
  // (in2 conjugated for dotpcn/dotpcnr).  Each also runs with zeros in every
  // third element of in1 (gZc), whose products must add nothing
  printf("\ndotpn / dotpnr / dotpcn / dotpcnr  (each product normalized)\n");
  checkDotN("dotpn(complex,complex)", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpn(a, b, n, r[0]); }, false, TOL_EXACT);

  checkDotN("dotpn(complex,complex) by value", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { r[0] = EVM::dotpn(a, b, n); }, false, TOL_EXACT);

  checkDotN("dotpn(complex,complex) zeros", Zc, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpn(a, b, n, r[0]); }, false, TOL_EXACT);

  checkDotN("dotpnr(complex,complex)", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpnr(a, b, n, r[0]); }, false, TOL_RCP);

  checkDotN("dotpnr(complex,complex) by value", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { r[0] = EVM::dotpnr(a, b, n); }, false, TOL_RCP);

  checkDotN("dotpnr(complex,complex) zeros", Zc, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpnr(a, b, n, r[0]); }, false, TOL_RCP);

  checkDotN("dotpcn(complex,complex)", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpcn(a, b, n, r[0]); }, true, TOL_EXACT);

  checkDotN("dotpcn(complex,complex) by value", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { r[0] = EVM::dotpcn(a, b, n); }, true, TOL_EXACT);

  checkDotN("dotpcn(complex,complex) zeros", Zc, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpcn(a, b, n, r[0]); }, true, TOL_EXACT);

  checkDotN("dotpcnr(complex,complex)", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpcnr(a, b, n, r[0]); }, true, TOL_RCP);

  checkDotN("dotpcnr(complex,complex) by value", Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { r[0] = EVM::dotpcnr(a, b, n); }, true, TOL_RCP);

  checkDotN("dotpcnr(complex,complex) zeros", Zc, Bc,
    [](const complex<float>* a, const complex<float>* b, int n,
       complex<float>* r) { EVM::dotpcnr(a, b, n, r[0]); }, true, TOL_RCP);

  //--- pow2 / pow4 / pow8 : out = in^2, in^4, in^8 --------------------------
  printf("\npow2 / pow4 / pow8\n");
  check<float, float, float>("pow2(float)", TOL_EXACT, Af, nullF,
    [](const float* a, const float*, int n, float* o) { EVM::pow2(a, n, o); },
    [](int k){ return powSq(toCdbl(gAf[k]), 1); });

  check<complex<float>, float, complex<float> >("pow2(complex)", TOL_EXACT,
    Ac, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::pow2(a, n, o); },
    [](int k){ return powSq(toCdbl(gAc[k]), 1); });

  check<float, float, float>("pow4(float)", TOL_EXACT, Af, nullF,
    [](const float* a, const float*, int n, float* o) { EVM::pow4(a, n, o); },
    [](int k){ return powSq(toCdbl(gAf[k]), 2); });

  check<complex<float>, float, complex<float> >("pow4(complex)", TOL_EXACT,
    Ac, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::pow4(a, n, o); },
    [](int k){ return powSq(toCdbl(gAc[k]), 2); });

  check<float, float, float>("pow8(float)", TOL_EXACT, Af, nullF,
    [](const float* a, const float*, int n, float* o) { EVM::pow8(a, n, o); },
    [](int k){ return powSq(toCdbl(gAf[k]), 3); });

  check<complex<float>, float, complex<float> >("pow8(complex)", TOL_EXACT,
    Ac, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::pow8(a, n, o); },
    [](int k){ return powSq(toCdbl(gAc[k]), 3); });

  //--- magSq / mag : out = |in|^2, |in| -------------------------------------
  printf("\nmagSq / mag\n");
  check<complex<float>, float, float>("magSq(complex)", TOL_EXACT, Ac, nullF,
    [](const complex<float>* a, const float*, int n, float* o)
      { EVM::magSq(a, n, o); },
    [](int k){ return cdbl(std::norm(toCdbl(gAc[k]))); });

  check<float, float, float>("magSq(float)", TOL_EXACT, Af, nullF,
    [](const float* a, const float*, int n, float* o) { EVM::magSq(a, n, o); },
    [](int k){ return powSq(toCdbl(gAf[k]), 1); });

  check<float, float, float>("mag(float)", TOL_EXACT, Af, nullF,
    [](const float* a, const float*, int n, float* o) { EVM::mag(a, n, o); },
    [](int k){ return cdbl(std::fabs(static_cast<double>(gAf[k]))); });

  check<float, float, float>("abs(float)", TOL_EXACT, Af, nullF,
    [](const float* a, const float*, int n, float* o) { EVM::abs(a, n, o); },
    [](int k){ return cdbl(std::fabs(static_cast<double>(gAf[k]))); });

  check<complex<float>, float, float>("mag(complex)", TOL_EXACT, Zc, nullF,
    [](const complex<float>* a, const float*, int n, float* o)
      { EVM::mag(a, n, o); },
    [](int k){ return cdbl(std::abs(toCdbl(gZc[k]))); });

  check<complex<float>, float, float>("magr(complex)", TOL_RCP, Zc, nullF,
    [](const complex<float>* a, const float*, int n, float* o)
      { EVM::magr(a, n, o); },
    [](int k){ return cdbl(std::abs(toCdbl(gZc[k]))); });

  //--- nmag / nmagr / nmagnr : out = in/|in| --------------------------------
  // run on gZc: an element whose squared magnitude is below FLT_MIN gives 0
  printf("\nnmag / nmagr / nmagnr\n");
  check<complex<float>, float, complex<float> >("nmag(complex)", TOL_EXACT,
    Zc, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::nmag(a, n, o); },
    [](int k)
    {
      const cdbl z = toCdbl(gZc[k]);
      return (std::norm(z) >= FLT_MIN) ? z/std::abs(z) : cdbl(0.0);
    });

  check<complex<float>, float, complex<float> >("nmagr(complex)", TOL_RCP,
    Zc, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::nmagr(a, n, o); },
    [](int k)
    {
      const cdbl z = toCdbl(gZc[k]);
      return (std::norm(z) >= FLT_MIN) ? z/std::abs(z) : cdbl(0.0);
    });

  check<complex<float>, float, complex<float> >("nmagnr(complex)", TOL_NR,
    Zc, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::nmagnr(a, n, o); },
    [](int k)
    {
      const cdbl z = toCdbl(gZc[k]);
      return (std::norm(z) >= FLT_MIN) ? z/std::abs(z) : cdbl(0.0);
    });

  //--- magSum / magrSum / magSqSum : out = w1*|in|(^2) + w2*out --------------
  // in2 holds the values out starts with (copied in by the lambda, so "in
  // place on in2" is the plain case again); the error is measured against
  // |w1|*|in|(^2) + |w2*out| since the sum can cancel.  out may not be on top
  // of the complex input, so that mode is left out
  printf("\nmagSum / magrSum / magSqSum\n");
  {
    const float w1 = gSf, w2 = 1.0f - 0.5f*gSf;
    check<complex<float>, float, float>("magSum(complex)", TOL_EXACT, Zc, Af,
      [w1, w2](const complex<float>* a, const float* b, int n, float* o)
      {
        for(int k=0; k<n; ++k) o[k] = b[k];
        EVM::magSum(a, n, o, w1, w2);
      },
      [w1, w2](int k)
      {
        return cdbl(w1*std::abs(toCdbl(gZc[k])) +
                    w2*static_cast<double>(gAf[k]));
      },
      [w1, w2](int k)
      {
        return std::fabs(w1)*std::abs(toCdbl(gZc[k])) +
               std::fabs(w2*static_cast<double>(gAf[k]));
      }, false);

    check<complex<float>, float, float>("magrSum(complex)", TOL_RCP, Zc, Af,
      [w1, w2](const complex<float>* a, const float* b, int n, float* o)
      {
        for(int k=0; k<n; ++k) o[k] = b[k];
        EVM::magrSum(a, n, o, w1, w2);
      },
      [w1, w2](int k)
      {
        return cdbl(w1*std::abs(toCdbl(gZc[k])) +
                    w2*static_cast<double>(gAf[k]));
      },
      [w1, w2](int k)
      {
        return std::fabs(w1)*std::abs(toCdbl(gZc[k])) +
               std::fabs(w2*static_cast<double>(gAf[k]));
      }, false);

    check<complex<float>, float, float>("magSqSum(complex)", TOL_EXACT, Zc, Af,
      [w1, w2](const complex<float>* a, const float* b, int n, float* o)
      {
        for(int k=0; k<n; ++k) o[k] = b[k];
        EVM::magSqSum(a, n, o, w1, w2);
      },
      [w1, w2](int k)
      {
        return cdbl(w1*std::norm(toCdbl(gZc[k])) +
                    w2*static_cast<double>(gAf[k]));
      },
      [w1, w2](int k)
      {
        return std::fabs(w1)*std::norm(toCdbl(gZc[k])) +
               std::fabs(w2*static_cast<double>(gAf[k]));
      }, false);

  }

  //--- magBias / magrBias : out = |in| + bias --------------------------------
  // the error is measured against |in| + |bias|, since a negative bias can
  // cancel the magnitude
  printf("\nmagBias / magrBias\n");
  check<complex<float>, float, float>("magBias(complex)", TOL_EXACT, Zc, nullF,
    [](const complex<float>* a, const float*, int n, float* o)
      { EVM::magBias(a, n, o, gSf); },
    [](int k){ return cdbl(std::abs(toCdbl(gZc[k])) + gSf); },
    [](int k){ return std::abs(toCdbl(gZc[k])) + std::fabs(gSf); });

  check<complex<float>, float, float>("magrBias(complex)", TOL_RCP, Zc, nullF,
    [](const complex<float>* a, const float*, int n, float* o)
      { EVM::magrBias(a, n, o, gSf); },
    [](int k){ return cdbl(std::abs(toCdbl(gZc[k])) + gSf); },
    [](int k){ return std::abs(toCdbl(gZc[k])) + std::fabs(gSf); });

  //--- magSqMult / magSqMults : out = |in1|^2 * in2 (* scale) ----------------
  printf("\nmagSqMult / magSqMults\n");
  check<complex<float>, float, float>("magSqMult(complex,float)", TOL_EXACT,
    Ac, Bf,
    [](const complex<float>* a, const float* b, int n, float* o)
      { EVM::magSqMult(a, b, n, o); },
    [](int k)
    {
      return cdbl(std::norm(toCdbl(gAc[k]))) * toCdbl(gBf[k]);
    });

  check<float, float, float>("magSqMult(float,float)", TOL_EXACT, Af, Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::magSqMult(a, b, n, o); },
    [](int k)
    {
      return toCdbl(gAf[k]) * toCdbl(gAf[k]) * toCdbl(gBf[k]);
    });

  check<complex<float>, float, float>("magSqMults(complex,float)", TOL_EXACT,
    Ac, Bf,
    [](const complex<float>* a, const float* b, int n, float* o)
      { EVM::magSqMults(a, b, n, gSf, o); },
    [](int k)
    {
      return cdbl(std::norm(toCdbl(gAc[k]))) * toCdbl(gBf[k]) *
             toCdbl(gSf);
    });

  check<float, float, float>("magSqMults(float,float)", TOL_EXACT, Af, Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::magSqMults(a, b, n, gSf, o); },
    [](int k)
    {
      return toCdbl(gAf[k]) * toCdbl(gAf[k]) * toCdbl(gBf[k]) *
             toCdbl(gSf);
    });

  //--- magSqDiv / magSqDivr : out = |in1|^2 / max(in2, eps) ------------------
  printf("\nmagSqDiv / magSqDivr\n");
  check<complex<float>, float, float>("magSqDiv(complex,float)", TOL_EXACT,
    Ac, Bf,
    [](const complex<float>* a, const float* b, int n, float* o)
      { EVM::magSqDiv(a, b, n, o); },
    [](int k)
    {
      return cdbl(std::norm(toCdbl(gAc[k]))) /
             toCdbl(gBf[k]);
    });

  check<float, float, float>("magSqDiv(float,float)", TOL_EXACT, Af, Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::magSqDiv(a, b, n, o); },
    [](int k)
    {
      return toCdbl(gAf[k]) * toCdbl(gAf[k]) /
             toCdbl(gBf[k]);
    });

  check<complex<float>, float, float>("magSqDiv(complex,float) eps", TOL_EXACT,
    Ac, Bf,
    [](const complex<float>* a, const float* b, int n, float* o)
      { EVM::magSqDiv(a, b, n, o, GEPS); },
    [](int k)
    {
      return cdbl(std::norm(toCdbl(gAc[k]))) /
             toCdbl(std::max(gBf[k], GEPS));
    });

  check<float, float, float>("magSqDiv(float,float) eps", TOL_EXACT, Af, Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::magSqDiv(a, b, n, o, GEPS); },
    [](int k)
    {
      return toCdbl(gAf[k]) * toCdbl(gAf[k]) /
             toCdbl(std::max(gBf[k], GEPS));
    });

  check<complex<float>, float, float>("magSqDivr(complex,float)", TOL_RCP,
    Ac, Bf,
    [](const complex<float>* a, const float* b, int n, float* o)
      { EVM::magSqDivr(a, b, n, o); },
    [](int k)
    {
      return cdbl(std::norm(toCdbl(gAc[k]))) /
             toCdbl(gBf[k]);
    });

  check<float, float, float>("magSqDivr(float,float)", TOL_RCP, Af, Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::magSqDivr(a, b, n, o); },
    [](int k)
    {
      return toCdbl(gAf[k]) * toCdbl(gAf[k]) /
             toCdbl(gBf[k]);
    });

  check<complex<float>, float, float>("magSqDivr(complex,float) eps", TOL_RCP,
    Ac, Bf,
    [](const complex<float>* a, const float* b, int n, float* o)
      { EVM::magSqDivr(a, b, n, o, GEPS); },
    [](int k)
    {
      return cdbl(std::norm(toCdbl(gAc[k]))) /
             toCdbl(std::max(gBf[k], GEPS));
    });

  check<float, float, float>("magSqDivr(float,float) eps", TOL_RCP, Af, Bf,
    [](const float* a, const float* b, int n, float* o)
      { EVM::magSqDivr(a, b, n, o, GEPS); },
    [](int k)
    {
      return toCdbl(gAf[k]) * toCdbl(gAf[k]) /
             toCdbl(std::max(gBf[k], GEPS));
    });

  //--- magSqDiff : out = |in1 - in2|^2 ---------------------------------------
  // the second run swaps the inputs, so its in-place mode puts out on top of
  // in2.  The error is measured against |in1|^2 + |in2|^2, since the
  // difference can cancel
  printf("\nmagSqDiff\n");
  check<complex<float>, complex<float>, float>("magSqDiff(complex,complex)",
    TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, float* o)
      { EVM::magSqDiff(a, b, n, o); },
    [](int k){ return cdbl(std::norm(toCdbl(gAc[k]) - toCdbl(gBc[k]))); },
    [](int k)
      { return std::norm(toCdbl(gAc[k])) + std::norm(toCdbl(gBc[k])); });

  check<complex<float>, complex<float>, float>(
    "magSqDiff(complex,complex) a=in2", TOL_EXACT, Ac, Bc,
    [](const complex<float>* a, const complex<float>* b, int n, float* o)
      { EVM::magSqDiff(b, a, n, o); },
    [](int k){ return cdbl(std::norm(toCdbl(gAc[k]) - toCdbl(gBc[k]))); },
    [](int k)
      { return std::norm(toCdbl(gAc[k])) + std::norm(toCdbl(gBc[k])); });

  //--- magnitude special values -------------------------------------------
  // icpx's fast math turns 1.0f/sqrtf(x) in the default versions into an
  // rsqrt approximation and a Newton-Raphson step, which gives NaN for an
  // overflowing |in|^2, so the normalizations are only checked on the SIMD
  // paths with icpx.  magr's AVX-512 version may give NaN for an overflowing
  // |in|^2 (the other paths give INF)
  printf("\nmagnitude special values\n");
  {
    const std::string lvl = EVM::getOptLevel();
#if defined(__INTEL_LLVM_COMPILER)
    const bool normSpecial = (lvl == "AVX512" || lvl == "AVX2" ||
                              lvl == "AVX");
#else
    const bool normSpecial = true;
#endif
    checkMagSpecial("magSq(complex) special", TOL_EXACT, MS_SQ,
      [](const complex<float>* a, int n, float* o, complex<float>* c)
        { (void)o; (void)c; EVM::magSq(a, n, o); });
    checkMagSpecial("mag(complex) special", TOL_EXACT, MS_MAG,
      [](const complex<float>* a, int n, float* o, complex<float>* c)
        { (void)o; (void)c; EVM::mag(a, n, o); });
    checkMagSpecial("magr(complex) special", TOL_RCP, MS_MAG,
      [](const complex<float>* a, int n, float* o, complex<float>* c)
        { (void)o; (void)c; EVM::magr(a, n, o); }, lvl == "AVX512");
    checkMagSpecial("magBias(complex) special", TOL_EXACT, MS_BIAS,
      [](const complex<float>* a, int n, float* o, complex<float>* c)
        { (void)o; (void)c; EVM::magBias(a, n, o, 0.5f); });
    checkMagSpecial("magrBias(complex) special", TOL_RCP, MS_BIAS,
      [](const complex<float>* a, int n, float* o, complex<float>* c)
        { (void)o; (void)c; EVM::magrBias(a, n, o, 0.5f); });
    checkMagSpecial("magSum(complex) special", TOL_EXACT, MS_SUM,
      [](const complex<float>* a, int n, float* o, complex<float>* c)
        { (void)o; (void)c; EVM::magSum(a, n, o, 0.75f, 0.25f); });
    checkMagSpecial("magrSum(complex) special", TOL_RCP, MS_SUM,
      [](const complex<float>* a, int n, float* o, complex<float>* c)
        { (void)o; (void)c; EVM::magrSum(a, n, o, 0.75f, 0.25f); });
    if(normSpecial)
    {
      checkMagSpecial("nmag(complex) special", TOL_EXACT, MS_NORM,
        [](const complex<float>* a, int n, float* o, complex<float>* c)
          { (void)o; (void)c; EVM::nmag(a, n, c); });
      checkMagSpecial("nmagr(complex) special", TOL_RCP, MS_NORM,
        [](const complex<float>* a, int n, float* o, complex<float>* c)
          { (void)o; (void)c; EVM::nmagr(a, n, c); });
      checkMagSpecial("nmagnr(complex) special", TOL_NR, MS_NORM,
        [](const complex<float>* a, int n, float* o, complex<float>* c)
          { (void)o; (void)c; EVM::nmagnr(a, n, c); });
    }
    else
      printf("  %-34s skipped (icpx default versions)\n",
             "nmag ... nmagnr special");
  }

  //--- sum, sumSqDiff : sums of a buffer ------------------------------------
  // sum(float) and sum(complex) add the elements, sumSqDiff |in1-in2|^2
  printf("\nsum / sumSqDiff\n");
  {
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float inf = std::numeric_limits<float>::infinity();
    checkSum<float>("sum(float)", Af, Bf,
      [](const float* a, const float*, int n, float* r)
        { EVM::sum(a, n, r[0]); },
      [](const float x, const float){ return toCdbl(x); }, 1, 0,
      [=](const int pass, cdbl& want)
        {
          want = (pass == 0) ? cdbl(nan, 0.0) : cdbl(-inf, 0.0);
          return (pass == 0) ? nan : -inf;
        });
    checkSum<complex<float> >("sum(complex)", Ac, Bc,
      [](const complex<float>* a, const complex<float>*, int n,
         complex<float>* r) { EVM::sum(a, n, r[0]); },
      [](const complex<float> x, const complex<float>){ return toCdbl(x); },
      1, 0,
      [=](const int pass, cdbl& want)
        {
          want = (pass == 0) ? cdbl(nan, inf) : cdbl(-inf, nan);
          return (pass == 0) ? complex<float>(nan, inf) :
                               complex<float>(-inf, nan);
        });
    checkSum<float>("sumSqDiff(complex,complex)", Ac, Bc,
      [](const complex<float>* a, const complex<float>* b, int n, float* r)
        { EVM::sumSqDiff(a, b, n, r[0]); },
      [](const complex<float> x, const complex<float> y)
        { return cdbl(std::norm(toCdbl(x) - toCdbl(y))); }, 2, 3,
      [=](const int pass, cdbl& want)
        {
          want = (pass == 0) ? cdbl(nan, 0.0) : cdbl(inf, 0.0);
          return (pass == 0) ? complex<float>(nan, 0.5f) :
                               complex<float>(1.0f, -inf);
        });
  }

  //--- wsum : out = in1 + w*in2 and out = w1*in1 + w2*in2 -------------------
  // the first weights are gSf and gSc, the second ones W2F and W2C.  The
  // error is measured against |in1| + |w*in2| (|w1*in1| + |w2*in2|), since
  // the terms can cancel, as for multSum.  The complex overloads with real
  // weights run the real versions on twice the length
  printf("\nwsum\n");
  {
    static const float W2F = -0.73f;
    static const complex<float> W2C(0.37f, -1.1f);
    check<float, float, float>("wsum(float,float,w)", TOL_EXACT, Af, Bf,
      [](const float* a, const float* b, int n, float* o)
        { EVM::wsum(a, b, n, gSf, o); },
      [](int k){ return toCdbl(gAf[k]) +
                        static_cast<double>(gSf)*toCdbl(gBf[k]); },
      [](int k){ return std::abs(toCdbl(gAf[k])) +
                        std::abs(static_cast<double>(gSf)*toCdbl(gBf[k])); });

    check<complex<float>, complex<float>, complex<float> >(
      "wsum(complex,complex,w)", TOL_EXACT, Ac, Bc,
      [](const complex<float>* a, const complex<float>* b, int n,
         complex<float>* o) { EVM::wsum(a, b, n, gSf, o); },
      [](int k){ return toCdbl(gAc[k]) +
                        static_cast<double>(gSf)*toCdbl(gBc[k]); },
      [](int k){ return std::abs(toCdbl(gAc[k])) +
                        std::abs(static_cast<double>(gSf)*toCdbl(gBc[k])); });

    check<complex<float>, complex<float>, complex<float> >(
      "wsum(complex,complex,cw)", TOL_EXACT, Ac, Bc,
      [](const complex<float>* a, const complex<float>* b, int n,
         complex<float>* o) { EVM::wsum(a, b, n, gSc, o); },
      [](int k){ return toCdbl(gAc[k]) + toCdbl(gSc)*toCdbl(gBc[k]); },
      [](int k){ return std::abs(toCdbl(gAc[k])) +
                        std::abs(toCdbl(gSc)*toCdbl(gBc[k])); });

    check<float, float, float>("wsum(float,float,w1,w2)", TOL_EXACT, Af, Bf,
      [](const float* a, const float* b, int n, float* o)
        { EVM::wsum(a, b, n, gSf, W2F, o); },
      [](int k){ return static_cast<double>(gSf)*toCdbl(gAf[k]) +
                        static_cast<double>(W2F)*toCdbl(gBf[k]); },
      [](int k){ return std::abs(static_cast<double>(gSf)*toCdbl(gAf[k])) +
                        std::abs(static_cast<double>(W2F)*toCdbl(gBf[k])); });

    check<complex<float>, complex<float>, complex<float> >(
      "wsum(complex,complex,w1,w2)", TOL_EXACT, Ac, Bc,
      [](const complex<float>* a, const complex<float>* b, int n,
         complex<float>* o) { EVM::wsum(a, b, n, gSf, W2F, o); },
      [](int k){ return static_cast<double>(gSf)*toCdbl(gAc[k]) +
                        static_cast<double>(W2F)*toCdbl(gBc[k]); },
      [](int k){ return std::abs(static_cast<double>(gSf)*toCdbl(gAc[k])) +
                        std::abs(static_cast<double>(W2F)*toCdbl(gBc[k])); });

    check<complex<float>, complex<float>, complex<float> >(
      "wsum(complex,complex,cw1,cw2)", TOL_EXACT, Ac, Bc,
      [](const complex<float>* a, const complex<float>* b, int n,
         complex<float>* o) { EVM::wsum(a, b, n, gSc, W2C, o); },
      [](int k){ return toCdbl(gSc)*toCdbl(gAc[k]) +
                        toCdbl(W2C)*toCdbl(gBc[k]); },
      [](int k){ return std::abs(toCdbl(gSc)*toCdbl(gAc[k])) +
                        std::abs(toCdbl(W2C)*toCdbl(gBc[k])); });
  }

  //--- cnvrtInt*ToFloat : out = static_cast<float>(in) ----------------------
  // the conversions must give exactly the float the cast gives (int8/int16
  // are exact, int32 rounds to nearest); "in place" goes through the void*
  // wrapper, whose output is wider than its input
  printf("\ninteger to float conversion\n");
  check<int8_t, float, float>("cnvrtInt8ToFloat", TOL_CMP, gI8.data(), nullF,
    [](const int8_t* a, const float*, int n, float* o)
      { EVM::cnvrtInt8ToFloat(a, n, o); },
    [](int k){ return cdbl(static_cast<float>(gI8[k])); });

  check<int8_t, float, float>("cnvrtInt8ToFloat in place", TOL_CMP,
    gI8.data(), nullF,
    [](const int8_t* a, const float*, int n, float* o)
      { runCnvrtInPlace(a, n, o, [](void* p, int m)
          { EVM::cnvrtInt8ToFloat(p, m); }); },
    [](int k){ return cdbl(static_cast<float>(gI8[k])); });

  checkCnvrtMinMax("cnvrtInt8ToFloat with min/max", gI8.data(),
    [](const int8_t* a, int n, float* o, float& mn, float& mx, bool ni)
      { EVM::cnvrtInt8ToFloat(a, n, o, mn, mx, ni); });

  check<int16_t, float, float>("cnvrtInt16ToFloat", TOL_CMP, gI16.data(),
    nullF,
    [](const int16_t* a, const float*, int n, float* o)
      { EVM::cnvrtInt16ToFloat(a, n, o); },
    [](int k){ return cdbl(static_cast<float>(gI16[k])); });

  check<int16_t, float, float>("cnvrtInt16ToFloat in place", TOL_CMP,
    gI16.data(), nullF,
    [](const int16_t* a, const float*, int n, float* o)
      { runCnvrtInPlace(a, n, o, [](void* p, int m)
          { EVM::cnvrtInt16ToFloat(p, m); }); },
    [](int k){ return cdbl(static_cast<float>(gI16[k])); });

  checkCnvrtMinMax("cnvrtInt16ToFloat with min/max", gI16.data(),
    [](const int16_t* a, int n, float* o, float& mn, float& mx, bool ni)
      { EVM::cnvrtInt16ToFloat(a, n, o, mn, mx, ni); });

  check<int32_t, float, float>("cnvrtInt32ToFloat", TOL_CMP, gI32.data(),
    nullF,
    [](const int32_t* a, const float*, int n, float* o)
      { EVM::cnvrtInt32ToFloat(a, n, o); },
    [](int k){ return cdbl(static_cast<float>(gI32[k])); });

  check<int32_t, float, float>("cnvrtInt32ToFloat in place", TOL_CMP,
    gI32.data(), nullF,
    [](const int32_t* a, const float*, int n, float* o)
      { runCnvrtInPlace(a, n, o, [](void* p, int m)
          { EVM::cnvrtInt32ToFloat(p, m); }); },
    [](int k){ return cdbl(static_cast<float>(gI32[k])); });

  checkCnvrtMinMax("cnvrtInt32ToFloat with min/max", gI32.data(),
    [](const int32_t* a, int n, float* o, float& mn, float& mx, bool ni)
      { EVM::cnvrtInt32ToFloat(a, n, o, mn, mx, ni); });

  //--- polToCart / cartToPol : polar (mag, phase) <-> cartesian -------------
  // the polar inputs are gAc as (magnitude, phase), the cartesian inputs are
  // gZc so that atan2(0, 0) is covered.  The cartesian error is measured
  // against the magnitude, since one part can be near zero.  polToCart,
  // cartToPol and cartToPolRad always use evm_trig.h (in their SIMD
  // versions), polToCartRad uses SVML on the AVX-512 and AVX2 paths when
  // EVM_SVML is defined (icpx without -DDISABLE_SVML) and evm_trig.h
  // otherwise
  printf("\npolar / cartesian conversion  (sin/cos of radians: %s)\n",
#if defined(EVM_SVML)
         "SVML, evm_trig.h on AVX");
#else
         "evm_trig.h");
#endif
  check<complex<float>, float, complex<float> >("polToCart(complex)",
    TOL_TRIG, Ac, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::polToCart(a, n, o); },
    [](int k){ return polarRef(gAc[k], 2.0*M_PI); },
    [](int k){ return std::fabs(static_cast<double>(gAc[k].real())); });

  check<complex<float>, float, complex<float> >("polToCartRad(complex)",
    TOL_TRIG, Ac, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::polToCartRad(a, n, o); },
    [](int k){ return polarRef(gAc[k], 1.0); },
    [](int k){ return std::fabs(static_cast<double>(gAc[k].real())); });

  check<complex<float>, float, complex<float> >("polToCart(complex) large phase",
    TOL_TRIG, gLc.data(), nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::polToCart(a, n, o); },
    [](int k){ return polarRefCyc(gLc[k]); },
    [](int k){ return std::fabs(static_cast<double>(gLc[k].real())); });

  check<complex<float>, float, complex<float> >("polToCartRad(complex) large x",
    TOL_TRIG, gLr.data(), nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::polToCartRad(a, n, o); },
    [](int k){ return polarRef(gLr[k], 1.0); },
    [](int k){ return std::fabs(static_cast<double>(gLr[k].real())); });

  // polToCarts: the phases (cycles) times the scale, rounded to float as the
  // kernels round them (volatile, so the product is not kept in double),
  // then whole cycles removed exactly
  auto scaledRef = [](const complex<float> v, const float sc)
  {
    volatile float ph = sc*v.imag();
    return polarRefCyc(complex<float>(v.real(), ph));
  };
  check<complex<float>, float, complex<float> >("polToCarts(complex)",
    TOL_TRIG, Ac, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::polToCarts(a, n, 0.37f, o); },
    [&](int k){ return scaledRef(gAc[k], 0.37f); },
    [](int k){ return std::fabs(static_cast<double>(gAc[k].real())); });
  check<complex<float>, float, complex<float> >("polToCarts(complex) large",
    TOL_TRIG, gLc.data(), nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::polToCarts(a, n, -2.5f, o); },
    [&](int k){ return scaledRef(gLc[k], -2.5f); },
    [](int k){ return std::fabs(static_cast<double>(gLc[k].real())); });

  // the overloads with the magnitudes and phases in separate buffers (gAc's,
  // and gLc's and gLr's large phases, whose magnitudes are gAc's)
  {
    std::vector<float> am(MAXLEN), ap(MAXLEN), lc(MAXLEN), lr(MAXLEN);
    for(int i=0; i<MAXLEN; ++i)
    {
      am[i] = gAc[i].real();
      ap[i] = gAc[i].imag();
      lc[i] = gLc[i].imag();
      lr[i] = gLr[i].imag();
    }
    checkPolarSplit("polToCart(mag, phase)", am.data(), ap.data(),
      [](const float* m, const float* ph, int n, complex<float>* o)
        { EVM::polToCart(m, ph, n, o); },
      [](int k){ return polarRef(gAc[k], 2.0*M_PI); });
    checkPolarSplit("polToCart(mag, phase) large phase", am.data(),
      lc.data(),
      [](const float* m, const float* ph, int n, complex<float>* o)
        { EVM::polToCart(m, ph, n, o); },
      [](int k){ return polarRefCyc(gLc[k]); });
    checkPolarSplit("polToCartRad(mag, phase)", am.data(), ap.data(),
      [](const float* m, const float* ph, int n, complex<float>* o)
        { EVM::polToCartRad(m, ph, n, o); },
      [](int k){ return polarRef(gAc[k], 1.0); });
    checkPolarSplit("polToCartRad(mag, phase) large x", am.data(), lr.data(),
      [](const float* m, const float* ph, int n, complex<float>* o)
        { EVM::polToCartRad(m, ph, n, o); },
      [](int k){ return polarRef(gLr[k], 1.0); });
    checkPolarSplit("polToCarts(mag, phase)", am.data(), ap.data(),
      [](const float* m, const float* ph, int n, complex<float>* o)
        { EVM::polToCarts(m, ph, n, 0.37f, o); },
      [&](int k){ return scaledRef(gAc[k], 0.37f); });
    checkPolarSplit("polToCarts(mag, phase) large", am.data(), lc.data(),
      [](const float* m, const float* ph, int n, complex<float>* o)
        { EVM::polToCarts(m, ph, n, -2.5f, o); },
      [&](int k){ return scaledRef(gLc[k], -2.5f); });
  }

  check<complex<float>, float, complex<float> >("cartToPol(complex)",
    TOL_TRIG, Zc, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::cartToPol(a, n, o); },
    [](int k){ return cdbl(std::abs(toCdbl(gZc[k])),
                           std::arg(toCdbl(gZc[k]))/(2.0*M_PI)); });

  check<complex<float>, float, complex<float> >("cartToPolRad(complex)",
    TOL_TRIG, Zc, nullF,
    [](const complex<float>* a, const float*, int n, complex<float>* o)
      { EVM::cartToPolRad(a, n, o); },
    [](int k){ return cdbl(std::abs(toCdbl(gZc[k])),
                           std::arg(toCdbl(gZc[k]))); });

  // special values: polar inputs (magnitude, phase) and cartesian (re, im).
  // The magnitudes leave out infinities: an infinite magnitude at a quarter
  // cycle is inf*0 = NaN where the sine or cosine is exactly 0 (evm_trig.h)
  // but +-inf where it is only near 0 (sincosf), and neither is meaningful.
  // Phases stay below 8.5e37 cycles, where the AVX2/AVX polToCart overflow
  // to NaN
  {
    const float inf = std::numeric_limits<float>::infinity();
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float mag[] = {0.0f, -0.0f, 1.0f, -1.0f, 0.5f, nan, 1.0e-40f,
                         3.4e38f, -3.4e38f};
    const float ph[] = {0.0f, -0.0f, 0.25f, -0.5f, 1.0f, -1.0f, 0.125f, inf,
                        -inf, nan, 1.0e-40f, 1.0e37f, -1.0e30f};
    const float cart[] = {0.0f, -0.0f, 1.0f, -1.0f, inf, -inf, nan, 1.0e-40f,
                          -1.0e-40f, 3.4e38f, -3.4e38f, 0.5f};
    const int nm = sizeof(mag)/sizeof(mag[0]);
    const int np = sizeof(ph)/sizeof(ph[0]);
    const int nc = sizeof(cart)/sizeof(cart[0]);
    // magnitude times the sine or cosine; an infinite or NaN phase gives NaN.
    // Phases in cycles are reduced exactly to the nearest quarter cycle first,
    // so quarter cycles give exact zeros
    auto polWant = [](const complex<float> v, const int part, const bool cyc)
    {
      const double m = v.real(), p = v.imag();
      if(!std::isfinite(p)) return std::nan("");
      if(!cyc) return m*(part ? std::sin(p) : std::cos(p));
      const double r = p - std::nearbyint(p);
      const double q = std::nearbyint(4.0*r), t = r - q/4.0;
      const double s0 = std::sin(2.0*M_PI*t), c0 = std::cos(2.0*M_PI*t);
      const int quad = (static_cast<int>(q)%4 + 4)%4;
      const double sv[4] = {s0, c0, -s0, -c0}, cv[4] = {c0, -s0, -c0, s0};
      return m*(part ? sv[quad] : cv[quad]);
    };
    auto polSize = [](const complex<float> v, const int)
      { return std::fabs(static_cast<double>(v.real())); };
    // |z| as the kernels compute it (|z|^2 overflows to inf past 1.8e19)
    auto cartWant = [](const complex<float> v, const int part, const bool cyc)
    {
      const double x = v.real(), y = v.imag();
      if(part == 0)
        return std::sqrt(static_cast<double>(static_cast<float>(x*x + y*y)));
      return std::atan2(y, x)/(cyc ? 2.0*M_PI : 1.0);
    };
    auto cartSize = [](const complex<float> v, const int part)
    {
      const double x = v.real(), y = v.imag();
      return part ? std::fabs(std::atan2(y, x))
                  : std::sqrt(static_cast<double>(static_cast<float>(x*x + y*y)));
    };
    checkPolarSpecial("polToCart special values", mag, nm, ph, np,
      [](const complex<float>* a, int n, complex<float>* o)
        { EVM::polToCart(a, n, o); },
      [&](const complex<float> v, int part){ return polWant(v, part, true); },
      polSize, false);
    checkPolarSpecial("polToCartRad special values", mag, nm, ph, np,
      [](const complex<float>* a, int n, complex<float>* o)
        { EVM::polToCartRad(a, n, o); },
      [&](const complex<float> v, int part){ return polWant(v, part, false); },
      polSize, false);
    // the split overloads get the same grid, split into magnitudes and
    // phases first (so the in-place run only covers the complex ones), and
    // the scaled ones a scale of -1.25, the phase times the scale rounded
    // to float first as in the kernels
    auto split = [](const complex<float>* a, const int n,
                    std::vector<float>& m, std::vector<float>& p)
    {
      m.resize(n); p.resize(n);
      for(int k=0; k<n; ++k) { m[k] = a[k].real(); p[k] = a[k].imag(); }
    };
    auto scWant = [&](const complex<float> v, const int part)
    {
      volatile float p = -1.25f*v.imag();
      return polWant(complex<float>(v.real(), p), part, true);
    };
    checkPolarSpecial("polToCart(mag, phase) special values", mag, nm, ph,
      np, [&](const complex<float>* a, int n, complex<float>* o)
        { std::vector<float> m, p; split(a, n, m, p);
          EVM::polToCart(m.data(), p.data(), n, o); },
      [&](const complex<float> v, int part){ return polWant(v, part, true); },
      polSize, false);
    checkPolarSpecial("polToCartRad(mag, phase) special", mag, nm, ph, np,
      [&](const complex<float>* a, int n, complex<float>* o)
        { std::vector<float> m, p; split(a, n, m, p);
          EVM::polToCartRad(m.data(), p.data(), n, o); },
      [&](const complex<float> v, int part){ return polWant(v, part, false); },
      polSize, false);
    checkPolarSpecial("polToCarts(complex) special values", mag, nm, ph, np,
      [](const complex<float>* a, int n, complex<float>* o)
        { EVM::polToCarts(a, n, -1.25f, o); },
      scWant, polSize, false);
    checkPolarSpecial("polToCarts(mag, phase) special", mag, nm, ph, np,
      [&](const complex<float>* a, int n, complex<float>* o)
        { std::vector<float> m, p; split(a, n, m, p);
          EVM::polToCarts(m.data(), p.data(), n, -1.25f, o); },
      scWant, polSize, false);
    checkPolarSpecial("cartToPol special values", cart, nc, cart, nc,
      [](const complex<float>* a, int n, complex<float>* o)
        { EVM::cartToPol(a, n, o); },
      [&](const complex<float> v, int part){ return cartWant(v, part, true); },
      cartSize, true);
    checkPolarSpecial("cartToPolRad special values", cart, nc, cart, nc,
      [](const complex<float>* a, int n, complex<float>* o)
        { EVM::cartToPolRad(a, n, o); },
      [&](const complex<float> v, int part){ return cartWant(v, part, false); },
      cartSize, true);
  }

  //--- exp / log ------------------------------------------------------------
  // exp runs on gEf ([-88, 88], results from denormals to near FLT_MAX), log
  // on gLf (FLT_MIN to 2^125, a quarter of them near 1).  Both use SVML only
  // on the AVX path when EVM_SVML is defined, and evm_trig.h otherwise.
  // TOL_EXPLOG passes SVML's
  // worst case with icpx's default fast fp-model (4 ulp) but not a wrong
  // polynomial coefficient that TOL_TRIG would let through
  printf("\nexp / log  (%s)\n",
#if defined(EVM_SVML)
         "SVML on AVX, evm_trig.h on AVX-512/AVX2");
#else
         "evm_trig.h");
#endif
  check<float, float, float>("exp(float)", TOL_EXPLOG, gEf.data(), nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::exp(a, n, o); },
    [](int k){ return cdbl(std::exp(static_cast<double>(gEf[k]))); });

  check<float, float, float>("log(float)", TOL_EXPLOG, gLf.data(), nullF,
    [](const float* a, const float*, int n, float* o)
      { EVM::log(a, n, o); },
    [](int k){ return cdbl(std::log(static_cast<double>(gLf[k]))); });

  // icpx's fast math in the default version gives NaN for log of a negative
  // denormal under DAZ when it runs in place, where its vectorized loop and
  // the SIMD paths read the input as -0 and give -INF
  {
    const std::string lvl = EVM::getOptLevel();
#if defined(__INTEL_LLVM_COMPILER)
    const bool negDenNaN = !(lvl == "AVX512" || lvl == "AVX2" ||
                             lvl == "AVX");
#else
    const bool negDenNaN = false;
#endif
    checkExpLogSpecial("exp(float) special", false, false,
      [](const float* a, int n, float* o){ EVM::exp(a, n, o); });
    checkExpLogSpecial("log(float) special", true, negDenNaN,
      [](const float* a, int n, float* o){ EVM::log(a, n, o); });
  }

  //--- tune / tunecx / tunesp : out = in * e^(j*(phInit + k*phIncr)) -------
  // tune and tunecx advance a phasor with complex multiplies, so their error
  // may ramp up along the buffer, by 1e-8 per element on the SIMD paths
  // (whose phasors round once per 32 elements, 16 on AVX) and 2e-7 in the
  // default versions (once per element, from a float phase increment).
  // tunesp's phasors come from double precision phases, so its error does
  // not ramp up.  The phases cover an increment near pi radians (or half a
  // cycle) and a large phInit
  printf("\nNCO tuning\n");
  {
    typedef complex<float> C;
    const std::string lvl = EVM::getOptLevel();
    const bool simd = (lvl == "AVX512" || lvl == "AVX2" || lvl == "AVX");
    const double drift = simd ? 1.0e-8 : 2.0e-7;
    const std::vector<double> rad = {0.7, 0.37, -2.9, 3.1, 1000.0, -1.0e-3};
    const std::vector<double> cyc = {0.3, 0.0123, -0.45, 0.49, 123.25, 0.4321};
    checkTune<C>("tune(complex)", false, false, true, drift, rad,
      [](const C* a, int n, double p, double d, C* o, C*)
        { return EVM::tune(a, n, p, d, o); });
    checkTune<float>("tune(real)", false, false, false, drift, rad,
      [](const float* a, int n, double p, double d, C* o, C*)
        { return EVM::tune(a, n, p, d, o); });
    checkTune<C>("tune(complex) phOut", false, true, true, drift, rad,
      [](const C* a, int n, double p, double d, C* o, C* ph)
        { return EVM::tune(a, n, p, d, o, ph); });
    checkTune<float>("tune(real) phOut", false, true, false, drift, rad,
      [](const float* a, int n, double p, double d, C* o, C* ph)
        { return EVM::tune(a, n, p, d, o, ph); });
    checkTune<C>("tunecx(complex)", true, false, true, drift, cyc,
      [](const C* a, int n, double p, double d, C* o, C*)
        { return EVM::tunecx(a, n, p, d, o); });
    checkTune<float>("tunecx(real)", true, false, false, drift, cyc,
      [](const float* a, int n, double p, double d, C* o, C*)
        { return EVM::tunecx(a, n, p, d, o); });
    checkTune<C>("tunecx(complex) phOut", true, true, true, drift, cyc,
      [](const C* a, int n, double p, double d, C* o, C* ph)
        { return EVM::tunecx(a, n, p, d, o, ph); });
    checkTune<float>("tunecx(real) phOut", true, true, false, drift, cyc,
      [](const float* a, int n, double p, double d, C* o, C* ph)
        { return EVM::tunecx(a, n, p, d, o, ph); });
    checkTune<C>("tunesp(complex)", true, false, false, 0.0, cyc,
      [](const C* a, int n, double p, double d, C* o, C*)
        { return EVM::tunesp(a, n, p, d, o); });
    checkTune<float>("tunesp(real)", true, false, false, 0.0, cyc,
      [](const float* a, int n, double p, double d, C* o, C*)
        { return EVM::tunesp(a, n, p, d, o); });
    checkTune<C>("tunesp(complex) phOut", true, true, false, 0.0, cyc,
      [](const C* a, int n, double p, double d, C* o, C* ph)
        { return EVM::tunesp(a, n, p, d, o, ph); });
    checkTune<float>("tunesp(real) phOut", true, true, false, 0.0, cyc,
      [](const float* a, int n, double p, double d, C* o, C* ph)
        { return EVM::tunesp(a, n, p, d, o, ph); });
  }

  //--- cxexp / cxexpsi / cxexpso / cxexpss : os*e^(j*(is*in + eps)) ---------
  printf("\ncomplex exponentials  (sin/cos of radians: %s)\n",
#if defined(EVM_SVML)
         "SVML, evm_trig.h on AVX");
#else
         "evm_trig.h");
#endif
  {
    typedef complex<float> C;
    std::vector<float> xs(MAXLEN), xl(MAXLEN);
    for(int i=0; i<MAXLEN; ++i)
    {
      xs[i] = 3.0f*gAf[i];        // |x| <= 6
      xl[i] = gLr[i].imag();      // up to 1e7, past the 65536 where
    }                             // evm_trig.h hands lanes to sincosf
    checkCxexp("cxexp", xs, 1.0, 1.0, -0.0f, true,
      [](const float* a, int n, C* o) { EVM::cxexp(a, n, o); });
    checkCxexp("cxexp large x", xl, 1.0, 1.0, -0.0f, true,
      [](const float* a, int n, C* o) { EVM::cxexp(a, n, o); });
    checkCxexp("cxexpsi", xs, 2.5, 1.0, -0.0f, true,
      [](const float* a, int n, C* o) { EVM::cxexpsi(a, n, 2.5f, o); });
    checkCxexp("cxexpsi eps", xs, -0.75, 1.0, -1.25f, false,
      [](const float* a, int n, C* o)
        { EVM::cxexpsi(a, n, -0.75f, o, -1.25f); });
    checkCxexp("cxexpsi eps = 0", xs, 2.5, 1.0, 0.0f, false,
      [](const float* a, int n, C* o) { EVM::cxexpsi(a, n, 2.5f, o, 0.0f); });
    checkCxexp("cxexpso", xs, 1.0, 3.0, -0.0f, true,
      [](const float* a, int n, C* o) { EVM::cxexpso(a, n, 3.0f, o); });
    checkCxexp("cxexpso large x", xl, 1.0, -0.5, -0.0f, true,
      [](const float* a, int n, C* o) { EVM::cxexpso(a, n, -0.5f, o); });
    checkCxexp("cxexpso eps", xs, 1.0, -0.5, 0.3f, false,
      [](const float* a, int n, C* o) { EVM::cxexpso(a, n, -0.5f, o, 0.3f); });
    checkCxexp("cxexpss", xs, 2.5, -0.5, -0.0f, true,
      [](const float* a, int n, C* o) { EVM::cxexpss(a, n, 2.5f, -0.5f, o); });
    checkCxexp("cxexpss eps", xs, -0.75, 3.0, 0.3f, false,
      [](const float* a, int n, C* o)
        { EVM::cxexpss(a, n, -0.75f, 3.0f, o, 0.3f); });
  }

  //--- atan2 / atan2s / multcatan2 / multcatan2s : angles in radians -------
  // atan2 runs on gZc, so that the zeros of all four signs are covered, and
  // multcatan2 on gZc and gCc, a second time with the buffers passed the
  // other way round so that the output also goes on top of in2.  The error
  // of multcatan2 is measured against max(|angle|, 1) times |scale|: the
  // product is rounded in float, so its angle is only within a few 2^-24 rad,
  // which is large next to a small angle.  The SIMD versions always use
  // evm_trig.h's atan2, as cartToPolRad does
  printf("\narctangents  (atan2: evm_trig.h)\n");
  {
    typedef complex<float> C;
    const C* Cc = gCc.data();
    auto arg = [](const C v) { return std::arg(toCdbl(v)); };
    // angle of a*conj(b), from the exact product
    auto argMult = [](const C a, const C b)
    {
      const double ar = a.real(), ai = a.imag();
      const double br = b.real(), bi = b.imag();
      return std::atan2(ai*br - ar*bi, ar*br + ai*bi);
    };
    check<C, float, float>("atan2(complex)", TOL_TRIG, Zc, nullF,
      [](const C* a, const float*, int n, float* o) { EVM::atan2(a, n, o); },
      [&](int k){ return cdbl(arg(gZc[k]), 0.0); });
    check<C, float, C>("atan2(complex) complex out", TOL_TRIG, Zc, nullF,
      [](const C* a, const float*, int n, C* o) { EVM::atan2(a, n, o); },
      [&](int k){ return cdbl(arg(gZc[k]), 0.0); });
    check<C, float, float>("atan2s(complex)", TOL_TRIG, Zc, nullF,
      [](const C* a, const float*, int n, float* o)
        { EVM::atan2s(a, n, o, gSf); },
      [&](int k){ return cdbl(gSf*arg(gZc[k]), 0.0); });
    check<C, float, C>("atan2s(complex) complex out", TOL_TRIG, Zc, nullF,
      [](const C* a, const float*, int n, C* o) { EVM::atan2s(a, n, o, gSf); },
      [&](int k){ return cdbl(gSf*arg(gZc[k]), 0.0); });
    check<C, C, float>("multcatan2", TOL_TRIG, Zc, Cc,
      [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2(a, b, n, o); },
      [&](int k){ return cdbl(argMult(gZc[k], gCc[k]), 0.0); },
      [&](int k){ return std::max(std::fabs(argMult(gZc[k], gCc[k])), 1.0); });
    check<C, C, float>("multcatan2 (output on in2)", TOL_TRIG, Cc, Zc,
      [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2(b, a, n, o); },
      [&](int k){ return cdbl(argMult(gZc[k], gCc[k]), 0.0); },
      [&](int k){ return std::max(std::fabs(argMult(gZc[k], gCc[k])), 1.0); });
    check<C, C, float>("multcatan2s", TOL_TRIG, Zc, Cc,
      [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2s(a, b, n, o, gSf); },
      [&](int k){ return cdbl(gSf*argMult(gZc[k], gCc[k]), 0.0); },
      [&](int k){ return std::fabs(gSf)*
                         std::max(std::fabs(argMult(gZc[k], gCc[k])), 1.0); });
    check<C, C, float>("multcatan2s (output on in2)", TOL_TRIG, Cc, Zc,
      [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2s(b, a, n, o, gSf); },
      [&](int k){ return cdbl(gSf*argMult(gZc[k], gCc[k]), 0.0); },
      [&](int k){ return std::fabs(gSf)*
                         std::max(std::fabs(argMult(gZc[k], gCc[k])), 1.0); });
    // zeroAngle: the angle is +0 wherever the product is zero (every third
    // element of gZc is a zero of one of the four signs, which otherwise
    // gives +-0 or +-pi)
    auto argMultZ = [&](const C a, const C b)
    {
      const double ar = a.real(), ai = a.imag();
      const double br = b.real(), bi = b.imag();
      return (ar*br + ai*bi == 0.0 && ai*br - ar*bi == 0.0) ? 0.0 :
          argMult(a, b);
    };
    check<C, C, float>("multcatan2 zeroAngle", TOL_TRIG, Zc, Cc,
      [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2(a, b, n, o, true); },
      [&](int k){ return cdbl(argMultZ(gZc[k], gCc[k]), 0.0); },
      [&](int k){ return std::max(std::fabs(argMultZ(gZc[k], gCc[k])), 1.0); });
    check<C, C, float>("multcatan2s zeroAngle (out on in2)", TOL_TRIG, Cc, Zc,
      [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2s(b, a, n, o, gSf, true); },
      [&](int k){ return cdbl(gSf*argMultZ(gZc[k], gCc[k]), 0.0); },
      [&](int k){ return std::fabs(gSf)*
                         std::max(std::fabs(argMultZ(gZc[k], gCc[k])), 1.0); });
    // in1 = +-in2: the products in the imaginary part are equal, so it must
    // be exactly +0 on every path and the angle +0 or +pi (not -pi).  Built
    // with FMA on the command line (-march=native, -xHost), the compiler may
    // fuse the default version's multiply and subtraction, so there the
    // angles only need to be within 1e-7 of 0 and near pi of either sign
    {
      ++gTests;
#if defined(__FMA__)
      // the default version runs wherever no AVX level is picked, also
      // when getOptLevel() reports an SSE level
      const std::string lvl = EVM::getOptLevel();
      const bool fused = !(lvl == "AVX512" || lvl == "AVX2" ||
                           lvl == "AVX");
#else
      const bool fused = false;
#endif
      std::vector<C> neg(MAXLEN);
      std::vector<float> o0(MAXLEN), o1(MAXLEN);
      for(int i=0; i<MAXLEN; ++i) neg[i] = -gCc[i];
      EVM::multcatan2(Cc, Cc, MAXLEN, o0.data());
      EVM::multcatan2(Cc, neg.data(), MAXLEN, o1.data());
      int nbad = 0;
      for(int i=0; i<MAXLEN; ++i)
      {
        const bool ok = fused ?
            std::fabs(o0[i]) <= 1.0e-7f && std::fabs(o1[i]) > 3.0f :
            floatBits(o0[i]) == 0 && o1[i] > 3.0f;
        if(ok || ++nbad > MAX_REPORT) continue;
        printf("      element %-5d in1=(%g,%g) in2=in1: %.9g, in2=-in1: %.9g\n",
               i, gCc[i].real(), gCc[i].imag(), o0[i], o1[i]);
      }
      if(nbad) ++gFailed;
      printf("  %-34s %s  [%d values each, %s", "multcatan2 of in1 = +-in2",
             nbad ? "FAIL" : "PASS", MAXLEN,
             fused ? "fused default, within 1e-7" : "+0 and +pi exact");
      if(nbad) printf(", %d wrong", nbad);
      printf("]\n");
    }
    // special values, with a negative scale for the scaled ones, which flips
    // the sign of a zero angle
    checkAtan2Special<float>("atan2 special values", false, 1.0, false,
      [](const C* a, const C*, int n, float* o) { EVM::atan2(a, n, o); });
    checkAtan2Special<C>("atan2 complex out special values", false, 1.0,
      false,
      [](const C* a, const C*, int n, C* o) { EVM::atan2(a, n, o); });
    checkAtan2Special<float>("atan2s special values", false, -1.25,
      false,
      [](const C* a, const C*, int n, float* o)
        { EVM::atan2s(a, n, o, -1.25f); });
    checkAtan2Special<C>("atan2s complex out special values", false, -1.25,
      false,
      [](const C* a, const C*, int n, C* o) { EVM::atan2s(a, n, o, -1.25f); });
    checkAtan2Special<float>("multcatan2 special values", true, 1.0, false,
      [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2(a, b, n, o); });
    checkAtan2Special<float>("multcatan2s special values", true, -1.25,
      false,
      [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2s(a, b, n, o, -1.25f); });
    checkAtan2Special<float>("multcatan2 zeroAngle special", true, 1.0, true,
      [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2(a, b, n, o, true); });
    checkAtan2Special<float>("multcatan2s zeroAngle special", true, -1.25,
      true, [](const C* a, const C* b, int n, float* o)
        { EVM::multcatan2s(a, b, n, o, -1.25f, true); });
  }

  //--- fmmod / fmmods : scale*e^(j*2*pi*(phscale*in + phstart + i*phincr)) -
  // phases in cycles, computed in double; the SIMD versions use evm_trig.h's
  // sincosCycDbl.  The inputs are 3*gAf (up to 6 cycles) and, for the large
  // input run, the phases of gLr (up to 1e7 cycles)
  printf("\nFM modulation  (sin/cos of cycles: evm_trig.h)\n");
  {
    typedef complex<float> C;
    std::vector<float> xs(MAXLEN), xl(MAXLEN);
    for(int i=0; i<MAXLEN; ++i)
    {
      xs[i] = 3.0f*gAf[i];
      xl[i] = gLr[i].imag();
    }
    checkFmmod("fmmod", xs, 0.37, 12345.678, 0.1234567, 1.0,
      [](const float* a, int n, double ps, double p0, double p1, C* o)
        { return EVM::fmmod(a, n, ps, p0, p1, o); });
    checkFmmod("fmmod negative incr", xs, -2.5, -0.3, -0.0371, 1.0,
      [](const float* a, int n, double ps, double p0, double p1, C* o)
        { return EVM::fmmod(a, n, ps, p0, p1, o); });
    checkFmmod("fmmod large input", xl, 1.0, 0.0, 0.25, 1.0,
      [](const float* a, int n, double ps, double p0, double p1, C* o)
        { return EVM::fmmod(a, n, ps, p0, p1, o); });
    checkFmmod("fmmods", xs, 0.37, 12345.678, 0.1234567, -1.7,
      [](const float* a, int n, double ps, double p0, double p1, C* o)
        { return EVM::fmmods(a, n, ps, p0, p1, o, -1.7f); });
    checkFmmod("fmmods large phstart", xs, 1.5, 1.0e6 + 0.3, 3.7, 0.5,
      [](const float* a, int n, double ps, double p0, double p1, C* o)
        { return EVM::fmmods(a, n, ps, p0, p1, o, 0.5f); });
  }

  //--- getSize sanity -------------------------------------------------------
  printf("\nhelpers\n");
  {
    ++gTests;
    bool ok = true;
    const size_t fPer = EVM::getSize<float>(1);
    const size_t cPer = EVM::getSize<complex<float> >(1);
    // getSize must round up, never down, and never shrink a request
    for(size_t want=0; want<=257; ++want)
    {
      if(EVM::getSize<float>(want) < want)          ok = false;
      if(EVM::getSize<complex<float> >(want) < want) ok = false;
    }
    if(fPer == 0 || cPer == 0) ok = false;
    if(!ok) ++gFailed;
    printf("  %-34s %s  [rounds up, float unit %zu, complex unit %zu]\n",
           "getSize<T>", ok ? "PASS" : "FAIL", fPer, cPer);
  }

  //--- summary --------------------------------------------------------------
  printf("\n======================================================================\n");
  printf("%d tests, %d passed, %d FAILED\n", gTests, gTests-gFailed, gFailed);
  if(gFailed)
    printf("reproduce this run with:  ./evm_unit_test %u\n", seed);
  printf("======================================================================\n");
  return gFailed ? 1 : 0;
}
