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

// Relative error tolerances by algorithm class.  The measured worst case for
// each class on AVX2/AVX/scalar is given for reference.
static const double TOL_EXACT = 1.0e-6;   // true multiply/divide   (~1.7e-07)
static const double TOL_NR    = 5.0e-6;   // reciprocal + 1x Newton (~2.5e-07)
static const double TOL_RCP   = 2.0e-3;   // bare reciprocal        (~2.9e-04)
static const double TOL_CMP   = 0.0;      // comparison/selection   (exact)
static const double TOL_TRIG  = 2.0e-6;   // sin/cos/atan2          (~8.1e-07)

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
// small input lands there).
//
// magnitude(k), if given, takes the place of |baseline|.  multSum passes
// |in1*in2| + |in3|: its result can cancel to near zero, where the rounding
// error of the product is large next to the result but not next to the terms
// (the dot-products are measured the same way).
//==============================================================================
template<class TA, class TB, class TO, class FEVM, class FREF,
         class FMAG = std::nullptr_t>
static void check(const char* name, const double tol,
                  const TA* const aBuf, const TB* const bBuf,
                  FEVM runEvm, FREF baseline, FMAG magnitude = nullptr)
{
  ++gTests;
  const bool cplx = IsComplex<TO>::value;

  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0;
  double worstRel = 0.0;

  std::vector<TO> out(MAXLEN + GUARD);
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
    if(mode == 3 && !cplxToReal) continue;
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
//
// The baseline is the double precision sum of a[k]*b[k].  Summed in float in
// any order, each component of the result is within (m+1)*2^-24 times the sum
// of the magnitudes of its m product terms (the standard forward error bound;
// m = 2*len when both inputs are complex, whose real and imaginary parts each
// add two products per element, otherwise m = len).  The error is measured
// against that magnitude sum rather than the result, which cancellation can
// bring near zero.  The inputs are copied into NaN-filled buffers, in2 one
// element off in1's alignment, so a read past the end shows up as a NaN
// result.  Negative and zero lengths must give exactly 0.
//==============================================================================
template<class TA, class TB, class FEVM>
static void checkDot(const char* name, const TA* const aBuf,
                     const TB* const bBuf, FEVM runEvm,
                     const bool conjB = false)
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
      const cdbl y = conjB ? std::conj(toCdbl(bBuf[k])) : toCdbl(bBuf[k]);
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
  // against the magnitude, since one part can be near zero.  The SIMD
  // versions need SVML, so they are only tested with the Intel compiler
  printf("\npolar / cartesian conversion  (SIMD versions need SVML: %s)\n",
#if defined(EVM_SVML)
         "yes");
#else
         "no, default version only");
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
