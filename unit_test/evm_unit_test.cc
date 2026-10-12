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
#include <complex>
#include <vector>
#include <random>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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

// Most failing elements to print per function before summarising the rest.
static const int MAX_REPORT = 8;

//==============================================================================
// Random input data
//==============================================================================
static std::mt19937                 gRng;
static std::vector<float>           gAf, gBf;   // real input vectors
static std::vector<complex<float> > gAc, gBc;   // complex input vectors
static float                        gSf;        // real scale factor
static complex<float>               gSc;        // complex scale factor

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

static void generateData(unsigned int seed)
{
  gRng.seed(seed);
  gAf.resize(MAXLEN);  gBf.resize(MAXLEN);
  gAc.resize(MAXLEN);  gBc.resize(MAXLEN);
  for(int i=0; i<MAXLEN; ++i)
  {
    gAf[i] = randReal();
    gBf[i] = randDivisorReal();
    gAc[i] = randComplex();
    gBc[i] = randDivisorComplex();
  }
  gSf = randDivisorReal();
  gSc = randDivisorComplex();
}

//==============================================================================
// Small helpers
//==============================================================================
template<class T> struct IsComplex            { static const bool value = false; };
template<> struct IsComplex<complex<float> >  { static const bool value = true;  };

static cdbl toCdbl(const float v)          { return cdbl(v, 0.0); }
static cdbl toCdbl(const complex<float> v) { return cdbl(v.real(), v.imag()); }

template<class T> static T sentinel();
template<> float sentinel<float>()
{
  return std::numeric_limits<float>::quiet_NaN();
}
template<> complex<float> sentinel<complex<float> >()
{
  const float n = std::numeric_limits<float>::quiet_NaN();
  return complex<float>(n, n);
}

static bool isSentinel(const float v) { return std::isnan(v); }
static bool isSentinel(const complex<float> v)
{
  return std::isnan(v.real()) && std::isnan(v.imag());
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
//==============================================================================
template<class TA, class TB, class TO, class FEVM, class FREF>
static void check(const char* name, const double tol,
                  const TA* const aBuf, const TB* const bBuf,
                  FEVM runEvm, FREF baseline)
{
  ++gTests;
  const bool cplx = IsComplex<TO>::value;

  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0;
  double worstRel = 0.0;

  std::vector<TO> out(MAXLEN + GUARD);

  // mode 0: separate buffers.  mode 1: out aliases in1.  mode 2: out aliases
  // in2.  The header states the output buffer may be the same as an input
  // buffer, so the aliasing modes are part of the contract, not a bonus.
  for(int mode=0; mode<3; ++mode)
  {
    if(mode == 1 && !std::is_same<TA, TO>::value) continue;
    if(mode == 2 && (!std::is_same<TB, TO>::value || bBuf == NULL)) continue;
    const char* modeName = (mode==0) ? "" :
                           (mode==1) ? " [in-place on in1]" : " [in-place on in2]";

    for(size_t li=0; li<gLengths.size(); ++li)
    {
      const int len   = gLengths[li];
      const int valid = (len > 0) ? len : 0;

      for(int k=0; k<MAXLEN+GUARD; ++k) out[k] = sentinel<TO>();

      if(mode == 0)
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
        const cdbl   got   = toCdbl(out[k]);
        const cdbl   want  = baseline(k);
        const double adiff = std::abs(got - want);
        const double mag   = std::abs(want);
        const double rdiff = (mag > 0.0) ? adiff/mag : adiff;
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
          if(isSentinel(out[k]))
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
      for(int k=valid; k<MAXLEN+GUARD; ++k)
      {
        if(isSentinel(out[k])) continue;
        ++nstray;
        if(nshown < MAX_REPORT)
        {
          ++nshown;
          printf("      len=%-5d element %-5d WROTE PAST END evm=%s%s\n",
                 len, k, format(toCdbl(out[k]), cplx).c_str(), modeName);
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
//
// The baseline is the double precision sum of a[k]*b[k].  Summed in float in
// any order, each component of the result is within (m+1)*2^-24 times the sum
// of the magnitudes of its m product terms (the standard forward error bound;
// m = len for real inputs and 2*len for complex, whose real and imaginary
// parts each add two products per element).  The error is measured against
// that magnitude sum rather than the result, which cancellation can bring
// near zero.  The inputs are copied into NaN-filled buffers, in2 one element
// off in1's alignment, so a read past the end shows up as a NaN result.
// Negative and zero lengths must give exactly 0.
//==============================================================================
template<class T, class FEVM>
static void checkDot(const char* name, const T* const aBuf,
                     const T* const bBuf, FEVM runEvm)
{
  ++gTests;
  const bool   cplx  = IsComplex<T>::value;
  const int    terms = cplx ? 2 : 1;   // products per element per component
  const double u     = std::ldexp(1.0, -24);
  std::vector<T> a(MAXLEN + GUARD), b(MAXLEN + GUARD);

  long   compared = 0;
  int    nbad = 0, nshown = 0, nstray = 0;
  double worstRel = 0.0;

  for(size_t li=0; li<gLengths.size(); ++li)
  {
    const int len   = gLengths[li];
    const int valid = (len > 0) ? len : 0;

    for(int k=0; k<MAXLEN+GUARD; ++k) a[k] = b[k] = sentinel<T>();
    for(int k=0; k<valid; ++k) { a[k] = aBuf[k]; b[k+1] = bBuf[k]; }
    T res[2] = {sentinel<T>(), sentinel<T>()};
    runEvm(a.data(), b.data()+1, len, res);
    const cdbl got = toCdbl(res[0]);

    // baseline, and the magnitude sums that scale each component's error
    cdbl   want  = 0.0;
    double absRe = 0.0, absIm = 0.0;
    for(int k=0; k<valid; ++k)
    {
      const cdbl x = toCdbl(aBuf[k]), y = toCdbl(bBuf[k]);
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
         name, pass ? "PASS" : "FAIL", compared, cplx ? "2*" : "", worstRel);
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
  const float*          nullF = NULL;

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
