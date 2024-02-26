/* 
 * Copyright (c) 2013 Nick Xenias
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

// implement an equalizer for complex data

#ifndef EQUALIZER_H
#define EQUALIZER_H

#include <vector>
#include "constel.h"
#include <cmath>
#include <mkl.h>

class Equalizer
{
  public:
    Equalizer(int eqLen=21, float mu=.001);
    Equalizer(const complex_8* initTaps, int eqLen, float mu=.001);
    inline void reset(int eqLen, float mu);
    inline void reset(const complex_8* initTaps, int eqLen, float mu);
    inline int train(complex_8* vin, complex_8* vtrain, int len);
    inline void equalize(complex_8* vin, int inLen, complex_8* vtrain,
        int trLen, complex_8* vout);  
    inline void equalize(complex_8* vin, int inLen, complex_8* vout);
    inline void eqNoAdapt(complex_8* vin, int inLen, complex_8* vout);
    inline void setMu(float mu){_mu=mu;};
    inline complex_8 quantize8psk(complex_8 inPt)
        {return Constel::PSK8[_quantize8psk(inPt)];};
    inline void getTaps(complex_8* taps)
        {memcpy(taps,&_eqTaps[0],_eqLen*_dSize);};

  private:
    int _eqLen;       // equalizer length
    int _eqLenM1;     // _eqLen-1 to save computation time
    float _mu;        // adaptive algorithm stepsize
    const int _dSize; // size of a data element in bytes (sizeof(complex_8))
    complex_8 _err;   // error between soft and hard decision
    vector<complex_8> _eqTaps;    // equalizer taps
    vector<complex_8> _mtrxAHA;   // covariance matrix for LS algorithm
    vector<complex_8> _vctrAHsd;  // vector AH x training
    
    inline int _quantize8psk(complex_8 inPt);
};

Equalizer::Equalizer(int eqLen, float mu) : _dSize(sizeof(complex_8))
{
  reset(eqLen, mu);
}

Equalizer::Equalizer(const complex_8* initTaps, int eqLen, float mu) : _dSize(sizeof(complex_8))
{
  reset(initTaps, eqLen, mu);
}

// reset EQ taps to center spike initialization
inline void Equalizer::reset(int eqLen, float mu)
{
  _mu = mu;
  _eqLen = eqLen;
  _eqLenM1 = _eqLen-1;
  _eqTaps.resize(_eqLen);
  fill(_eqTaps.begin(),_eqTaps.end(),complex_8(0,0));
  // center spike initialization
  _eqTaps[_eqLen/2] = complex_8(1,0);
}

// reset EQ taps to user-specified initialization vector
inline void Equalizer::reset(const complex_8* initTaps, int eqLen, float mu)
{
  _mu = mu;
  _eqLen = eqLen;
  _eqLenM1 = _eqLen-1;
  _eqTaps.resize(_eqLen);
  memcpy(&_eqTaps[0],initTaps,_eqLen*_dSize);
}

// compute optimal equalizer using LS method.  Return delay value that
// resulted in optimal eq taps
inline int Equalizer::train(complex_8* vin, complex_8* vtrain, int len)
{
  int idx = 0;                  // index for packed matrix AHA creation
  int multElems = len-_eqLen+1; // number of elements in the dotproducts for AHA
  int info;                     // MKL method error flag
  int delay = -1;               // optimal delay
  int incr = 1;                 // for MKL methods access by reference
  char uplo = 'L';              // all MKL triangular matrices will be lower
  char trans = 'N';             // solve system of equations Af=y for f
  char diag = 'N';              // whether matrix can be assumed unit triangular
  double mse, mmse;             // mean squared and minimum mean squared error
  complex_8 cxPt;               // interim result for MSE calculation
  float scale = (1.0/multElems);
  double tapScale = 0;
  for(int i=0; i<len; ++i) tapScale += abs2(vin[i]);
  tapScale =  len/tapScale;
  // compute covariance matrix AHA where:
  //  A (N-m+1 x m) =  [ x(m-1)  x(m-2)  . . . x(0)   ]
  //                   [ x(m)    x(m-1)  . . . x(1)   ]
  //                   [   .                          ]
  //                   [   .                          ]
  //                   [   .                          ]
  //                   [ x(N-1)  x(N-2)  . . . x(N-m) ]
  //
  //  AH (m x N-m+1) = [ x*(m-1)  x*(m)   x*(m+1)  . . .  x*(N-1) ]
  //                   [ x*(m-2)  x*(m-1) x*(m)    . . .  x*(N-2) ]
  //                   [   .                                      ]
  //                   [   .                                      ]
  //                   [   .                                      ]
  //                   [ x*(0)    x*(1)   x*(2)    . . .  x*(N-m) ]
  //
  // and x is the elements of vin, m is _eqLen, and N is the training length
  // (len).
  // setup matrix AHA in MKL lower packed format (since it is hermitian)
  _mtrxAHA.resize((_eqLen*(_eqLen+1))/2);
  _vctrAHsd.resize(_eqLen);
  for(int i=0; i<_eqLen; ++i)
  {
    for(int j=i; j<_eqLen; ++j)
    {
      // (size, x, xincr, y, yincr, result)
      cblas_cdotc_sub(multElems,&vin[_eqLenM1-j],1,&vin[_eqLenM1-i],1,
          &_mtrxAHA[idx++]);
      _mtrxAHA[idx-1] *= tapScale*scale;
    }
  }
  // perform Cholesky decomposition on AHA
  // (uplo, order of AHA, AHA, info)
  cpptrf(&uplo,&_eqLen,reinterpret_cast<MKL_Complex8*>(&_mtrxAHA[0]),&info);
  if(info < 0)
  {
    cerr<<"[Equalizer] Error: illegal value encountered during Cholesky"
        <<" facotrization, matrix can not be factorized."<<endl;
  }
  else if(info > 0)
  {
    cerr<<"[Equalizer] Error: matrix is not positive definite."<<endl;
  }
  // for each delay value, d, compute AHs(d) and solve system of linear
  // equations to find optimal equalizer taps
  mmse = -1;  // set to impossible value
  for(int d=0; d<_eqLen; ++d)
  {
    // compute AHsd
    for(int i=0; i<_eqLen; ++i)
    {
      // (size, x, xincr, y, yincr, result)
      cblas_cdotc_sub(multElems,&vin[_eqLenM1-i],1,&vtrain[_eqLenM1-d],1,
          &_vctrAHsd[i]);
      _vctrAHsd[i] *= scale;
    }
    // solve for equalizer taps
    // (uplo, system description, unit triang, order AHA, AHA, sd, incr sd)
    ctpsv(&uplo,&trans,&diag,&_eqLen,
        reinterpret_cast<MKL_Complex8*>(&_mtrxAHA[0]),
        reinterpret_cast<MKL_Complex8*>(&_vctrAHsd[0]),&incr);
    // scale taps by variance of vin
    //for(int i=0; i<_eqLen; ++i) _vctrAHsd[i]*=tapScale;
    // compute MSE resulting from these taps
    mse = 0;
    for(int i=_eqLenM1; i<len; ++i)
    {
      cxPt.re = 0;
      cxPt.im = 0;
      for(int j=0; j<_eqLen; ++j) cxPt += _vctrAHsd[j]*vin[i-j];
      mse += abs2(vtrain[i]-cxPt);
    }
    // save taps if this is the lowest MSE
    if(mse < mmse || mmse < 0)
    {
      mmse = mse;
      delay = d;
      memcpy(&_eqTaps[0],&_vctrAHsd[0],_eqLen*_dSize);
    }
  }
  // return the delay that result in the MMSE taps
  return delay;
}

// equalize with training sequence
inline void Equalizer::equalize(complex_8* vin, int inLen,
        complex_8* vtrain, int trLen, complex_8* vout)
{
  memset(vout,0,inLen*_dSize);
  // can not compute first output until we have _eqLen-1 inputs
  // adapt forwards over training portion
  for(int i=_eqLenM1; i<trLen; ++i)
  {
    // compute new output point and error term
    for(int j=0; j<_eqLen; ++j) vout[i] += _eqTaps[j]*vin[i-j];
    _err = vtrain[i]-vout[i];
    // update each tap of the eq
    for(int j=0; j<_eqLen; ++j) _eqTaps[j] += (_mu*_err*~vin[i-j]);
  }

  // adapt forwards over data portion
  for(int i=trLen; i<inLen; ++i)
  {
    // compute new output point and error term
    for(int j=0; j<_eqLen; ++j) vout[i] += _eqTaps[j]*vin[i-j];
    _err = Constel::PSK8[_quantize8psk(vout[i])]-vout[i];
    // update each tap of the eq
    for(int j=0; j<_eqLen; ++j) _eqTaps[j] += (_mu*_err*~vin[i-j]);
  }
}

// equalize without training sequence
inline void Equalizer::equalize(complex_8* vin, int inLen, complex_8* vout)
{
  memset(vout,0,inLen*_dSize);
  // adapt forwards over data portion
  for(int i=_eqLenM1; i<inLen; ++i)
  {
    // compute new output point and error term
    for(int j=0; j<_eqLen; ++j) vout[i] += _eqTaps[j]*vin[i-j];
    _err = Constel::PSK8[_quantize8psk(vout[i])]-vout[i];
    // update each tap of the eq
    for(int j=0; j<_eqLen; ++j) _eqTaps[j] += (_mu*_err*~vin[i-j]);
  }
}

// equalize using the current EQ taps without adaptive updates
inline void Equalizer::eqNoAdapt(complex_8* vin, int inLen, complex_8* vout)
{
  memset(vout,0,inLen*_dSize);
  // can not compute first output until we have _eqLen-1 inputs
  // equalize
  for(int i=_eqLenM1; i<inLen; ++i)
  {
    for(int j=0; j<_eqLen; ++j) vout[i] += _eqTaps[j]*vin[i-j];
  }
}

inline int Equalizer::_quantize8psk(complex_8 inPt)
{
  return int( ((int(floor(inPt.re))>>Constel::PSK8_SHFT1)&Constel::PSK8_MASK1)
             |((int(floor(inPt.im))>>Constel::PSK8_SHFT2)&Constel::PSK8_MASK2)
             |((int(floor(fabsf(inPt.re)-fabsf(inPt.im)))>>Constel::PSK8_SHFT3)&
              Constel::PSK8_MASK3));
}

#endif
