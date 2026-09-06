/* 
 * Copyright (c) 2022 Nick Xenias
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
//  Name:     polyphase.h
//
//  Purpose:  Generates a polyphase filter bank from a FIR input filter
//
//  Created:  2022/11/05
// 
//  Description:
//    Generates a polyphase filter bank from a FIR input filter.  For various
//    classes of windowed sinc filters, this class can generate the canonical
//    filter then divide it into a filter bank.  All filters in the polyphase
//    filter bank will be aligned to start on an address boundary equal to the
//    SSE or AVX register size for the given architecture so as to result in
//    efficient loads/stores.  Each filter in the bank will also be 0-padded
//    so the length of the filter is a multiple of the SSE or AVX register
//    to result in efficient dot-products with the filter.  By default, the
//    filters in the bank will be time-reversed so convolution can be achieved
//    with a time-domain dot-product.  Filters are stored internally in
//    contiguous memory (STL vector) but can be accessed either by index or
//    desired interpolation phase.
//
//    When generating windowed sinc filters, the following windowing functions
//    may be applied to the time domain taps:
//
//    Name   Sidelobe  Equivalent    Rolloff    Window Type 
//           Max (dB)  Bandwidth   (dB/Octave) 
//    NONE    -13        1.00           6       Rectangle
//    HANN    -32        1.50          18       Hann (cosine bell) 
//    HAMM    -43        1.36           6       Hamming (bell on pedestal)
//    BH61    -61        1.61           6       Blackman-Harris 3 weight  
//    BH67    -67        1.71           6       Blackman-Harris 3 weight optimal
//    BH74    -74        1.79           6       Blackman-Harris 4 weight  
//    BH92    -92        2.00           6       Blackman-Harris 4 weight optimal
//
//    Additional methods have been added to create square root raised cosine
//    filters (SRRC) which are often used for pulse shaping digital data.
//    For the SRRC methods, the rolloff factor, beta, is specified in place
//    of the window type.
//
//==============================================================================

#ifndef POLYPHASE_H
#define POLYPHASE_H

#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <complex>  // for Complex8 data type
#include "aligned_allocator.h"
#include "config.h"

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::complex;

#ifndef two_pi
#define two_pi 2*M_PI
#endif

#ifndef Complex8
#define Complex8 complex<float>
#endif

template<typename T>
class Polyphase
{
  public:
    //==============//
    // Constructors //
    //==============//
    // Default constructor creates a small windowed sinc filter bank, the total
    // number of taps in the canonical (high rate) windowed sinc filter will be
    // tapsPerFilt*numFilters but may contain some 0-valued taps on the edges
    // to ensure the filter is centered on time T=0 and that one of the taps is
    // at time T=0
    //    numFilters = desired number of filters in the polyphase bank, 
    //    tapsPerFilt = number of teps per filter in each low rate filter in
    //        the bank
    //    bw = bandwidth of the sinc interpolation filter as a fraction of
    //        the input sampling rate, i.e. at the rate of the polyphase filter
    //        not the upsampled canonical filter (a bw of .95 means 95% of the
    //        input bw, a bw > 1 is allowed but may result in aliasing)
    //    win = window function for the sinc interpolation filter, see
    //        genfilter.h for choices as this class is actually used to
    //        generate the canonical filter
    Polyphase(const int numFilters=3, const int tapsPerFilt=3,
        const double bw=1.0, const GenFilter::WinType win=GenFilter::HANN);

    // Constructor for generating a SRRC filter bank, the total number of taps
    // in the canonical (high rate) SRRC filter will be tapsPerFilt*numFilters
    // but may contain some 0-valued taps on the edges to ensure the filter is
    // centered on time T=0 and that one of the taps is at T=0
    //    numFilters = desired number of filters in the polyphase bank, 
    //    tapsPerFilt = number of teps per filter in each low rate filter in
    //        the bank
    //    bw = bandwidth of the sinc interpolation filter as a fraction of
    //        the input sampling rate, i.e. at the rate of the polyphase filter
    //        not the upsampled canonical filter (a bw of .95 means 95% of the
    //        input bw, a bw > 1 is allowed but may result in aliasing)
    //    rolloff = rollof factor for SRRC pulse shape, must be between 0 and 1
    //        inclusive, a rolloff of 0 is a perfect sinc, a rolloff of 1 has
    //        100% excess bw from rolloff 0
    Polyphase(const int numFilters, const int tapsPerFilt,
        const double bw, const double rolloff);
    
    // Constructor for generating a filter bank from a user defined canonical
    // filter, this method requires the filter to span the time T=0 though it
    // does not require the filter to specifically have a tap fall on that time
    //    filter = canonical filter taps to generate bank from
    //    len = number of taps in the filter
    //    delay = filter delay in seconds, i.e. -1*time of first tap
    //    sampInterval = sampling interval in seconds between filter taps
    //    numFilters = number of filters to split the canonical filter into
    //    initPhase = desired interpolation phase of the first filter in the
    //        bank, filters will be ordered in increasing phase increments of
    //	      1.0/numFilters cycles
    Polyphase(T const * const filter, const int len, const double delay,
        const double sampInterval, const int numFilters,
        const double initPhase);

  private:
    int _tapsPerFilter;   // number of taps in each low-rate filter
    int _sampsPerFilter;  // size in samples of each low-rate filter with any
                          // additional padding for SIMD register size
    int _numFilters;      // number of filters in the bank, interpolation phase
                          // resolution will be 1/_numFilters
    double _initPhase;    // interpolation phase in cycles of the first filter
                          // in the bank, in the range (-1 : 1)
    bool _extraFilt;      // true if an extra filter with interpolation phase
                          // _initPhase+1 cycles was added to the end of the
                          // bank
    ALIGNED_VECTOR(T) _vfilts;  // bank of low-rate interpolation filters,
                                // stored back to back with _sampsPerFilt
                                // samples of size(T) between the start of one
                                // filter and the start of the next
    
    // Splits a high-rate canonical filter into a polyphase bank of low-rate
    // filters that can be used to interpolate all fractional sample offsets
    // in the range [initPhase : 1+initPhase).  This method requires initPhase
    // to be in the range (-1 : 1) and requires the filter to span time T=0
    // (i.e. ts <= 0 and ts+len*sampInterval >= 0
    //    filter = canonical filter taps to generate bank from
    //    len = number of taps in the filter
    //    ts = time in seconds of the first tap of the canonical filter,
    //        typically will be a negative number, and if the filter is
    //        centered on time T=0, will be -delay of the filter bank
    //    sampInterval = sampling interval in seconds between filter taps
    //    numFilters = number of filters to split the canonical filter into
    //    initPhase = desired interpolation phase of the first filter in the
    //        bank, must be in the range (-1 : 1), filters will be ordered in
    //        increasing phase increments of 1.0/numFilters cycles
    //    evmPad = if true, pad each filter in the bank with enough zeroes so
    //        it's length is a multiple of the SSE/AVX register size
    //    Polyphase::_genFilterBank = whole sample delay of the low-rate
    //        polyphase filters (i.e. the sample index the bank will
    //       interpolate from)
    int _generateFilterBank(T const * const filter, const int len,
        const float ts, const float sampInterval, const int numFilters,
        const double initPhase=-.5, const bool evmPad=false);
};

// Splits a high-rate canonical filter into a polyphase bank of low-rate
// filters that can be used to interpolate all fractional sample offsets
// in the range [initPhase : 1+initPhase).  This method requires initPhase
// to be in the range (-1 : 1) and requires the filter to span time T=0
// (i.e. ts <= 0 and ts+len*sampInterval >= 0
//    filter = canonical filter taps to generate bank from
//    len = number of taps in the filter
//    ts = time in seconds of the first tap of the canonical filter,
//        typically will be a negative number, and if the filter is
//        centered on time T=0, will be -delay of the filter bank
//    sampInterval = sampling interval in seconds between filter taps
//    numFilters = number of filters to split the canonical filter into
//    initPhase = desired interpolation phase of the first filter in the
//        bank, must be in the range (-1 : 1), filters will be ordered in
//        increasing phase increments of 1.0/numFilters cycles
//    evmPad = if true, pad each filter in the bank with enough zeroes so
//        it's length is a multiple of the SSE/AVX register size
//    Polyphase::_genFilterBank = whole sample delay of the low-rate
//        polyphase filters (i.e. the sample index the bank will
//        interpolate from)
template<typename T>
int Polyphase::_generateFilterBank(T const * const filter, const int len,
    const bool isCX, const float ts, const float sampInterval,
    const int numFilters, const double initPhase, const bool evmPad)
{
  // enforce restrictions on initPhase and ts
  if(initPhase <= -1.0 || initPhase >= 1.0)
  {
    throw std::range_error("[Polyphase::_generateFilterBank] Invalid value "
        "specified for initial filter phase, phase must be in the range "
        "(-1 : 1) cycles, non-inclusive");
  }
  if(ts > 0 || ts+len*sampInterval < 0)
  {
    throw std::range_error("[Polyphase::_generateFilterBank] Invalid start "
      "specified for canonical filter, filter must span the time T=0, so "
      "start time ts must be in the range "
      "[-len*sampInterval : len*sampInterval]");
  }
  // the basic idea is to split the high rate canonical filter into
  // <numFilters> lower rate interpolation filters (i.e. downsample by
  // numFilters).  Each lower rate filter will have the same whole sample delay
  // as the others but with increasing interpolation phase, meaning that given
  // a single buffer of input data of length _tapsPerFilter, they can be used
  // to interpolate samples delay+initPhase to delay+initPhase+1 in increments
  // of 1/numFilters cycles.
  //
  // In order to achieve the desired initPhase of the bank and to ensure the
  // legnth of the canonical filter is an integer multiple of numFilters for
  // downsampling, we may need to zero-pad both the beginning and end of the
  // canonical filter.  Note that adding zeroes to the end of the canonical
  // filter will increase the interpolation phase of the first filter in the
  // bank by 1/numFilters cycles for each zero (since we time-reverse the
  // filter for convolution prior to downsampling), but the only way to
  // decrease the interpolation phase of the first filter is to increase the
  // delay of all the low rate filters by 1 sample, which will decrease the
  // initial interpolation phase by 1 cycle. 
  _numFilters = numFilters; // split the bank into this many filters
  // compute the lower-rate sampling interval of each filter in the bank
  const double lrInterval = sampInterval*_numFilters;
  // compute the taps per filter and the minimum amount of zero-padding that
  // must be done so the length of the canonical filter is divisible by the
  // number of filters in the bank (i.e. the downsampling factor)
  _tapsPerFilter = len/_numFilters;
  if(len%_numFilters) ++_tapsPerFilter;
  int zpad = _tapsPerFilter*_numFilters-len;
  // compute the start time (phi) and interpolation phase (theta) in cycles of
  // the first filter in the bank if we were to downsample the canonical filter
  // as-is, this will be:
  //    phi = (ts+(L-1)*Thr)/Tlr
  //    theta = phi-floor(phi)
  // where ts is the time in seconds of the first sample in the canonical
  // filter (i.e. -delay), Thr is the sampling interval of the canonical (i.e.
  // high-rate) filter and Tlr is the sampling interval of the downsampled
  // (i.e. low-rate) filters equal to numFilters*sampInterval
  double phi = (ts+(len-1)*sampInterval)/lrInterval;
  // since we've guaranteed the filter spans time T=0, ts is guaranteed to be
  // negative and phi is guaranteed to be positive, which means we can use
  // trunc() as a substitute for floor()
  int nd = trunc(phi); // whole-sample delay of the low-rate filters
  double theta = phi-nd;
  // since our desired value of theta is initPhase and the filter will be
  // time-reversed for convolution, we need to adjust the low-rate delay and/or
  // add zeroes to the end of the canonical filter to adjust the interpolation
  // phase.  To get theta==initPhase we need:
  //    nzeroes = round((initPhase-theta)*numFilters)
  // zeroes added to the end of the canonical filter.  Since we can't take
  // samples away from the end of the filter, only add zeroes, we need to
  // ensure nzeroes is > 0, i.e. that theta < initPhase.  To achieve this,
  // we can reduce theta by 1 full cycle at a time by increasing the delay of
  // the low rate filter by 1 sample (i.e. shifting the reference sample we're
  // interpolating from forward by 1, for example interpolating sample index
  // 2+.7 is the same as interpolation sample index 3-.3)
  while(theta > initPhase)
  {
    ++nd;
    theta -= 1.0;
  }
  // now compute the number of zeroes required at the end of the filter to
  // increase theta to initPhase
  const int nzreq = int(round((initPhase-theta)*_numFilters));
  // if the number of required zeroes is greater than the minimum number of
  // zeroes necessary for the filter length to be an integer multiple of
  // _numFilters, increase the minimum number by _numFilters, we'll add nzreq
  // zeroes to the end of the filter and zpad-nzreq to the beginning to get
  // an integer multiple of _numFilters taps while also setting the initial
  // interpolation phase to initPhase
  if(nzreq > zpad) zpad += _numFilters;
  _tapsPerFilter = (len+zpad)/_numFilters;
  newTs = ts-(zpad-nzreq)*sampInterval; // compute new ts with start pad
  // recompute theta using the new canonical filter start time with padding
  phi = (newTs+(len+zpad-1)*sampInterval)/lrInterval;
  nd = trunc(phi);
  theta = phi-nd;
  // theta should now be either equal to initPhase (if initPhase was positive)
  // or initPhase+1 (if initPhase was negative), if theta is initPhase+1,
  // reduce it by 1 cycle by increasing the filter delay (shifting the low
  // rate filter back by a sample)
  if(initPhase < 0)
  {
    ++nd;
    theta -= 1.0;
  }
  assert(initPhase-theta < 1/_numFilters);
  _initPhase = theta;
  // copy the canonical filter into a buffer with the desired front and back
  // zero padding
  vector<T> vcanonical(len+zpad);
  memset(&vcanonical[0], 0, sizeof(T)*vcanonical.size();)
  memcpy(&vcanonical[zpad-nzreq], filter, len*sizeof(T));
  // compute how many samples to pad each filter by to ensure length is a
  // multiple of the AMX/AVX/SSE register size
  _vfilts.resize();
  // now generate the filter bank in time-reversed order by starting with the
  // last sample, note that we've guaranteed that len+zpad is an integer
  // multiple of _tapsPerFilter
  for(int i=0; i<_numFilters; ++i) // for each filter in the bank
  {
    for(int j=0; j<_tapsPerFilter)
  }

  return nd;
}

#endif // POLYPHASE_H
