/* 
 * Copyright (c) 2010 Nick Xenias
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
//  Name:     genfilter.h
//
//  Purpose:  Collection of methods to generate real or complex filters using
//            the windowed sinc method
//
//  Created:  2010/10/18
// 
//  Description:
//    Creates various GLP FIR filters, returns time domain taps as STL vector.
//    All filters are created as a prototype lowpass (sinc) impulse response
//    then modulated to the desired frequency.  A tapering window is also
//    applied to the impulse response to achieve different levels of stopband
//    rejection and excess BW.  All filters will be symmetric about t=0.
//
//    The following windowing functions may be applied to the time domain taps
//    bia the provided WinType enumeration:
//
//    Window Type                           Sidelobe  Equivalent    Rolloff
//                                          Max (dB)  Bandwidth   (dB/Octave) 
//    Rectangle (NONE)                        -13        1.00           6
//    Hann cosine bell (HANN)                 -32        1.50          18
//    Hamming bell on pedestal (HAMM)         -43        1.36           6
//    Blackman-Harris 3 weight (BH61)         -61        1.61           6
//    Blackman-Harris 3 weight optimal (BH67) -67        1.71           6
//    Blackman-Harris 4 weight (BH74)         -74        1.79           6
//    Blackman-Harris 4 weight optimal (BH92) -92        2.00           6
//
//    Additional methods have been added to create square root raised cosine
//    filters (SRRC) which are often used for pulse shaping digital data.
//    For the SRRC methods, the rolloff factor, beta, is specified in place
//    of the window type.
//
//    This namespace also provides methods to generate windowed sinc
//    interpolation filters which can be used to interpolate a specific
//    fractional sample delay (i.e. can be used to delay an input buffer by
//    a fraction of the sample interval).  These methods may be useful for
//    creating polyphase filter banks for things like resamplers.
//
//==============================================================================

#ifndef GEN_FILTER_H
#define GEN_FILTER_H

#include <vector>
#include <math.h>
#include <string>
#include <stdexcept>
#include <algorithm> // std::transform
#include <complex>

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::complex;

namespace GenFilter
{
  // mathematical constants
  const float pi_f = 3.141592654f;
  const double pi_d = 3.141592653589793238462643e+0;
  const float two_pi_f = 6.283185308f;
  const double two_pi_d = 6.283185307179586476925287e+0;
  
  // window type enumeration
  enum WinType {UNKNOWN=-1,NONE,SRRC,HAMM,HANN,BH61,BH67,BH74,BH92};
  
  // convert string to window type, case insensitive, only first 4 letters
  // of input string will be parsed, if the input string does not match
  // a valid window type, GenFilter::UNKNOWN will be returned.
  //    winStr = string representation of tapering window type (i.e. "BH92")
  //    GenFilter::str2win = WinType enumeration of tapering window
  inline WinType str2win(string winStr);
  
  // convert a GenFilter::WinType enumeration value to its string
  // representation.
  //    win = window type enumeration value
  //    GenFilter::win2str = string representation of window type
  inline string win2str(WinType win);
  
  // creates the desired tapering window
  //    vwin = vector containing window after call
  //    len = length in taps of the window, time t=0 assumed to be len/2.0
  //    win = window name to create, first 4 letters only, can be upper 
  //              or lower case
  template <typename T>
  inline void createWindow(vector<double, T> &vwin, int len, WinType win);
  template <typename T>
  inline void createWindow(vector<float,T> &vwin, int len, WinType win);
  
  // returns the tap of the desired window at a specific phase (delay).  The
  // phase should be specified between -.5 and .5 inclusive where -.5 and .5
  // would be the points where the window value approaches 0 on either end
  // (since these are tapering windows)
  //    phase = phase in cycles between -.5 and .5 to generate window tap for
  //    win = tapering window function
  //    GenFilter::windowTap = value of the window at the specified phase
  inline double windowTap(const double phase, const WinType win);
  
  // generate a lowpass filter with real-valued taps, time domain taps are
  // computed directly as a sinc then windowed
  //    taps = impulse response length in samples
  //    cutoff = LPF cutoff frequency in Hz
  //    xdelta = filter sampling interval in seconds
  //    win = window to apply to impulse response
  //    filtOut = filter impulse response
  //    GenFilter::createLPF = filter delay in seconds
  template <typename T>
  inline double createLPF(int taps, double cutoff, double xdelta,
      WinType win, vector<float,T> &filtOut)
  {
    // check validity of cutoff and xdelta
    if(cutoff > 1/(2*xdelta) || cutoff < 0)
    {
      throw std::logic_error("[GenFilter::createLPF] Invalid cutoff frequency "
          "specified");
    }
    // setup filter parameters
    filtOut.resize(taps);
    double delay = double(taps-1)/2.0;  // starting tap sample index
    double tTs;         // time value t/Ts
    vector<double> vwin(taps,0);  // window function
    createWindow(vwin,taps,win);
    // compute the windowed sinc
    for(int i=0; i<taps; ++i)
    {
      tTs = (i-delay)*xdelta*2*cutoff;
      if(tTs == 0)  // center tap (t=0)
      {
        filtOut[i] = (2*cutoff*xdelta)*vwin[i]; // BW/Fs
      }
      else // all other taps
      {
        // BW/Fs factor cancels
        filtOut[i] = sin(pi_d*tTs)/(pi_d*(i-delay))*vwin[i];
      }
    }
    return delay*xdelta;
  }
  
  // generate a Lowpass filter with complex-valued taps
  // this method just provided for consistency, shouldn't need it in most
  // practical applications as it is still a purely real filter, just stored
  // in a complex buffer with imaginary part of each element set to 0
  //    taps = impulse response length in samples
  //    cutoff = LPF cutoff frequency in Hz
  //    xdelta = filter sampling interval in seconds
  //    win = window to apply to impulse response
  //    filtOut = filter impulse response
  //    GenFilter::createLPF = filter delay in seconds
  template <typename T>
  inline double createLPF(int taps, double cutoff, double xdelta,
      WinType win, vector<complex<float>,T> &filtOut)
  {
    // generate a real valued filter, then set imaginary part of each tap to 0
    vector<float> reTaps;
    double delay = createLPF(taps,cutoff,xdelta,win,reTaps);
    filtOut.resize(reTaps.size());
    for(unsigned int i=0; i<reTaps.size(); ++i)
    {
      filtOut[i].real(reTaps[i]);
      filtOut[i].imag(0);
    }
    return delay;
  }
  
  // generate a Bandpass filter with real-valued taps,
  // this filter is generated by creating a lowpass filter then modulating
  // it up to the correct frequency with a cosine
  //    taps = impulse response length in samples
  //    cutoff1 = 1st BPF cutoff frequency in Hz
  //    cutoff2 = 2nd BPF cutoff frequency in Hz
  //    xdelta = filter sampling interval in seconds
  //    win = window to apply to impulse response
  //    filtOut = filter impulse response
  //    GenFilter::createBPF = filter delay in seconds
  template <typename T>
  inline double createBPF(int taps, double cutoff1, double cutoff2,
      double xdelta, WinType win, vector<float,T> &filtOut)
  {
    // check validity of cutoff and xdelta
    if(cutoff2 > 1/(2*xdelta) || cutoff2 <= cutoff1 || cutoff1 < 0)
    {
      throw std::logic_error("[GenFilter::createBPF] Invalid cutoff frequency "
          "specified");
    }
    // setup filter parameters, LPF cutoff will be 1/2 BPF bandwidth
    double delay = createLPF(taps,(cutoff2-cutoff1)/2,xdelta,win,filtOut);
    // modulate filter to desired CF with cos so center of passband has phase=0
    double phase = 0;
    // 2*pi*Ts*Fc (2's cancel)
    double phaseIncr = pi_d*xdelta*(cutoff2+cutoff1);
    // to maintain proper delay, need to make t=0 modulated with cos(0),
    // fmod(delay,xdelta) should work but was not giving the right answer
    // when compiled with certain optimization flags, so instead a direct
    // even/odd filter tap check is used.  Also have to set start using
    // integer math in a branch rather than ceil or floor as agressive compiler
    // optimizations like --no-prec-div may cause the ceil and floor methods
    // to return the wrong answer for certain odd values of taps (ex. 255)
    int start = taps/2;
    if(taps%2 == 0) phase = phaseIncr/2;
    else
    {
      phase = 0;
      start = (taps-1)/2;
    }
    for(int i=start; i<taps; ++i)
    {
      filtOut[i] *= 2*cos(phase); // scale by 2 for LPF -> BPF conversion
      phase += phaseIncr;
    }
    if(taps%2 == 0)
    {
      phase = -phaseIncr/2;
      --start; // need to start just prior to t=0
    }
    else phase = 0;
    for(int i=start; i>=0; --i)
    {
      filtOut[i] *= 2*cos(phase); // scale by 2 for LPF -> BPF conversion
      phase -= phaseIncr;
    }
    if(taps%2!=0) filtOut[(taps-1)/2] *= .5; // t=0 tap got scaled twice
    return delay;
  }
  
  // generate a Bandpass filter with complex-valued taps,
  // this filter is constructed by generating a lowpass filter then
  // modulating it to the correct frequency with a complex exponential (if
  // singleSided==true) or a cosine (if singleSided==false)
  //    taps = impulse response length in samples
  //    cutoff1 = 1st BPF cutoff frequency in Hz
  //    cutoff2 = 2nd BPF cutoff frequency in Hz
  //    xdelta = filter sampling interval in seconds
  //    win = window to apply to impulse response
  //    filtOut = filter impulse response
  //    singleSided = whether to create a single sided frequency response
  //    GenFilter::createBPF = filter delay in seconds
  template <typename T>
  inline double createBPF(int taps, double cutoff1, double cutoff2,
      double xdelta, WinType win, vector<complex<float>,T> &filtOut,
      bool singleSided=true)
  {
    // check validity of cutoff and xdelta
    if(cutoff2 > 1/(2*xdelta) || cutoff2 <= cutoff1 || cutoff1 < -1/(2*xdelta))
    {
      throw std::logic_error("[GenFilter::createBPF] Invalid cutoff frequency "
          "specified");
    }
    // check validity of cutoff and xdelta
    double delay = createLPF(taps,(cutoff2-cutoff1)/2,xdelta,win,filtOut);
    // modulate filter to desired CF with complex exponential, center of
    // passband will have phase=0
    double phase = 0;
    // 2*pi*Ts*Fc (2's cancel)
    double phaseIncr = pi_d*xdelta*(cutoff2+cutoff1);
    int start = taps/2;
    if(singleSided) // modulate with complex sinusoid
    {
      // to maintain proper delay, need to make t=0 modulated with cos(0),
      // fmod(delay,xdelta) should work but was not giving the right answer
      // when compiled with certain optimization flags, so instead a direct
      // even/odd filter tap check is used.  Also have to set start using
      // integer math in a branch rather than ceil/floor as agressive compiler
      // optimizations like --no-prec-div may cause the ceil and floor methods
      // to return the wrong answer for certain odd values of taps (ex. 255)
      if(taps%2 == 0) phase = phaseIncr/2;
      else
      {
        phase = 0;
        start = (taps-1)/2;
      }
      for(int i=start; i<taps; ++i)
      {
        filtOut[i] *= complex<float>(cos(phase),sin(phase));
        phase += phaseIncr;
      }
      if(taps%2 == 0)
      {
        phase = -phaseIncr/2;
        --start; // need to start just prior to t=0
      }
      else phase = 0;
      for(int i=start; i>=0; --i)
      {
        filtOut[i] *= complex<float>(cos(phase),sin(phase));
        phase -= phaseIncr;
      }
    }
    else  // double sided, modulate with cosine
    {
      // to maintain proper delay, need to make t=0 modulated with cos(0),
      // fmod(delay,xdelta) should work but was not giving the right answer
      // when compiled with certain optimization flags, so instead a direct
      // even/odd filter tap check is used.  Also have to set start using
      // integer math in a branch rather than ceil/floor as agressive compiler
      // optimizations like --no-prec-div may cause the ceil and floor methods
      // to return the wrong answer for certain odd values of taps (ex. 255)
      if(taps%2 == 0) phase = phaseIncr/2;
      else
      {
        phase = 0;
        start = (taps-1)/2;
      }
      for(int i=start; i<taps; ++i)
      {
        // scale by 2 for LPF -> BPF conversion
        filtOut[i].real(filtOut[i].real()*2*cos(phase));
        phase += phaseIncr;
      }
      if(taps%2 == 0)
      {
        phase = -phaseIncr/2;
        --start; // need to start just prior to t=0
      }
      else phase = 0;
      for(int i=start; i>=0; --i)
      {
        // scale by 2 for LPF -> BPF conversion
        filtOut[i].real(filtOut[i].real()*2*cos(phase));
        phase -= phaseIncr;
      }
      // t=0 tap got scaled twice
      if(taps%2!=0) filtOut[(taps-1)/2].real(filtOut[(taps-1)/2].real()*.5);
    }
    return delay;
  }
  
  // generate a Highpass filter with real-valued taps, this filter is first
  // created as a prototype lowpass filter then every other tap is negated
  // to shift the filter to the nyquist rate
  //    taps = impulse response length in samples
  //    cutoff = HPF cutoff frequency in Hz
  //    xdelta = filter sampling interval in seconds
  //    win = window to apply to impulse response
  //    filtOut = filter impulse response
  //    GenFilter::createHPF = filter delay in seconds
  template <typename T>
  inline double createHPF(int taps, double cutoff, double xdelta,
      WinType win, vector<float,T> &filtOut)
  {
    // check validity of cutoff and xdelta
    if(cutoff > 1/(2*xdelta) || cutoff < 0)
    {
      throw std::logic_error("[GenFilter::createHPF] Invalid cutoff frequency "
          "specified");
    }
    // create prototype lowpass, reflect cutoff freq around Fs/4
    double delay = createLPF(taps,1.0/(2*xdelta)-cutoff,xdelta,win,filtOut);
    // negate every other sample to shift up to Nyquist
    // to maintain proper delay, need to make t=0 modulated with
    // cos(0),sin(0)
    int start = taps/2+1;
    // manually set start with integer math (rather than using ceil or floor)
    // since use of the -no-prec-div compiler flag can cause innacurate results
    // for odd filter lengths and very low values of xdelta
    if(taps%2!=0) start = (taps-1)/2+1;
    for(int i=start; i<taps; i+=2) filtOut[i]= -filtOut[i];
    // for odd length filters, floor(delay/xdelta) will be tap t=0 which
    // should not be negated, for even length filters, floor(delay/xdelta) will
    // be last tap before t=0, which we want to be negated
    start = (taps/2)-1;
    if(taps%2!=0) start = (taps-1)/2-1;
    for(int i=start; i>=0; i-=2) filtOut[i]= -filtOut[i];
    return delay;
  }
  
  // generate a Highpass filter with complex-valued taps
  // this method just provided for consistency, shouldn't need it in most
  // practical applications as it is still a purely real filter, just stored
  // in a complex buffer with imaginary part of each element set to 0
  //    taps = impulse response length in samples
  //    cutoff = HPF cutoff frequency in Hz
  //    xdelta = filter sampling interval in seconds
  //    win = window to apply to impulse response
  //    filtOut = filter impulse response
  //    GenFilter::createHPF = filter delay in seconds
  template <typename T>
  inline double createHPF(int taps, double cutoff, double xdelta,
      WinType win, vector<complex<float>,T> &filtOut)
  {
    vector<float> reTaps;
    double delay = createHPF(taps,cutoff,xdelta,win,reTaps);
    filtOut.resize(reTaps.size());
    for(int i=0; i<reTaps.size(); ++i)
    {
      filtOut[i].real(reTaps[i]);
      filtOut[i].imag(0);
    }
    return delay;
  }
  
  // generate a lowpass square root raised cosine pulse shaping filter with
  // real-value taps, the time variable t/Ts from the SRRC equation is
  // equivalent to current tap*xdelta*2*cutoff
  //    taps = impulse response length in samples
  //    cutoff = 1/2 the 3dB bandwidth of the SRRC (1/2 baudrate of signal to
  //             be pulse shaped)
  //    xdelta = filter sampling interval in seconds
  //    beta = SRRC rolloff factor
  //    filtOut = filter impulse response
  //    GenFilter::createLSRRC = filter delay in seconds
  template <typename T>
  inline double createLSRRC(int taps, double cutoff, double xdelta, double beta,
      vector<float,T> &filtOut)
  {
    // check validity of cutoff and xdelta
    if(cutoff > 1/(2*xdelta) || cutoff < 0)
    {
      throw std::logic_error("[GenFilter::createLSRRC] Invalid cutoff frequency"
          " specified");
    }
    // check validity of beta
    if(beta < 0 || beta > 1)
    {
      throw std::logic_error("[GenFilter::createLSRRC] Invalid SRRC rolloff "
          "factor specified");
    }
    double delay = double(taps-1)/2.0;  // starting tap sample index
    double tTs;         // time value t/Ts
    double edge = 0.0;  // time value that would result in divide by 0
    if(beta > 0) edge = 1/(4*beta);
    filtOut.resize(taps);
    // compute SRRC taps, avoiding divide by 0
    for(int i=0; i<taps; ++i)
    {
      tTs = (i-delay)*xdelta*2*cutoff;
      if(tTs == 0)  // center tap (t=0)
      {
        filtOut[i] = 1-beta+4*beta/pi_d;
      }
      else if(fabs(tTs) == edge) // divide by zero in SRRC equation
      {
        filtOut[i] = beta/sqrt(2.0)*((1+2/pi_d)*sin(pi_d*edge)+
            (1-2/pi_d)*cos(pi_d*edge));
      }
      else  // all other taps
      {
        filtOut[i] = (sin(pi_d*tTs*(1-beta))+
            4*beta*tTs*cos(pi_d*tTs*(1+beta)))/
            (pi_d*tTs*(1-(16*beta*beta*tTs*tTs)));
      }
    }
    return delay*xdelta;  // convert from sample to seconds
  }
  
  // generate a lowpass square root raised cosine pulse shaping filter with
  // complex-value taps, this method just provided for consistency, shouldn't
  // need it in most practical applications as it is still a purely real
  // filter, just stored in a complex buffer with imaginary part of each
  // element set to 0
  //    taps = impulse response length in samples
  //    cutoff = 1/2 the 3dB bandwidth of the SRRC (1/2 baudrate of signal to
  //             be pulse shaped)
  //    xdelta = filter sampling interval in seconds
  //    beta = SRRC rolloff factor
  //    filtOut = filter impulse response
  //    GenFilter::createLSRRC = filter delay in seconds
  template <typename T>
  inline double createLSRRC(int taps, double cutoff, double xdelta, double beta,
      vector<complex<float>,T> &filtOut)
  {
    // generate a real valued filter, then set imaginary part of each tap to 0
    vector<float> reTaps;
    double delay = createLSRRC(taps,cutoff,xdelta,beta,reTaps);
    filtOut.resize(reTaps.size());
    for(int i=0; i<reTaps.size(); ++i)
    {
      filtOut[i].real(reTaps[i]);
      filtOut[i].imag(0);
    }
    return delay;
  }
  
  // generate a bandpass square root raised cosine pulse shaping filter with
  // real-value taps, the time variable t/Ts from the SRRC equation is
  // equivalent to current tap*xdelta*(cutoff2-cutoff1), this filter is
  // generated by creating a prototype lowpass SRRC and modulating it up
  // to the desired frequency with a cosine
  //    taps = impulse response length in samples
  //    cutoff = 1/2 the 3dB bandwidth of the SRRC (1/2 baudrate of signal to
  //             be pulse shaped)
  //    xdelta = filter sampling interval in seconds
  //    beta = SRRC rolloff factor
  //    filtOut = filter impulse response
  //    GenFilter::createLSRRC = filter delay in seconds
  template <typename T>
  inline double createBSRRC(int taps, double cutoff1, double cutoff2,
      double xdelta, double beta, vector<float,T> &filtOut)
  {
    // check validity of cutoff and xdelta
    if(cutoff2 > 1/(2*xdelta) || cutoff2 <= cutoff1 || cutoff1 < 0)
    {
      throw std::logic_error("[GenFilter::createBSRRC] Invalid cutoff "
          "frequency specified");
    }
    // check validity of beta
    if(beta < 0 || beta > 1)
    {
      throw std::logic_error("[GenFilter::createBSRRC] Invalid SRRC rolloff "
          "factor specified");
    }
    // setup filter parameters, LPF cutoff will be 1/2 BPF bandwidth
    double delay = createLSRRC(taps,(cutoff2-cutoff1)/2,xdelta,beta,filtOut);
    // modulate filter to desired CF with cos so center of passband has phase=0
    double phase = 0;
    // 2*pi*Ts*Fc (2's cancel)
    double phaseIncr = pi_d*xdelta*(cutoff2+cutoff1);
    // to maintain proper delay, need to make t=0 modulated with cos(0),
    // fmod(delay,xdelta) should work but was not giving the right answer
    // when compiled with certain optimization flags, so instead a direct
    // even/odd filter tap check is used.  Also have to set start using
    // integer math in a branch rather than ceil/floor as agressive compiler
    // optimizations like --no-prec-div may cause the ceil and floor methods
    // to return the wrong answer for certain odd values of taps (ex. 255)
    int start = taps/2;
    if(taps%2 == 0) phase = phaseIncr/2;
    else
    {
      phase = 0;
      start = (taps-1)/2;
    }
    for(int i=start; i<taps; ++i)
    {
      filtOut[i] *= 2*cos(phase); // scale by 2 for LPF -> BPF conversion
      phase += phaseIncr;
    }
    if(taps%2 == 0)
    {
      phase = -phaseIncr/2;
      --start; // need to start just prior to t=0
    }
    else phase = 0;
    for(int i=start; i>=0; --i)
    {
      filtOut[i] *= 2*cos(phase); // scale by 2 for LPF -> BPF conversion
      phase -= phaseIncr;
    }
    if(taps%2!=0) filtOut[(taps-1)/2] *= .5; // t=0 tap got scaled twice
    return delay;
  }
  
  // generate a bandpass square root raised cosine pulse shaping filter with
  // real-value taps, the time variable t/Ts from the SRRC equation is
  // equivalent to current tap*xdelta*(cutoff2-cutoff1), this filter is
  // generated by creating a prototype lowpass SRRC and modulating it up
  // to the desired frequency with either a complex exponential (if
  // singleSided==true) or a cosine (if singleSided==false)
  //    taps = impulse response length in samples
  //    cutoff = 1/2 the 3dB bandwidth of the SRRC (1/2 baudrate of signal to
  //             be pulse shaped)
  //    xdelta = filter sampling interval in seconds
  //    beta = SRRC rolloff factor
  //    filtOut = filter impulse response
  //    singleSided = whether to create a single sided frequency response
  //    GenFilter::createLSRRC = filter delay in seconds
  template <typename T>
  inline double createBSRRC(int taps, double cutoff1, double cutoff2,
      double xdelta, double beta, vector<complex<float>,T> &filtOut,
      bool singleSided=true)
  {
    // check validity of cutoff and xdelta
    if(cutoff2 > 1/(2*xdelta) || cutoff2 <= cutoff1 || cutoff1 < -1/(2*xdelta))
    {
      throw std::logic_error("[GenFilter::createBSRRC] Invalid cutoff "
          "frequency specified");
    }
    // check validity of beta
    if(beta < 0 || beta > 1)
    {
      throw std::logic_error("[GenFilter::createBSRRC] Invalid SRRC rolloff "
          "factor specified");
    }
    // check validity of cutoff and xdelta
    double delay = createLSRRC(taps,(cutoff2-cutoff1)/2,xdelta,beta,filtOut);
    // modulate filter to desired CF, center of passband will have phase=0
    double phase = 0;
    // 2*pi*Ts*Fc (2's cancel)
    double phaseIncr = pi_d*xdelta*(cutoff2+cutoff1);
    int start = taps/2;
    if(singleSided) // modulate with complex sinusoid
    {
      // to maintain proper delay, need to make t=0 modulated with cos(0),
      // fmod(delay,xdelta) should work but was not giving the right answer
      // when compiled with certain optimization flags, so instead a direct
      // even/odd filter tap check is used.  Also have to set start using
      // integer math in a branch rather than ceil/floor as agressive compiler
      // optimizations like --no-prec-div may cause the ceil and floor methods
      // to return the wrong answer for certain odd values of taps (ex. 255)
      if(taps%2 == 0) phase = phaseIncr/2;
      else
      {
        phase = 0;
        start = (taps-1)/2;
      }
      for(int i=start; i<taps; ++i)
      {
        filtOut[i] *= complex<float>(cos(phase),sin(phase));
        phase += phaseIncr;
      }
      if(taps%2 == 0)
      {
        phase = -phaseIncr/2;
        --start; // need to start just prior to t=0
      }
      else phase = 0;
      for(int i=start; i>=0; --i)
      {
        filtOut[i] *= complex<float>(cos(phase),sin(phase));
        phase -= phaseIncr;
      }
    }
    else  // double sided, modulate with cosine
    {
      // to maintain proper delay, need to make t=0 modulated with cos(0),
      // fmod(delay,xdelta) should work but was not giving the right answer
      // when compiled with certain optimization flags, so instead a direct
      // even/odd filter tap check is used.  Also have to set start using
      // integer math in a branch rather than ceil/floor as agressive compiler
      // optimizations like --no-prec-div may cause the ceil and floor methods
      // to return the wrong answer for certain odd values of taps (ex. 255)
      if(taps%2 == 0) phase = phaseIncr/2;
      else
      {
        phase = 0;
        start = (taps-1)/2;
      }
      for(int i=start; i<taps; ++i)
      {
        // scale by 2 for LPF -> BPF conversion
        filtOut[i].real(filtOut[i].real()*2*cos(phase));
        phase += phaseIncr;
      }
      if(taps%2 == 0)
      {
        phase = -phaseIncr/2;
        --start; // need to start just prior to t=0
      }
      else phase = 0;
      for(int i=start; i>=0; --i)
      {
        // scale by 2 for LPF -> BPF conversion
        filtOut[i].real(filtOut[i].real()*2*cos(phase));
        phase -= phaseIncr;
      }
      // t=0 tap got scaled twice
      if(taps%2!=0) filtOut[(taps-1)/2].real(filtOut[(taps-1)/2].real()*.5);
    }
    return delay;
  }
  
  // generate a highpass square root raised cosine pulse shaping filter with
  // real-value taps, this filter is generating a prototype lowpass SRRC
  // then negating alternate taps to modulate it up to the nyquist rate
  //    taps = impulse response length in samples
  //    cutoff = 3 dB cutoff frequency in Hz
  //    xdelta = filter sampling interval in seconds
  //    beta = SRRC rolloff factor
  //    filtOut = filter impulse response
  //    GenFilter::createHSRRC = filter delay in seconds
  template <typename T>
  inline double createHSRRC(int taps, double cutoff, double xdelta, double beta,
      vector<float,T> &filtOut)
  {
    // check validity of cutoff and xdelta
    if(cutoff > 1/(2*xdelta) || cutoff < 0)
    {
      throw std::logic_error("[GenFilter::createHSRRC] Invalid cutoff "
          "frequency specified");
    }
    // check validity of beta
    if(beta < 0 || beta > 1)
    {
      throw std::logic_error("[GenFilter::createHSRRC] Invalid SRRC rolloff "
          "factor specified");
    }
    // generate prototype lowpass SRRC, reflect cutoff freq around Fs/4
    double delay = createLSRRC(taps,1.0/(2*xdelta)-cutoff,xdelta,beta,filtOut);
    // negate every other sample to shift up to Nyquist
    // to maintain proper delay, need to make t=0 modulated with
    // cos(0),sin(0)
    // have to set start using integer math in a branch rather than ceil/floor
    // as agressive compiler optimizations like --no-prec-div may cause the
    // ceil and floor methods to return the wrong answer for certain odd values
    // of taps (ex. 255)
    int start = taps/2+1;
    if(taps%2!=0) start = (taps-1)/2+1;
    for(int i=start; i<taps; i+=2) filtOut[i]= -filtOut[i];
    // for odd length filters, floor(delay/xdelta) will be tap t=0 which
    // should not be negated, for even length filters, floor(delay/xdelta) will
    // be last tap before t=0, which we want to be negated
    start = (taps/2)-1;
    if(taps%2!=0) start = (taps-1)/2-1;
    for(int i=start; i>=0; i-=2) filtOut[i]= -filtOut[i];
    return delay;
  }
  
  // generate a highpass square root raised cosine pulse shaping filter with
  // complex-value taps, this method just provided for consistency, shouldn't
  // need it in most practical applications as it is still a purely real
  // filter, just stored in a complex buffer with imaginary part of each
  // element set to 0
  //    taps = impulse response length in samples
  //    cutoff = 3 dB cutoff frequency in Hz
  //    xdelta = filter sampling interval in seconds
  //    beta = SRRC rolloff factor
  //    filtOut = filter impulse response
  //    GenFilter::createHSRRC = filter delay in seconds
  template <typename T>
  inline double createHSRRC(int taps, double cutoff, double xdelta, double beta,
      vector<complex<float>,T> &filtOut)
  {
    // generate a real valued filter, then set imaginary part of each tap to 0
    vector<float> reTaps;
    double delay = createHSRRC(taps,cutoff,xdelta,beta,reTaps);
    filtOut.resize(reTaps.size());
    for(int i=0; i<reTaps.size(); ++i)
    {
      filtOut[i].real(reTaps[i]);
      filtOut[i].imag(0);
    }
    return delay;
  }
  
  // Generate a lowpass interpolation filter with real-valued taps. Time domain
  // taps are computed directly as a sinc then windowed.  The taps may
  // optionally be time-reversed so convolution can be computed with a
  // time domain dot product.  The method returns the integer sample delay
  // of the filter (i.e. the number of samples required after your
  // interpolation reference sample required for interpolation).  For example,
  // if you'd like to interpolate the sample k+phi where k is the integer
  // sample index and phi is the interpolation phase in cycles (fractional
  // sample index), this method will return the number of samples, n, after
  // sample k that are required for the interpolation.  This can be used to
  // compute the index in the input sequence, d, where the interpolation dot
  // product should start via:
  //    d = k-(taps-1-n)
  // Thus the dot product (x[d] x[d+1] ... x[d+taps-1])*filtOut will interpolate
  // input sequence x at sample k+phi.
  //    taps = impulse response length in samples, equivalent to number of
  //        periods in the sinc if BW=1
  //    interpPhase = interpolation phase in cycles, should be between -1 and 1
  //    bw = bandwidth of the interpolation filter as a fraction of the sampling
  //        rate, i.e. .9 is equivalent to a BW of 90% of the input rate, this
  //        specifies the 3dB cutoff point, if no additional filtering is
  //        desired beyond interpolation, this should be set to 1.0
  //    win = window to apply to impulse response
  //    reverse = whether or not to create the filter time-reversed, set to
  //        true if convolution will be done via time domain dot-product
  //    repeat = whether or not to put 2 copies of each tap consecutively in
  //        filtOut, this may be useful for computing efficient complex*real
  //        dot products using SSE optimizations, setting this input to true
  //        will result in filtOut being size 2*taps
  //    filtOut = filter impulse response
  //    GenFilter::createInterpLPF = integer sample group delay to achieve this
  //        interpolation
  template <typename T>
  inline int createInterpFilt(int taps, double interpPhase, double bw,
      WinType win, bool reverse, bool repeat, vector<float,T> &filtOut)
  {
    int delay;    // integer sample delay to achieve desired interpolation phase
    double phase; // phase in cycles of first sample of the interp sinc
    double x;     // phase of current interp tap in radians
    const double smallNum = 1e-10;  // set taps smaller than this to 0
    vector<float> vtmp(taps);
    vector<float> vwin(taps);
    
    if(repeat) filtOut.resize(2*taps);
    else filtOut.resize(taps);
    
    // filter must have at least 2 taps or its not interpolating, just scaling
    if(taps < 2)
    {
      throw std::logic_error("[GenFilter::createInterpFilt] Number of taps to "
          "small, filter must have at least 2 taps.");
    }
    
    bw = fabs(bw); // bandwidth can not be negative
    
    // if number of taps is odd, the interpolation filters will be periodic
    // about the range [-.5 : .5) cycles and anything outside this range will
    // involve shifting the delay forward or back and using the corresponding
    // filter from that range
    if(taps%2)
    {
      delay = (taps-1)/2+int(floor(interpPhase+.5));
    }
    // if number of taps is even, the interpolation filters will be periodic
    // about the range [0 : 1) cycles and anything outside this range will
    // involve shifting the delay forward or back and using the corresponding
    // filter from that range
    else
    {
      delay = taps/2+int(floor(interpPhase));
    }
    // delay specifies the starting period of the sinc for interpolation and
    // interpPhase is the fractional offset into that starting period
    phase = -delay+interpPhase;
    
    // each sample in the filter will move the phase forward by one cycle
    for(int i=0; i<taps; ++i)
    {
      x = pi_d*(phase+i);
      if(x == 0) // sinc(0) = 1
      {
        vtmp[i] = bw*windowTap(0,win);
      }
      else // all other taps
      {
        vtmp[i] = (sin(bw*x)/x)*windowTap((phase+i)/taps,win);
      }
      // set taps that are really small to 0
      if(fabs(vtmp[i]) < smallNum) vtmp[i] = 0;
    }
    // time reverse and repeat if desired
    if(repeat && reverse)
    {
      for(int i=0; i<taps; ++i)
      {
        filtOut[2*i] = vtmp[taps-1-i];
        filtOut[2*i+1] = filtOut[2*i];
      }
    }
    else if(repeat)
    {
      for(int i=0; i<taps; ++i)
      {
        filtOut[2*i] = vtmp[i];
        filtOut[2*i+1] = filtOut[2*i];
      }
    }
    else if(reverse)
    {
      for(int i=0; i<taps; ++i) filtOut[i] = vtmp[taps-1-i];
    }
    else
    {
      for(int i=0; i<taps; ++i) filtOut[i] = vtmp[i];
    }
    return delay;
  }
  
  // Generate a bandpass interpolation filter with real-valued taps. Time
  // domain interpolation filter taps are computed directly as a sinc (via
  // createInterpFilt) then windowed and modulated to the desired center
  // frequency with a cosine.  See createInterpFilt for an explanation of the
  // delay return value and how to compute dot-product index from this value.
  //    taps = impulse response length in samples
  //    interpPhase = interpolation phase in cycles, should be between -1 and 1
  //    cutoff1 = 1st BPF cutoff frequency in Hz
  //    cutoff2 = 2nd BPF cutoff frequency in Hz
  //    xdelta = sampling interval of filter
  //    win = window to apply to impulse response
  //    reverse = whether or not to create the filter time-reversed, set to
  //        true if convolution will be done via time domain dot-product
  //    repeat = whether or not to put 2 copies of each tap consecutively in
  //        filtOut, this may be useful for computing efficient complex*real
  //        dot products using SSE optimizations, setting this input to true
  //        will result in filtOut being size 2*taps
  //    filtOut = filter impulse response
  //    GenFilter::createInterpBPF = integer sample group delay to achieve this
  //        interpolation
  template <typename T>
  inline int createInterpBPF(int taps, double interpPhase, double cutoff1,
      double cutoff2, double xdelta, WinType win, bool reverse, bool repeat,
      vector<float,T> &filtOut)
  {
    // check validity of cutoff and xdelta
    if(cutoff2 > 1/(2*xdelta) || cutoff2 <= cutoff1 || cutoff1 < 0)
    {
      throw std::logic_error("[GenFilter::createInterpBPF] Invalid cutoff "
          "frequency specified");
    }
    // setup filter parameters, interp LPF BW is given as a fraction of
    // sampling rate so need to convert our cutoffs from Hz
    int delay = createInterpFilt(taps,interpPhase,(cutoff2-cutoff1)*xdelta,win,
        reverse,repeat,filtOut);
    // generate cosine to modulate filter, repeat and reverse taps as necessary
    // 2*pi*Ts*Fc (2's cancel)
    double phaseIncr = pi_d*xdelta*(cutoff2+cutoff1);
    double phase = (-delay+interpPhase)*phaseIncr;
    if(reverse && repeat)
    {
      for(int i=0; i<taps; ++i)
      {
        filtOut[2*(taps-1-i)] *= 2*cos(phase); // scale by 2 for LPF -> BPF conversion
        filtOut[2*(taps-1-i)+1] = filtOut[2*(taps-1-i)];
        phase += phaseIncr;
      }
    }
    else if(reverse)
    {
      for(int i=0; i<taps; ++i)
      {
        filtOut[taps-1-i] *= 2*cos(phase); // scale by 2 for LPF -> BPF conversion
        phase += phaseIncr;
      }
    }
    else if(repeat)
    {
      for(int i=0; i<2*taps; i+=2)
      {
        filtOut[i] *= 2*cos(phase); // scale by 2 for LPF -> BPF conversion
        filtOut[i+1] *= filtOut[i];
        phase += phaseIncr;
      }
    }    
    else
    {
      for(int i=0; i<taps; ++i)
      {
        filtOut[i] *= 2*cos(phase); // scale by 2 for LPF -> BPF conversion
        phase += phaseIncr;
      }
    }
    return delay;
  }
  
  // Generate a bandpass interpolation filter with complex-valued taps. Time
  // domain interpolation filter taps are computed directly as a sinc (via
  // createInterpFilt) then windowed and modulated to the desired center
  // frequency with a complex exponential (if singleSided==true) or a cosine
  // (if singleSided==false).  See createInterpFilt for an explanation of the
  // delay return value and how to compute dot-product index from this value.
  //    taps = impulse response length in samples
  //    interpPhase = interpolation phase in cycles, should be between -1 and 1
  //    cutoff1 = 1st BPF cutoff frequency in Hz
  //    cutoff2 = 2nd BPF cutoff frequency in Hz
  //    xdelta = sampling interval of filter
  //    win = window to apply to impulse response
  //    reverse = whether or not to create the filter time-reversed, set to
  //        true if convolution will be done via time domain dot-product
  //    filtOut = filter impulse response
  //    singleSided = whether to create a single sided frequency response
  //    GenFilter::createInterpBPF = integer sample group delay to achieve this
  //        interpolation
  template <typename T>
  inline int createInterpBPF(int taps, double interpPhase, double cutoff1,
      double cutoff2, double xdelta, WinType win, bool reverse,
      vector<complex<float>,T> &filtOut, bool singleSided=true)
  {
    // check validity of cutoff and xdelta
    if(cutoff2 > 1/(2*xdelta) || cutoff2 < cutoff1 || cutoff1 < -1/(2*xdelta))
    {
      throw std::logic_error("[GenFilter::createBPF] Invalid cutoff frequency "
          "specified");
    }
    // setup filter parameters, interp LPF BW is given as a fraction of
    // sampling rate so need to convert our cutoffs from Hz
    vector<float> vtmp;
    int delay = createInterpFilt(taps,interpPhase,(cutoff2-cutoff1)*xdelta,win,
        false,false,vtmp);
    filtOut.resize(vtmp.size());
    // generate cosine to modulate filter, repeat and reverse taps as necessary
    // 2*pi*Ts*Fc (2's cancel)
    double phaseIncr = pi_d*xdelta*(cutoff2+cutoff1);
    double phase = (-delay+interpPhase)*phaseIncr;
    if(singleSided)
    {
      if(reverse)
      {
        for(int i=0; i<taps; ++i)
        {
          filtOut[taps-1-i] = vtmp[i]*complex<float>(cos(phase),sin(phase));
          phase += phaseIncr;
        }
      }
      else
      {
        for(int i=0; i<taps; ++i)
        {
          filtOut[i] = vtmp[i]*complex<float>(cos(phase),sin(phase));
          phase += phaseIncr;
        }
      }
    }
    else
    {
      if(reverse)
      {
        for(int i=0; i<taps; ++i)
        {
          // scale by 2 for LPF -> BPF conversion
          filtOut[taps-1-i].real(2*cos(phase)*vtmp[i]);
          filtOut[taps-1-i].imag(0);
          phase += phaseIncr;
        }
      }
      else
      {
        for(int i=0; i<taps; ++i)
        {
          // scale by 2 for LPF -> BPF conversion
          filtOut[i].real(2*cos(phase)*vtmp[i]);
          filtOut[i].imag(0);
          phase += phaseIncr;
        }
      }
    }
    return delay;
  }
  
  // Generate a lowpass interpolation filter with real-valued taps. Time domain
  // taps are computed directly as square-root-raised-cosine (SRRC).  The taps
  // may optionally be time-reversed so convolution can be computed with a
  // time domain dot product.  The method returns the integer sample delay
  // of the filter (i.e. the number of samples required after your
  // interpolation reference sample required for interpolation).  For example,
  // if you'd like to interpolate the sample k+phi where k is the integer
  // sample index and phi is the interpolation phase in cycles (fractional
  // sample index), this method will return the number of samples, n, after
  // sample k that are required for the interpolation.  This can be used to
  // compute the index in the input sequence, d, where the interpolation dot
  // product should start via:
  //    d = k-(taps-1-n)
  // Thus the dot product (x[d] x[d+1] ... x[d+taps-1])*filtOut will interpolate
  // input sequence x at sample k+phi.
  //    taps = impulse response length in samples, equivalent to number of
  //        periods in the SRRC if BW=1
  //    interpPhase = interpolation phase in cycles, should be between 1 and -1
  //    bw = bandwidth of the interpolation filter as a fraction of the sampling
  //        rate, i.e. .9 is equivalent to a BW of 90% of the input rate, this
  //        specifies the 3dB cutoff point, if no additional filtering is
  //        desired beyond interpolation, this should be set to 1.0
  //    beta = SRRC rolloff factor
  //    reverse = whether or not to create the filter time-reversed, set to
  //        true if convolution will be done via time domain dot-product
  //    repeat = whether or not to put 2 copies of each tap consecutively in
  //        filtOut, this may be useful for computing efficient complex*real
  //        dot products using SSE optimizations, setting this input to true
  //        will result in filtOut being size 2*taps
  //    filtOut = filter impulse response
  //    GenFilter::createInterpLSRRC = integer sample group delay to achieve
  //        this interpolation
  template <typename T>
  inline double createInterpSRRC(int taps, double interpPhase, double bw,
      double beta, bool reverse, bool repeat, vector<float,T> &filtOut)
  {
    int delay;    // integer sample delay to achieve desired interpolation phase
    double phase; // phase in cycles of first sample of the interp sinc
    double x;     // phase of current interp tap in radians
    double edge=0;// checks for edge cases where denominator of SRRC would be 0
    const double smallNum = 1e-10;  // set taps smaller than this to 0
    vector<float> vtmp(taps);
    vector<float> vwin(taps);
    
    if(repeat) filtOut.resize(2*taps);
    else filtOut.resize(taps);
    
    // filter must have at least 2 taps or its not interpolating, just scaling
    if(taps < 2)
    {
      throw std::logic_error("[GenFilter::createInterpFilt] Number of taps to "
          "small, filter must have at least 2 taps.");
    }
    
    bw = fabs(bw); // bandwidth can not be negative
    
    // check to make sure that rolloff factor beta is in bounds
    if(beta<0)
    {
      beta = 0;
      cerr<<"[GenFilter::createInterpFilt] WARNING: Requested rolloff factor "
          <<"(beta) is negative, instead using beta="<<beta<<endl;
    }
    else if(beta > 1)
    {
      beta = 1;
      cerr<<"[GenFilter::createInterpFilt] WARNING: Requested rolloff factor "
          <<"(beta) is >1, instead using beta="<<beta<<endl;
    }
 
    // avoid divide by 0
    if(beta>0) edge = 1/(4*beta);
    
    // if number of taps is odd, the interpolation filters will be periodic
    // about the range [-.5 : .5) cycles and anything outside this range will
    // involve shifting the delay forward or back and using the corresponding
    // filter from that range
    if(taps%2)
    {
      delay = (taps-1)/2+floor(interpPhase+.5);
    }
    // if number of taps is even, the interpolation filters will be periodic
    // about the range [0 : 1) cycles and anything outside this range will
    // involve shifting the delay forward or back and using the corresponding
    // filter from that range
    else
    {
      delay = taps/2+floor(interpPhase);
    }
    // delay specifies the starting period of the sinc for interpolation and
    // interpPhase is the fractional offset into that starting period
    phase = -delay+interpPhase;
    
    // each sample in the filter will move the phase forward by one cycle
    for(int i=0; i<taps; ++i)
    {
      x = (phase+i)*bw;
      if(x == 0) // center tap at t=0 (zero denominator)
      {
        vtmp[i] = 1-beta+4*beta/pi_d;
      }
      else if(abs(x)==edge)   // edge case resulting in zero denominator
      {
        vtmp[i] = beta/sqrt(2.0)*((1+2/pi_d)*sin(pi_d*edge)+
            (1-2/pi_d)*cos(pi_d*edge));
      }
      else // all other taps
      {
        vtmp[i] = (sin(pi_d*x*(1-beta))+
            4*beta*x*cos(pi_d*x*(1+beta)))/
            (pi_d*x*(1-(16*beta*beta*x*x)));
      }
      // set taps that are really small to 0
      if(fabs(vtmp[i]) < smallNum) vtmp[i] = 0;
    }
    // time reverse and repeat if desired
    if(repeat && reverse)
    {
      for(int i=0; i<taps; ++i)
      {
        filtOut[2*i] = vtmp[taps-1-i];
        filtOut[2*i+1] = filtOut[2*i];
      }
    }
    else if(repeat)
    {
      for(int i=0; i<taps; ++i)
      {
        filtOut[2*i] = vtmp[i];
        filtOut[2*i+1] = filtOut[2*i];
      }
    }
    else if(reverse)
    {
      for(int i=0; i<taps; ++i) filtOut[i] = vtmp[taps-1-i];
    }
    else
    {
      for(int i=0; i<taps; ++i) filtOut[i] = vtmp[i];
    }
    return delay;
  }
  
  // convert string to window type, case insensitive, only first 4 letters
  // of input string will be parsed, if the input string does not match
  // a valid window type, GenFilter::UNKNOWN will be returned.
  //    winStr = string representation of tapering window type (i.e. "BH92")
  //    GenFilter::str2win = WinType enumeration of tapering window
  inline WinType str2win(string winStr)
  {
    std::transform(winStr.begin(),winStr.end(),winStr.begin(),::toupper);
    if(winStr.find("NONE") == 0) return NONE;
    else if(winStr.find("HANN") == 0) return HANN;
    else if(winStr.find("HAMM") == 0) return HAMM;
    else if(winStr.find("BH61") == 0) return BH61;
    else if(winStr.find("BH67") == 0) return BH67;
    else if(winStr.find("BH74") == 0) return BH74;
    else if(winStr.find("BH92") == 0) return BH92;
    else if(winStr.find("SRRC") == 0) return SRRC;
    return UNKNOWN;
  }
  
  // convert a GenFilter::WinType enumeration value to its string
  // representation.  If win enumeration value is not in the list,
  // "UNKNOWN" will be returned.  Returned strings are the 4 letter
  // abbreviation of the window type in all capital letters.
  //    win = window type enumeration value
  //    GenFilter::win2str = string representation of window type
  inline string win2str(WinType win)
  {
    if(win == NONE) return "NONE";
    else if(win == HANN) return "HANN";
    else if(win == HAMM) return "HAMM";
    else if(win == BH61) return "BH61";
    else if(win == BH67) return "BH67";
    else if(win == BH74) return "BH74";
    else if(win == BH92) return "BH92";
    else if(win == SRRC) return "SRRC";
    return "UNKNOWN";
  }
  
  // returns the tap of the desired window at a specific phase (delay).  The
  // phase should be specified between -.5 and .5 inclusive where -.5 and .5
  // would be the points where the window value approaches 0 on either end
  // (since these are tapering windows)
  //    phase = phase in cycles between -.5 and .5 to generate window tap for
  //    win = tapering window function
  //    GenFilter::windowTap = value of the window at the specified phase
  inline double windowTap(const double phase, const WinType win)
  {
    double x = two_pi_d*phase; // convert to radians
    if(win == NONE) return 1.0;
    else if(win == HAMM) return 0.54+0.46*cos(x);
    else if(win == HANN) return 0.5*(1+cos(x));
    else if(win == BH61) return 0.44959+0.49364*cos(x)+0.05677*cos(x*2);
    else if(win == BH67) return 0.42323+0.49755*cos(x)+0.07922*cos(x*2);
    else if(win == BH74) return 0.40217+0.49703*cos(x)+0.09392*cos(x*2)+0.00183*cos(x*3);
    else if(win == BH92) return 0.35875+0.48829*cos(x)+0.14128*cos(x*2)+0.01168*cos(x*3);
    else
    {
      throw std::logic_error("[GenFilter::createWindow] Window type not "
          "recognized");
    }
    return 0.0;
  }
      
  // creates windows
  //    vwin = vector to hold window
  //    len = number of taps create window over
  //    win = window type (NONE, HAMM, HANN, BH61, BH67, BH74, BH92)
  template <typename T>
  inline void createWindow(vector<double,T> &vwin, int len, WinType win)
  {
    vwin.resize(len);
    double ph = double(len-1)/2.0;
    double x;
    if(win == NONE)
    {
      for(int i=0; i<len; ++i) vwin[i] = 1;
    }
    else if(win == HAMM)
    {
      for(int i=0; i<len; ++i)
      {
        x = two_pi_d*(i-ph)/(len-1);
        vwin[i] = 0.54+0.46*cos(x);
      }
    }
    else if(win == HANN)
    {
      for(int i=0; i<len; ++i)
      {
        x = two_pi_d*(i-ph)/(len-1);
        vwin[i] = 0.5*(1+cos(x));
      }
    }
    else if(win == BH61)
    {
      for(int i=0; i<len; ++i)
      {
        x = two_pi_d*(i-ph)/(len-1);
        vwin[i] = 0.44959+0.49364*cos(x)+0.05677*cos(x*2);
      }
    }
    else if(win == BH67)
    {
      for(int i=0; i<len; ++i)
      {
        x = two_pi_d*(i-ph)/(len-1);
        vwin[i] = 0.42323+0.49755*cos(x)+0.07922*cos(x*2);
      }
    }
    else if(win == BH74)
    {
      for(int i=0; i<len; ++i)
      {
        x = two_pi_d*(i-ph)/(len-1);
        vwin[i] = 0.40217+0.49703*cos(x)+0.09392*cos(x*2)+0.00183*cos(x*3);
      }
    }
    else if(win == BH92)
    {
      for(int i=0; i<len; ++i)
      {
        x = two_pi_d*(i-ph)/(len-1);
        vwin[i] = 0.35875+0.48829*cos(x)+0.14128*cos(x*2)+0.01168*cos(x*3);
      }
    }
    else
    {
      throw std::logic_error("[GenFilter::createWindow] Window type not "
          "recognized");
    }
    return;
  }
  
  // creates windows, overloaded method to convert format to single-precision
  // floating point
  //    vwin = vector to hold window
  //    len = number of taps create window over
  //    win = window type (NONE, HAMM, HANN, BH61, BH67, BH74, BH92)
  template <typename T>
  inline void createWindow(vector<float,T> &vwin, int len, WinType win)
  {
    vector<double> vtmp;
    createWindow(vtmp,len,win);
    vwin.resize(len);
    for(int i=0; i<len; ++i) vwin[i] = float(vtmp[i]);
  }
  
} // namespace GenFilter

#endif // GEN_FILTER_H
