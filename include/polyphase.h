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
#include <math.h>
#include <string>
#include <stdexcept>
#include <complex>  // for Complex8 data type
#include "evm.h"    // for numerical math constant pi_d

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::complex;

#ifndef Complex8
#define Complex8 complex<float>
#endif

class Polyphase
{
  public:
    //==============//
    // Constructors //
    //==============//
    // Default constructor creates a small windowed sinc filter bank, the total
    // number of taps in the canonical (high rate) windowed sinc filter will be
    // tapsPerFilt*numFilters but may contain some 0-valued taps on the edges to
    // ensure the filter is centered on time T=0 and that one of the taps is at
    // time T=0
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
    inline Polyphase(const int numFilters=3, const int tapsPerFilt=3,
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
    inline Polyphase(const int numFilters, const int tapsPerFilt,
        const double bw, const double rolloff);

    // Constructor for generating a filter bank from a user defined canonical
    // filter
    //    filter = canonical filter taps to generate bank from
    //    ntaps = number of taps in the filter
    //    isCX = whether taps are real or complex valued (both are assumed to
    //        be single precision floating point)
    //    numFilters = number of filters to split the canonical filter into
    //    centerTap = index of the center or time t=0 tap, may be fractional
    //        if the specific tap corresponding to t=0 (i.e. phase of 0) is
    //        not in the filter, must be >= 0 and < numFilters*tapsPerFilt
    //    firstPhase = desired interpolation phase of the first filter in the
    //        bank, filters will be ordered in increasing phase increments of
    //	      1.0/numFilters cycles
    inline Polyphase(char const * const filter, const int ntaps, const bool isCX,
        const int numFilters, const double centerTap, const double firstPhase);

    // Constructor for generating a windowed sinc filter bank

    // Constructor for generating a SRRC filter bank  
}

#endif // POLYPHASE_H
