/* 
 * Copyright (c) 2017 Nick Xenias
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
//  Purpose:  Generates a polyphase filter bank from an FIR input filter
//
//  Created:  2017/04/03
// 
//  Description:
//    Generates a polyphase filter bank from an FIR input filter.  For various
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
//    HANN    -32        1.50          18       Hanning (cosine bell) 
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
#include "complex_num.h"  // for Complex8 data type
#include "evm.h" // for numerical math constant pi_d

using std::vector;
using std::string;
using std::cerr;
using std::endl;

class Polyphase
{
  public:
    //==============//
    // Constructors //
    //==============//
    // Default constructor creates a small windowed sinc filter bank
    inline Polyphase();

    // Constructor for generating a filter bank from a user defined filter,
    // the length of the filter is assumed to be numFilters*numTaps.
    //    filter = canonical filter taps to generate bank from, assumed to be
    //        length numFilters*tapsPerFilt
    //    isCX = whether taps are real or complex format
    //    numFilters = number of filters to split the canonical filter into
    //    tapsPerFilt = number of taps per split filter (i.e. polyphase leg)
    //    centerTap = index of the center or time t=0 tap, may be fractional
    //        if the specific tap corresponding to t=0 (i.e. phase of 0) is
    //        not in the filter, must be >= 0 and < numFilters*tapsPerFilt
    //    firstPhase = interpolation phase of the first filter in the bank,
    //        filters will be ordered in increasing phase increments of
    //	      1.0/numFilters cycles
    inline Polyphase(const char *filter, const bool isCX, const int numFilters,
        const int tapsPerFilt, const double centerTap, const double firstPhase);

    // Constructor for generating a windowed sinc filter bank

    // Constructor for generating a SRRC filter bank  
}

#endif // POLYPHASE_H
