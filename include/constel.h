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
//==============================================================================
//  Name:     constel.h
//
//  Purpose:  Namespace for various PSK/QAM constellations and bit mappings
//
//  Created:  2013/03/17
//
//  Modifications:
//
//  Description:
//    Namespace containing various PSK and QAM constellations and support
//    structures for fast lookup and quantization.
//
//==============================================================================

#ifndef CONSTEL_H
#define CONSTEL_H

#include "complex_num.h"

namespace Constel
{
  // unit-energy symbol alphabets for various constellations, symbols are
  // stored in the array according to a grey coded bit mapping so that the
  // decimal representation of the bit map for each symbol is equal to the
  // array index for that symbol so the array can be used as a lookup table
  const Complex8 BPSK[2] = {
      Complex8(1,0),   // 0, 0b
      Complex8(-1,0)}; // pi, 1b

  const Complex8 QPSK[4] = {
      Complex8(.5*sqrt(2),.5*sqrt(2)),   // pi/4,  00b
      Complex8(.5*sqrt(2),-.5*sqrt(2)),  // 7pi/4, 01b
      Complex8(-.5*sqrt(2),.5*sqrt(2)),  // 3pi/4, 10b
      Complex8(-.5*sqrt(2),-.5*sqrt(2))};// 5pi/4, 11b

  const Complex8 PSK8[8] = {
      Complex8(cos(dpi/8),sin(dpi/8)),       // pi/8,   000b
      Complex8(cos(7*dpi/8),sin(7*dpi/8)),   // 7pi/8,  001b
      Complex8(cos(15*dpi/8),sin(15*dpi/8)), // 15pi/8, 010b
      Complex8(cos(9*dpi/8),sin(9*dpi/8)),   // 9pi/8,  011b
      Complex8(cos(3*dpi/8),sin(3*dpi/8)),   // 3pi/8,  100b
      Complex8(cos(5*dpi/8),sin(5*dpi/8)),   // 5pi/8,  101b
      Complex8(cos(13*dpi/8),sin(13*dpi/8)), // 13pi/8, 110b
      Complex8(cos(11*dpi/8),sin(11*dpi/8))};// 11pi/8, 111b

  const Complex8 QAM16[16] = {
      Complex8(.3*sqrt(10),.3*sqrt(10)),  // (+3,+3), 0000b
      Complex8(.3*sqrt(10),.1*sqrt(10)),  // (+3,+1), 0001b
      Complex8(.3*sqrt(10),-.3*sqrt(10)), // (+3,-3), 0010b
      Complex8(.3*sqrt(10),-.1*sqrt(10)), // (+3,-1), 0011b
 
      Complex8(.1*sqrt(10),.3*sqrt(10)),  // (+1,+3), 0100b
      Complex8(.1*sqrt(10),.1*sqrt(10)),  // (+1,+1), 0101b
      Complex8(.1*sqrt(10),-.3*sqrt(10)), // (+1,-3), 0110b
      Complex8(.1*sqrt(10),-.1*sqrt(10)), // (+1,-1), 0111b
 
      Complex8(-.3*sqrt(10),.3*sqrt(10)),  // (-3,+3), 1000b
      Complex8(-.3*sqrt(10),.1*sqrt(10)),  // (-3,+1), 1001b
      Complex8(-.3*sqrt(10),-.3*sqrt(10)), // (-3,-3), 1010b
      Complex8(-.3*sqrt(10),-.1*sqrt(10)), // (-3,-1), 1011b
 
      Complex8(-.1*sqrt(10),.3*sqrt(10)),  // (-1,+3), 1100b
      Complex8(-.1*sqrt(10),.1*sqrt(10)),  // (-1,+1), 1101b
      Complex8(-.1*sqrt(10),-.3*sqrt(10)), // (-1,-3), 1110b
      Complex8(-.1*sqrt(10),-.1*sqrt(10))};// (-1,-1), 1111b

  const Complex8 QAM64[64] = {
      Complex8(7.0/sqrt(42.0),7.0/sqrt(42.0)), // (+7,+7), 000 000b
      Complex8(7.0/sqrt(42.0),5.0/sqrt(42.0)), // (+7,+5), 000 001b
      Complex8(7.0/sqrt(42.0),1.0/sqrt(42.0)), // (+7,+1), 000 010b
      Complex8(7.0/sqrt(42.0),3.0/sqrt(42.0)), // (+7,+3), 000 011b
      Complex8(7.0/sqrt(42.0),-3.0/sqrt(42.0)), // (+7,-3), 000 100b
      Complex8(7.0/sqrt(42.0),-5.0/sqrt(42.0)), // (+7,-5), 000 101b
      Complex8(7.0/sqrt(42.0),-1.0/sqrt(42.0)), // (+7,-1), 000 110b
      Complex8(7.0/sqrt(42.0),-7.0/sqrt(42.0)), // (+7,-7), 000 111b
 
      Complex8(5.0/sqrt(42.0),7.0/sqrt(42.0)), // (+5,+7), 001 000b
      Complex8(5.0/sqrt(42.0),5.0/sqrt(42.0)), // (+5,+5), 001 001b
      Complex8(5.0/sqrt(42.0),1.0/sqrt(42.0)), // (+5,+1), 001 010b
      Complex8(5.0/sqrt(42.0),3.0/sqrt(42.0)), // (+5,+3), 001 011b
      Complex8(5.0/sqrt(42.0),-3.0/sqrt(42.0)), // (+5,-3), 001 100b
      Complex8(5.0/sqrt(42.0),-5.0/sqrt(42.0)), // (+5,-5), 001 101b
      Complex8(5.0/sqrt(42.0),-1.0/sqrt(42.0)), // (+5,-1), 001 110b
      Complex8(5.0/sqrt(42.0),-7.0/sqrt(42.0)), // (+5,-7), 001 111b
 
      Complex8(1.0/sqrt(42.0),7.0/sqrt(42.0)), // (+1,+7), 010 000b
      Complex8(1.0/sqrt(42.0),5.0/sqrt(42.0)), // (+1,+5), 010 001b
      Complex8(1.0/sqrt(42.0),1.0/sqrt(42.0)), // (+1,+1), 010 010b
      Complex8(1.0/sqrt(42.0),3.0/sqrt(42.0)), // (+1,+3), 010 011b
      Complex8(1.0/sqrt(42.0),-3.0/sqrt(42.0)), // (+1,-3), 010 100b
      Complex8(1.0/sqrt(42.0),-5.0/sqrt(42.0)), // (+1,-5), 010 101b
      Complex8(1.0/sqrt(42.0),-1.0/sqrt(42.0)), // (+1,-1), 010 110b
      Complex8(1.0/sqrt(42.0),-7.0/sqrt(42.0)), // (+1,-7), 010 111b
 
      Complex8(3.0/sqrt(42.0),7.0/sqrt(42.0)), // (+3,+7), 011 000b
      Complex8(3.0/sqrt(42.0),5.0/sqrt(42.0)), // (+3,+5), 011 001b
      Complex8(3.0/sqrt(42.0),1.0/sqrt(42.0)), // (+3,+1), 011 010b
      Complex8(3.0/sqrt(42.0),3.0/sqrt(42.0)), // (+3,+3), 011 011b
      Complex8(3.0/sqrt(42.0),-3.0/sqrt(42.0)), // (+3,-3), 011 100b
      Complex8(3.0/sqrt(42.0),-5.0/sqrt(42.0)), // (+3,-5), 011 101b
      Complex8(3.0/sqrt(42.0),-1.0/sqrt(42.0)), // (+3,-1), 011 110b
      Complex8(3.0/sqrt(42.0),-7.0/sqrt(42.0)), // (+3,-7), 011 111b
 
      Complex8(-3.0/sqrt(42.0),7.0/sqrt(42.0)), // (-3,+7), 100 000b
      Complex8(-3.0/sqrt(42.0),5.0/sqrt(42.0)), // (-3,+5), 100 001b
      Complex8(-3.0/sqrt(42.0),1.0/sqrt(42.0)), // (-3,+1), 100 010b
      Complex8(-3.0/sqrt(42.0),3.0/sqrt(42.0)), // (-3,+3), 100 011b
      Complex8(-3.0/sqrt(42.0),-3.0/sqrt(42.0)), // (-3,-3), 100 100b
      Complex8(-3.0/sqrt(42.0),-5.0/sqrt(42.0)), // (-3,-5), 100 101b
      Complex8(-3.0/sqrt(42.0),-1.0/sqrt(42.0)), // (-3,-1), 100 110b
      Complex8(-3.0/sqrt(42.0),-7.0/sqrt(42.0)), // (-3,-7), 100 111b
 
      Complex8(-5.0/sqrt(42.0),7.0/sqrt(42.0)), // (-5,+7), 101 000b
      Complex8(-5.0/sqrt(42.0),5.0/sqrt(42.0)), // (-5,+5), 101 001b
      Complex8(-5.0/sqrt(42.0),1.0/sqrt(42.0)), // (-5,+1), 101 010b
      Complex8(-5.0/sqrt(42.0),3.0/sqrt(42.0)), // (-5,+3), 101 011b
      Complex8(-5.0/sqrt(42.0),-3.0/sqrt(42.0)), // (-5,-3), 101 100b
      Complex8(-5.0/sqrt(42.0),-5.0/sqrt(42.0)), // (-5,-5), 101 101b
      Complex8(-5.0/sqrt(42.0),-1.0/sqrt(42.0)), // (-5,-1), 101 110b
      Complex8(-5.0/sqrt(42.0),-7.0/sqrt(42.0)), // (-5,-7), 101 111b
 
      Complex8(-1.0/sqrt(42.0),7.0/sqrt(42.0)), // (-1,+7), 110 000b
      Complex8(-1.0/sqrt(42.0),5.0/sqrt(42.0)), // (-1,+5), 110 001b
      Complex8(-1.0/sqrt(42.0),1.0/sqrt(42.0)), // (-1,+1), 110 010b
      Complex8(-1.0/sqrt(42.0),3.0/sqrt(42.0)), // (-1,+3), 110 011b
      Complex8(-1.0/sqrt(42.0),-3.0/sqrt(42.0)), // (-1,-3), 110 100b
      Complex8(-1.0/sqrt(42.0),-5.0/sqrt(42.0)), // (-1,-5), 110 101b
      Complex8(-1.0/sqrt(42.0),-1.0/sqrt(42.0)), // (-1,-1), 110 110b
      Complex8(-1.0/sqrt(42.0),-7.0/sqrt(42.0)), // (-1,-7), 110 111b
 
      Complex8(-7.0/sqrt(42.0),7.0/sqrt(42.0)), // (-7,+7), 111 000b
      Complex8(-7.0/sqrt(42.0),5.0/sqrt(42.0)), // (-7,+5), 111 001b
      Complex8(-7.0/sqrt(42.0),1.0/sqrt(42.0)), // (-7,+1), 111 010b
      Complex8(-7.0/sqrt(42.0),3.0/sqrt(42.0)), // (-7,+3), 111 011b
      Complex8(-7.0/sqrt(42.0),-3.0/sqrt(42.0)), // (-7,-3), 111 100b
      Complex8(-7.0/sqrt(42.0),-5.0/sqrt(42.0)), // (-7,-5), 111 101b
      Complex8(-7.0/sqrt(42.0),-1.0/sqrt(42.0)), // (-7,-1), 111 110b
      Complex8(-7.0/sqrt(42.0),-7.0/sqrt(42.0))};// (-7,-7), 111 111b
  
  // define some bit masks for fast quantization of soft symbols
  const int PSK8_MASK1=1; // 001b
  const int PSK8_MASK2=2; // 010b
  const int PSK8_MASK3=4; // 100b
  const int PSK8_SHFT1=sizeof(int)*8-1; // place int sign bit at 001b
  const int PSK8_SHFT2=sizeof(int)*8-2; // place int sign bit at 010b
  const int PSK8_SHFT3=sizeof(int)*8-3; // place int sign bit at 100b
  
  // Given an input buffer of bits, stored as chars with 1 bit per byte,
  // maps the bits to indices for the specified constellation so they can be
  // used as a lookup into the symbol buffer to perform bit to symbol mapping.
  // If the number of input bits does not result in an integet number of
  // symbols, the number of unused bits is returned
  //    bits = buffer of bits, with each bit stored in the LSb of a char
  //    nbits = number of bits (i.e. number of chars) in bits buffer
  //    idx = output indices of the constellation symbol in the corresponding
  //        constellation lookup buffer, must have room for
  //        nbits/<bits per symbol> integers
  inline int bytesToBpskIdx(const char *bits, const int nbits, int *idx)
  {
    for(int i=0; i<nbits; ++i) idx[i] = bits[i];
    return 0;
  }
  
  inline int bytesToQpskIdx(const char *bits, const int nbits, int *idx)
  {
    int bps = 2;
    int nout = nbits/bps;
    int j = 0;
    for(int i=0; i<nout*bps; i+=bps)
    {
      idx[j++] = (bits[i]<<1)|(bits[i+1]);
    }
    return nbits-(bps*nout);
  }
  
  inline int bytesTo8PskIdx(const char *bits, const int nbits, int *idx)
  {
    int bps = 3;
    int nout = nbits/bps;
    int j = 0;
    for(int i=0; i<nout*bps; i+=bps)
    {
      idx[j++] = (bits[i]<<2)|(bits[i+1]<<1)|(bits[i+2]);
    }
    return nbits-(bps*nout);
  }
  
  inline int bytesTo16QamIdx(const char *bits, const int nbits, int *idx)
  {
    int bps = 4;
    int nout = nbits/bps;
    int j = 0;
    for(int i=0; i<nout*bps; i+=bps)
    {
      idx[j++] = (bits[i]<<3)|(bits[i+1]<<2)|(bits[i+2]<<1)|(bits[i+3]);
    }
    return nbits-(bps*nout);
  }
  inline int bytesTo64QamIdx(const char *bits, const int nbits, int *idx)
  {
    int bps = 6;
    int nout = nbits/bps;
    int j = 0;
    for(int i=0; i<nout*bps; i+=bps)
    {
      idx[j++] = (bits[i]<<5)|(bits[i+1]<<4)|(bits[i+2]<<3)|
          (bits[i+3]<<2)|(bits[i+4]<<1)|(bits[i+5]);
    }
    return nbits-(bps*nout);
  }
  
  // Given an input buffer of bits, stored as chars with 1 bit per byte,
  // maps the bits to symbols.  If the number of input bits does not result
  // in an integet number of symbols, the number of unused bits is returned
  //    bits = buffer of bits, with each bit stored in the LSb of a char
  //    nbits = number of bits (i.e. number of chars) in bits buffer
  //    syms = output constellation symbols, must have room for
  //        nbits/<bits per symbol> Complex8 samples
  inline int bytesToBpsk(const char *bits, const int nbits, Complex8 *syms)
  {
    for(int i=0; i<nbits; ++i) syms[i] = BPSK[bits[i]];
    return 0;
  }
  
  inline int bytesToQpsk(const char *bits, const int nbits, Complex8 *syms)
  {
    int bps = 2;
    int nout = nbits/bps;
    int j = 0;
    for(int i=0; i<nout*bps; i+=bps)
    {
      syms[j++] = QPSK[(bits[i]<<1)|(bits[i+1])];
    }
    return nbits-(bps*nout);
  }
  
  inline int bytesTo8Psk(const char *bits, const int nbits, Complex8 *syms)
  {
    int bps = 3;
    int nout = nbits/bps;
    int j = 0;
    for(int i=0; i<nout*bps; i+=bps)
    {
      syms[j++] = PSK8[(bits[i]<<2)|(bits[i+1]<<1)|(bits[i+2])];
    }
    return nbits-(bps*nout);
  }
  
  inline int bytesTo16Qam(const char *bits, const int nbits, Complex8 *syms)
  {
    int bps = 4;
    int nout = nbits/bps;
    int j = 0;
    for(int i=0; i<nout*bps; i+=bps)
    {
      syms[j++] = QAM16[(bits[i]<<3)|(bits[i+1]<<2)|(bits[i+2]<<1)|(bits[i+3])];
    }
    return nbits-(bps*nout);
  }
  inline int bytesTo64Qam(const char *bits, const int nbits, Complex8 *syms)
  {
    int bps = 6;
    int nout = nbits/bps;
    int j = 0;
    for(int i=0; i<nout*bps; i+=bps)
    {
      syms[j++] = QAM64[(bits[i]<<5)|(bits[i+1]<<4)|(bits[i+2]<<3)|
          (bits[i+3]<<2)|(bits[i+4]<<1)|(bits[i+5])];
    }
    return nbits-(bps*nout);
  }
  
} // namespace Constel
#endif // CONSTEL_H
