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

#include <complex>
using std::complex;

#ifndef Complex8
#define Complex8 complex<float>
#endif

namespace Constel
{
  // unit-energy symbol alphabets for various constellations, symbols are
  // stored in the array according to a grey coded bit mapping so that the
  // decimal representation of the bit map for each symbol is equal to the
  // array index for that symbol so the array can be used as a lookup table
  //        |
  //  -- 1 -|- 0 --
  //        |
  static const Complex8 BPSK[2] = {
      Complex8(1,0),   // 0, 0b
      Complex8(-1,0)}; // pi, 1b

  //        |
  //   10   |   00
  //        |
  // -------|-------
  //        |
  //   11   |   01
  //        |
  static const Complex8 QPSK[4] = {
      Complex8(.5*sqrt(2),.5*sqrt(2)),   // pi/4,  00b
      Complex8(.5*sqrt(2),-.5*sqrt(2)),  // 7pi/4, 01b
      Complex8(-.5*sqrt(2),.5*sqrt(2)),  // 3pi/4, 10b
      Complex8(-.5*sqrt(2),-.5*sqrt(2))};// 5pi/4, 11b

  //            |
  //       101  |  100
  //            |
  //   001      |      000
  // -----------|-----------
  //   011      |      010
  //            |
  //       111  |  110
  //            |
  static const Complex8 PSK8[8] = {
      Complex8(cos(dpi/8),sin(dpi/8)),       // pi/8,   000b
      Complex8(cos(7*dpi/8),sin(7*dpi/8)),   // 7pi/8,  001b
      Complex8(cos(15*dpi/8),sin(15*dpi/8)), // 15pi/8, 010b
      Complex8(cos(9*dpi/8),sin(9*dpi/8)),   // 9pi/8,  011b
      Complex8(cos(3*dpi/8),sin(3*dpi/8)),   // 3pi/8,  100b
      Complex8(cos(5*dpi/8),sin(5*dpi/8)),   // 5pi/8,  101b
      Complex8(cos(13*dpi/8),sin(13*dpi/8)), // 13pi/8, 110b
      Complex8(cos(11*dpi/8),sin(11*dpi/8))};// 11pi/8, 111b

  //              |
  //  1000  1100  |  0100  0000
  //              |
  //  1001  1101  |  0101  0001
  // ---------------------------
  //  1011  1111  |  0111  0011
  //              |
  //  1010  1110  |  0110  0010
  //              |
  // Define decision boundary for 16-QAM
  const float QAM16_BOUND = .2*sqrt(10);  // boundary between syms +1 and +3
  static const Complex8 QAM16[16] = {
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

  //                              |
  //  111000 101000 100000 110000 | 010000 011000 001000 000000
  //                              |
  //  111001 101001 100001 110001 | 010001 011001 001001 000001
  //                              |
  //  111011 101011 100011 110011 | 010011 011011 001011 000011
  //                              |
  //  111010 101010 100010 110010 | 010010 011010 001010 000010
  // -----------------------------|-----------------------------
  //  111110 101110 100110 110110 | 010110 011110 001110 000110
  //                              |
  //  111100 101100 100100 110100 | 010100 011100 001100 000100
  //                              |
  //  111101 101101 100101 110101 | 010101 011101 001101 000101
  //                              |
  //  111111 101111 100111 110111 | 010111 011111 001111 000111
  //                              |
  // Define some symbol boundaries for unit energy constellation
  const float QAM64_BOUND1 = 2/sqrt(42.0); // Boundary between syms +1 and +3
  const float QAM64_BOUND2 = 4/sqrt(42.0); // Boundary between syms +3 and +5
  const float QAM64_BOUND3 = 6/sqrt(42.0); // Boundary between syms +5 and +7
  static const Complex8 QAM64[64] = {
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
  //    Constel::bitsToIdx = number of leftover bits that were not converted
  inline int bitsToBpskIdx(char const * const bits, const int nbits,
      int * const idx)
  {
    for(int i=0; i<nbits; ++i) idx[i] = bits[i];
    return 0;
  }
  
  inline int bitsToQpskIdx(char const * const bits, const int nbits,
      int * const idx)
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
  
  inline int bitsTo8PskIdx(char const * const bits, const int nbits,
      int * const idx)
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
  
  inline int bitsTo16QamIdx(char const * const bits, const int nbits,
      int * const idx)
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
  inline int bitsTo64QamIdx(char const * const bits, const int nbits,
      int * const idx)
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
  inline int bitsToBpsk(char const * const bits, const int nbits,
      Complex8 * const syms)
  {
    for(int i=0; i<nbits; ++i) syms[i] = BPSK[bits[i]];
    return 0;
  }
  
  inline int bitsToQpsk(char const * const bits, const int nbits,
      Complex8 * const syms)
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
  
  inline int bytesTo8Psk(char const * const bits, const int nbits,
      Complex8 * const syms)
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
  
  inline int bitsTo16Qam(char const * const bits, const int nbits,
      Complex8 * const syms)
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
  inline int bitsTo64Qam(char const * const bits, const int nbits,
      Complex8 * const syms)
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
  
  // Given a soft decision, return the index of the nearest point in the
  // specified constellation
  //    soft = soft decision to compute constellation index of
  //    Constel::indexOf = index of the nearest constellation point in the
  //        constellation buffer for the specified modulation
  inline int indexOfBPSK(const Complex8 soft)
  {
    return int(soft.real() < 0);
  }

  inline int indexOfQPSK(const Complex8 soft)
  {
    return (int(soft.real()<0)<<1 | int(soft.imag()<0));
  }

  inline int indexOf8PSK(const Complex8 soft)
  {
    return (soft.real()<0) | (soft.imag()<0)<<1 |
        (fabsf(soft.real())<fabsf(soft.imag()))<<2;
  }

  inline int indexOf16QAM(const Complex8 soft)
  {
    return (soft.real()<0)<<3 | (fabsf(soft.real())<QAM16_BOUND)<<2 |
        (soft.imag()<0)<<1 | (fabsf(soft.imag())<QAM16_BOUND);
  }

  inline int indexOf64QAM(const Complex8 soft)
  {
    return (soft.real()<0)<<5 | 


        (fabsf(soft.real())<QAM16_BOUND)<<2 |
        
        (soft.imag()<0)<<2 |
        (fabsf(soft.imag())<QAM16_BOUND);
  }

  // Given an input buffer of complex-valued soft decisions, quantizes to the
  // nearest conestallation point of the specified constellation
  //    soft = buffer of soft decisions
  //    len = number of soft decisions to quantize
  //    hard = buffer to hold output hard decisions
  inline void quantizeBPSK(Complex8 const * const soft, const int len,
      Complex8 * const hard)
  {
    for(int i=0; i<len; ++i) hard[i] = BPSK[indexOfBPSK(soft[i])];
  }

  inline void quantizeQPSK(Complex8 const * const soft, const int len,
      Complex8 * const hard)
  {
    for(int i=0; i<len; ++i) hard[i] = QPSK[indexOfQPSK(soft[i])];
  }

  inline void quantize8PSK(Complex8 const * const soft, const int len,
      Complex8 * const hard)
  {
    for(int i=0; i<len; ++i) hard[i] = PSK8[indexOf8PSK(soft[i])];
  }

  inline void quantize16QAM(Complex8 const * const soft, const int len,
      Complex8 * const hard)
  {
    for(int i=0; i<len; ++i) hard[i] = QAM16[indexOf16QAM(soft[i])];
  }

  inline void quantize64QAM(Complex8 const * const soft, const int len,
      Complex8 * const hard)
  {
    for(int i=0; i<len; ++i) hard[i] = QAM64[indexOf64QAM(soft[i])];
  }

  // Given an input buffer of soft decisions, quantizes to the nearest BPSK
  // symbol in the specified constellation and maps to bits.  Bits are output
  // as chars with 1 bit per byte, placed in the LSb of the char.
  //    soft = buffer of soft decisions to quantize
  //    len = number of soft decions to quantize
  //    bits = buffer of quantized bits, with each bit stored in the LSb of a
  //        char, must have room for len elements
  inline void quantizeBitsBPSK(Complex8 const * const soft, const int len,
      char * const bits)
  {
    for(int i=0; i<len; ++i) bits[i] = static_cast<char>(indexOfBPSK(soft[i]));
  }
  // Given an input buffer of soft decisions, quantizes to the nearest QPSK
  // symbol in the specified constellation and maps to bits.  Bits are output
  // as chars with 1 bit per byte, placed in the LSb of the char.
  //    soft = buffer of soft decisions to quantize
  //    len = number of soft decions to quantize
  //    bits = buffer of quantized bits, with each bit stored in the LSb of a
  //        char, must have room for 2*len elements
  inline void quantizeBitsQPSK(Complex8 const * const soft, const int len,
      char * bits)
  {
    char idx;
    for(int i=0; i<len; ++i)
    {
      idx = static_cast<char>(indexOfQPSK(soft[i]));
      *bits++ = idx>>1; // shift bit 2 to first position, no need to mask
      *bits++ = idx&1;  // mask off lowest bit
    }
  }
  // Given an input buffer of soft decisions, quantizes to the nearest 8-PSK
  // symbol in the specified constellation and maps to bits.  Bits are output
  // as chars with 1 bit per byte, placed in the LSb of the char.
  //    soft = buffer of soft decisions to quantize
  //    len = number of soft decions to quantize
  //    bits = buffer of quantized bits, with each bit stored in the LSb of a
  //        char, must have room for 3*len elements
  inline void quantizeBits8PSK(Complex8 const * const soft, const int len,
      char * bits)
  {
    char idx;
    for(int i=0; i<len; ++i)
    {
      idx = static_cast<char>(indexOf8PSK(soft[i]));
      *bits++ = idx>>2; // shift 3rd bit to lowest position, no need to mask
      *bits++ = (idx>>1)&1; // shift 2nd bit and mask
      *bits++ = idx&1;      // mask off lowest bit
    }
  }
  // Given an input buffer of soft decisions, quantizes to the nearest 16-QAM
  // symbol in the specified constellation and maps to bits.  Bits are output
  // as chars with 1 bit per byte, placed in the LSb of the char.
  //    soft = buffer of soft decisions to quantize
  //    len = number of soft decions to quantize
  //    bits = buffer of quantized bits, with each bit stored in the LSb of a
  //        char, must have room for 4*len elements
  inline void quantizeBits16QAM(Complex8 const * const soft, const int len,
      char * const bits)
  {
    char idx;
    for(int i=0; i<len; ++i)
    {
      idx = static_cast<char>(indexOf8PSK(soft[i]));
      *bits++ = idx>>3; // shift 4th bit to lowest position, no need to mask
      *bits++ = (idx>>2)&1; // shift 3rd bit and mask
      *bits++ = (idx>>1)&1; // shift 2nd bit and mask
      *bits++ = idx&1;      // mask off lowest bit
    }
  }

  // Given an input buffer of soft decisions, quantizes to the nearest 64-QAM
  // symbol in the specified constellation and maps to bits.  Bits are output
  // as chars with 1 bit per byte, placed in the LSb of the char.
  //    soft = buffer of soft decisions to quantize
  //    len = number of soft decions to quantize
  //    bits = buffer of quantized bits, with each bit stored in the LSb of a
  //        char, must have room for 6*len elements
  inline void quantizeBits64QAM(Complex8 const * const soft, const int len,
      char * const bits)
  {
    char idx;
    for(int i=0; i<len; ++i)
    {
      idx = static_cast<char>(indexOf8PSK(soft[i]));
      *bits++ = idx>>5; // shift 6th bit to lowest position, no need to mask
      *bits++ = (idx>>4)&1; // shift 3rd bit and mask
      *bits++ = (idx>>3)&1; // shift 2nd bit and mask
      *bits++ = (idx>>2)&1; // shift 3rd bit and mask
      *bits++ = (idx>>1)&1; // shift 2nd bit and mask
      *bits++ = idx&1;      // mask off lowest bit
    }
  }

} // namespace Constel
#endif // CONSTEL_H
