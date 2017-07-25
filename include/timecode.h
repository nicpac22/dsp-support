/*
 * vim: sw=2
 * Copyright (C) 2014 Nicholas Xenias
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
//==============================================================================
//
//  Program:  timecode.h
//
//  Purpose:  Basic class for doing arithmetic with time since an epoch
//  
//  Author:   NSX
//
//  Created:  2014/04/12
//
//  Description:
//    Very basic class for keeping track of time since an epoch and tagging
//    it to a sample index.  Provides basic arithmatic functions for
//    manipulating the timecode as well.
//    
//==============================================================================

#ifndef TIMECODE_H_
#define TIMECODE_H_

#include <math.h>
#include <ostream>
#include <stdint.h>

#define TIMECODE_H_TOLR 1e-10 // for external use to compare timecodes

class TimeCode
{
  public:
    double whole;   // whole seconds since epoch
    double frac;    // fractional second offset from whole
    int64_t samp;   // sample index associated with this timecode
    int status;     // 1 = good timecode, 0 = bad timecode
  
    //==============//
    // Constructors //
    //==============//
    inline TimeCode()
    {
      whole=0.0; frac=0.0; samp=0; status=0; return;
    }
    
    inline TimeCode(const double wh, const double fr, const int64_t idx,
        const int st=1)
    {
      whole=wh; frac=fr; samp=idx; status=st; return;
    }
    
    inline TimeCode(const double all, const int64_t idx, const int st=1)
    {
      whole=trunc(all); frac=all-trunc(all); samp=idx; status=st; return;
    }
    
    inline TimeCode(const TimeCode &tc)
    {
      whole=tc.whole; frac=tc.frac; samp=tc.samp; status=tc.status; return;
    }
    
    //===========//
    // Modifiers //
    //===========//
    // Corrects timecode to ensure the fractional part is positive and less
    // than 1, does not change the value of the timecode just its representation
    inline void correct()
    {
      if(fabs(this->frac) >= 1.0) // correct for too large a value
      {
        double newwh = trunc(this->frac);
        this->frac -= newwh;
        this->whole += newwh;
      }
      if(this->frac < 0) // correct for negative value
      {
        this->whole -= 1;
        this->frac = 1+this->frac;
      }
      return;
    }
    
    // Adjusts the timecode to reflect a specific sample index, given the
    // sampling interval.  This method will uset the timecode's current samp
    // field and adjust it to reflect the "newsamp"
    //    newsamp = whole sample index to adjust the timecode to reflect
    //    interval = sampling interval (time between samples) in seconds
    inline void adjust(const int64_t newsamp, const double interval)
    {
      // volatile might help maintain precision in following steps...
      volatile int64_t diff = newsamp - samp;
      volatile double rhsFr = diff*interval;
      volatile double rhsWh = trunc(rhsFr);
      rhsFr -= rhsWh;
      whole = whole+rhsWh+trunc(frac+rhsFr);
      frac = frac+rhsFr-trunc(frac+rhsFr);
      samp = newsamp;
      return;
    }
    
    //============//
    // Inspectors //
    //============//
    // Returns the timecode as a single double with whole and fractional parts
    // summed together, this can result in a major loss or precision
    inline double r8() const {return this->whole+this->frac;}
    
    // Returns the difference between this timecode and another in seconds
    inline double diff(const TimeCode &rhs) const
    {
      volatile double whDiff = this->whole-rhs.whole;
      return whDiff+(this->frac-rhs.frac);
    }
    
    //=================//
    // Unary Operators //
    //=================//
    // Negate both the whole and fractional seconds
    TimeCode operator- () const
    {
      return TimeCode(-this->whole,-this->frac,this->samp,this->status);
    }
    
    //==================//
    // Binary Operators //
    //==================//
    // Assignment
    TimeCode& operator= (const TimeCode &rhs)
    {
      this->whole=rhs.whole;
      this->frac=rhs.frac;
      this->samp=rhs.samp;
      this->status=rhs.status;
      return *this;
    }
    // Addition, TimeCode
    TimeCode operator+ (const TimeCode &rhs) const
    {
      double tmp = this->frac+rhs.frac;
      return TimeCode(this->whole+rhs.whole+trunc(tmp), tmp-trunc(tmp),
          this->samp, this->status);
    }
    // Addition, whole and fractional seconds combined in one double
    TimeCode operator+ (const double &rhs) const
    {
      double rhsWh = trunc(rhs);
      // volatile might help maintain precision in the next step
      volatile double rhsFr = rhs-rhsWh;
      double tmp = this->frac+rhsFr;
      return TimeCode(this->whole+rhsWh+trunc(tmp), tmp-trunc(tmp),
          this->samp, this->status);
    }
    // Addition assignment, TimeCode
    TimeCode& operator+= (const TimeCode &rhs)
    {
      double tmp = this->frac+rhs.frac;
      this->whole += (trunc(tmp)+rhs.whole);
      this->frac = tmp-trunc(tmp);
      return *this;
    }
    // Addition assignment, whole and fractional seconds combined in one double
    TimeCode& operator+= (const double &rhs)
    {
      double rhsWh = trunc(rhs);
      // volatile might help maintain precision in the next step...
      volatile double rhsFr = rhs-rhsWh;
      double tmp = this->frac+rhsFr;
      this->whole += (trunc(tmp)+rhsWh);
      this->frac = tmp-trunc(tmp);
      return *this;
    }
    // Subtraction, TimeCode
    TimeCode operator- (const TimeCode &rhs) const
    {
      double tmp = this->frac-rhs.frac;
      return TimeCode(this->whole-rhs.whole+trunc(tmp),tmp-trunc(tmp),
          this->samp,this->status);
    }
    // Subtraction, whole and fractional seconds combined in one double
    TimeCode operator- (const double &rhs) const
    {
      double rhsWh = trunc(rhs);
      // volatile might help maintain precision in the next step...
      volatile double rhsFr = rhs-rhsWh;
      double tmp = this->frac-rhsFr;
      return TimeCode(this->whole-rhsWh+trunc(tmp),tmp-trunc(tmp),this->samp,
          this->status);
    }
    // Subtraction assignment, TimeCode
    TimeCode& operator-= (const TimeCode &rhs)
    {
      double tmp = this->frac-rhs.frac;
      this->whole -= (rhs.whole-trunc(tmp));
      this->frac = tmp-trunc(tmp);
      return *this;
    }
    // Subtraction assignment, whole and fractional seconds combined
    TimeCode& operator-= (const double &rhs)
    {
      double rhsWh = trunc(rhs);
      // volatile might help maintain precision in the next step...
      volatile double rhsFr = rhs-rhsWh;
      double tmp = this->frac-rhsFr;
      this->whole -= (rhsWh-trunc(tmp));
      this->frac = tmp-trunc(tmp);
      return *this;
    }
    
    //===================//
    // Boolean Operators //
    //===================//
    // Greater than
    bool operator> (const TimeCode &rhs) const
    {
      return ((this->whole-rhs.whole)+(this->frac-rhs.frac)>0);
    }
    // Less than
    bool operator< (const TimeCode &rhs) const
    {
      return ((this->whole-rhs.whole)+(this->frac-rhs.frac)<0);
    }
};

// Output operator for the TimeCode class, note that this function does not
// need to be a friend of the class since the data members of TimeCode are
// public
inline std::ostream& operator<< (std::ostream& output, const TimeCode &tc)
{
  output<<"("<<tc.whole<<", "<<tc.frac<<", "<<tc.samp<<", "<<tc.status<<")";
  return output;
}

#endif // TIMECODE_H_
