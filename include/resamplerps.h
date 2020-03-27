//==============================================================================
//  Name:     resamplerps.h
//
//  Purpose:  Performs arbitrary rate resampling and pulse-shape filtering
//            with user-defined pulse shape
//
//  Created:  2020/03/27
//
//  Description:
//    Performs arbitrary rate resampling and pulse-shaping via a polyphase
//    filter bank.  The user may specify their own filter or choose from a
//    a selection of windowed sinc or SRRC pulse shapes to create the polyphase
//    filter bank.  The phase accumulator for the resampler is double precision
//    floating point but all filter taps and filtering is performed as single
//    precision.  Although resampling rate may be arbitrary, the phase
//    precision (and thus instantaneous timing error) of a any interpolated
//    sample is determined by the number of filters in the filter bank (more
//    filters = smaller maximum phase error but filter bank takes up more
//    memory).  Each interpolation phase (sample time) is rounded to the phase
//    of the nearest filter in the bank.
//
//    Filtering (interpolation) is performed in the time domain as a
//    dot-product between data and time-reversed filter so all filters in the
//    bank are produced time-reverse.  In addition to user-defined filters, a
//    basic sinc filter with windowing function can be applied with choices of
//    the following windows:
//
//    Name            Sidelobe Equivalent   Rolloff     Window Type 
//                    Max (dB) Bandwidth  (dB/Octave) 
//    GenFilter::NONE   -13      1.00        6        Rectangle (unwindowed sinc)
//    GenFilter::HANN   -32      1.50       18        Hann (cosine bell) 
//    GenFilter::HAMM   -43      1.36        6        Hamming (bell on pedestal)
//    GenFilter::BH61   -61      1.61        6        Blackman-Harris 3 weight  
//    GenFilter::BH67   -67      1.71        6        Blackman-Harris 3 weight optimal
//    GenFilter::BH74   -74      1.79        6        Blackman-Harris 4 weight  
//    GenFilter::BH92   -92      2.00        6        Blackman-Harris 4 weight optimal
//
//    In place of the window argument, the user may specify a rolloff factor
//    for a square-root-raised-cosine (SRRC) interpolator.
//
//    This resampler can be used to maintain precision time-tagging of the
//    resampled data.  Two methods are provided for this, getDesiredSamp() and
//    getActualSamp().  Both methods return the whole and fractional index of
//    the NEXT output sample to be interpolated as an offset from the start
//    of the input data.  The difference between the two methods is that
//    getDesired uses the phase value of the internal phase accumulator (the
//    phase that the interpolator is attempting to achieve) while getActual
//    uses the phase that the interpolator was able to achieve by rounding to
//    the closest filter in the bank to the desired phase.  For most
//    applications the values returned by getDesiredSamp() will be sufficient
//    as they reflect the resampler's phase value "on average" and the
//    instantaneous jitter of the actual vs. desired will be negligeable as
//    long as there are a decent number of filters in the bank (>1000 or so).
//
//    NOTES:
//
//    1. Methods to set the phase of the resampler expect their input in
//    cycles in the range [0:1).  Methods to get the phase of the resampler
//    will return their output in cycles in the range [0:1].  The "get" method
//    can return a phase equal to 1 due to the fact that it needs to choose
//    the filter closest to the desired phase (see note 2).    
//
//    2. Since the interpoaltor is implemented as a bank of pre-computed sincs
//    rather than generating sincs on the fly (computationally intensive),
//    the instantaneous phase of the resampler is limited to a discrete set
//    of values in the range [0:1/N:1] cycles where N is the number of
//    filters in the bank.  This means that an interpolated sample can
//    potentially have a maximum phase error of +- 1/(2*N) cycles
//    vs. the expected phase of a true arbitrary rate resampler.  The setPhase 
//    and getActualPhase methods will return the actual phase achieved.
//    This may differ slightly from the internal phase accumulator which
//    tracks the desired phase and can be accessed via getDesiredPhase().
//
//    3. Based on the compiler optimizations available on your platform (SSE2)
//    there may be miniscule differences in rounding error of the phase
//    accumulator if the ratio of input to output rates is an irrational number
//    or goes to more than 15 significant digits.  This rounding error can
//    manifest as a slight rate difference (on the order of 1e-15*input rate)
//    in the output data.  The rate differences stem from the fact that,
//    depending on optimization level, the phase accumulator may be incremented
//    by 2 or 4 times the normal per-sample phase increment which is trying to
//    represent an irrational number as a double precision float.  This can
//    cause a loss of precision in the lowest digit of the phase accumulator at
//    different times in the input signal.  Note that even when resampling to
//    an irrational rate, the getPhase() methods will appropriately reflect the
//    internal accumulator so timecode can be adjusted properly.
//
//==============================================================================

#ifndef RESAMPLERPS_H
#define RESAMPLERPS_H

#include <vector>
#include <stdexcept>
#include <iomanip>
#include <math.h>
#include "complex_num.h"
#include "sincbank.h"
#include "evm.h"
#include "genfilter.h"
#include "aligned_allocator.h"

#ifdef HAVE_CONFIG_H
#include "config.h" // automatically generated, contains arch specific #defines
#endif

#ifdef HAVE_SSE2
#include <emmintrin.h> // for efficient round and trunc
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

using std::vector;
using std::logic_error;
using std::max;

//==================//
//  Resampler Class //
//==================//
class Resampler
{
  public:
    //==============//
    // Constructors //
    //==============//
    // Default Constructor, use reset() method to setup resampler
    inline Resampler();
    
    // Windowed Sinc constructor
    //    numTaps = taps per leg of the polyphase interpolation filter bank
    //    numFilters = number of filters (legs) in the polyphase interpolation
    //        filter bank, the phase precision of the resampler in cycles will
    //        be the inverse of this number
    //    bw = bandwidth of the resampling filter as a fraction of the input
    //        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
    //    win = tapering window for sinc interpolation filter, options are
    //        NONE, HANN, HAMM, BH61, BH67, BH74, BH92
    //    complex = true if input data will be complex, false if it will
    //        be real
    //    startPhase = starting phase of interpolator in cycles of the input
    //        rate, will skip or zero pad data as needed
    inline Resampler(const int numTaps, const int numFilters, const double bw,
        const GenFilter::WinType win, const bool complex=false,
        const double startPhase=0.0);
    
    // SRRC constructor
    //    numTaps = taps per leg of the polyphase interpolation filter bank
    //    numFilters = number of filters (legs) in the polyphase interpolation
    //        filter bank, the phase precision of the resampler in cycles will
    //        be the inverse of this number
    //    bw = bandwidth of the resampling filter as a fraction of the input
    //        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
    //    rolloff = rolloff factor for SRRC
    //    complex = true if input data will be complex, false if it will
    //        be real
    //    startPhase = starting phase of interpolator in cycles of the input
    //        rate, will skip or zero pad data as needed
    inline Resampler(const int numTaps, const int numFilters, const double bw,
        const double rolloff, const bool complex=false,
        const double startPhase=0.0);
    
    //============//
    // Processors //
    //============//
    // Resamples the input data to the specified rate, the output buffer
    // must be sized large enough to hold the resampled output data (use
    // the getNumOutputSamps method for assistance).  The data size and
    // is determined by the "complex" flag which was set on either reset
    // or constructor.  It can be checked with isComplex().  If the
    // complex flag is set to true, input and output data are expected to
    // be complex float data type (4 bytes real, 4 bytes imaginary, i.e.
    // Complex8).  If the complex flag is false, input and output data are
    // expected to be float data type.
    //    in = input data buffer with either Complex8 or float values
    //    inLen = number of input elements to use (samples, not bytes)
    //    phaseIncr = phase increment in cycles of the interpolator, this
    //        should be (input sampling rate)/(output sampling rate)
    //    out = output data buffer, most be sized large enough to hold all
    //        available output samples prior to this method call
    //    outLen = number of resampled output samples that were produced
    //        (samples, not bytes)
    inline void resample(const char *in, const size_t inLen,
        const double phaseIncr, char *out, size_t &outLen);
    // Wrappers for resample method that throw an exception if input/output
    // data types don't match the complex flag
    inline void resample(const Complex8 *in, const size_t inLen,
        const double phaseIncr, Complex8 *out, size_t &outLen);
    inline void resample(const float *in, const size_t inLen,
        const double phaseincr, float *out, size_t &outLen);
    // Wrappers for resample method that use STL vectors
    template <typename T, typename U>
    inline void resample(const vector<Complex8,T> &vin, const double phaseIncr,
        vector<Complex8,U> &vout);
    template <typename T, typename U>
    inline void resample(const vector<float,T> &vin, const double phaseIncr,
        vector<float,U> &vout);
    
    // Flushes the internal saved data buffer and computes final output samples
    // by zero-padding.  This may be useful when the transient filter response
    // is desired.
    //    phaseIncr = phase increment in cycles of the interpolator, this
    //        should be (input sampling rate)/(output sampling rate)
    //    out = output data buffer, most be sized large enough to hold all
    //        available output samples prior to this method call
    //    outLen = number of resampled output points that were produced
    inline void flush(const double phaseIncr, char *out, size_t &outLen);
    // Wrappers for flush method that throw an exception if input/output
    // data types don't match the complex flag
    inline void flush(const double phaseIncr, Complex8 *out, size_t &outLen);
    inline void flush(const double phaseIncr, float *out, size_t &outLen);
    // Wrappers for flush method that uses STL vectors
    template <typename T>
    inline void flush(const double phaseIncr, vector<Complex8,T> &vout);
    template <typename T>
    inline void flush(const double phaseIncr, vector<float,T> &vout);
    
    //===========//
    // Modifiers //
    //===========//
    // Updates the resampler with a new windowed sinc filter bank and resets
    // all internal buffers and parameteres (including sample count)
    //    numTaps = taps per leg of the polyphase interpolation filter bank
    //    numFilters = number of filters (legs) in the polyphase interpolation
    //        filter bank, the phase precision of the resampler in cycles will
    //        be the inverse of this number
    //    bw = bandwidth of the resampling filter as a fraction of the input
    //        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
    //    win = tapering window for sinc interpolation filter, options are
    //        NONE, HANN, HAMM, BH61, BH67, BH74, BH92
    //    complex = true if input data will be complex, false if it will
    //        be real
    //    startPhase = starting phase of interpolator in cycles of the input
    //        rate, will skip or zero pad data as needed
    inline void reset(const int numTaps, const int numFilters, const double bw,
        const GenFilter::WinType win, const bool complex=false,
        const double startPhase=0.0);
    
    // Updates the resampler with a new SRRC filter bank and resets
    // all internal buffers and parameteres (including sample count)
    //    numTaps = taps per leg of the polyphase interpolation filter bank
    //    numFilters = number of filters (legs) in the polyphase interpolation
    //        filter bank, the phase precision of the resampler in cycles will
    //        be the inverse of this number
    //    bw = bandwidth of the resampling filter as a fraction of the input
    //        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
    //    rolloff = rolloff factor for SRRC
    //    complex = true if input data will be complex, false if it will
    //        be real
    //    startPhase = starting phase of interpolator in cycles of the input
    //        rate, will skip or zero pad data as needed
    inline void reset(const int numTaps, const int numFilters, const double bw,
        const double rolloff, const bool complex=false,
        const double startPhase=0.0);
    
    // Clears the internal save buffer and resets the resampler phase to
    // startPhase without reconstructing the filter bank
    //    startPhase = desired interpolation phase in cycles
    inline void clear(const double startPhase=0.0);
    
    // Sets the number of OpenMP threads and the threads blocktime (currently
    // hardcoded depending on number of threads).  If the HAVE_OPENMP macro
    // is not defined, threads will be set to 1.
    //    threads = number of OpenMP threads to allow the resampler to use
    inline void setThreads(const int threads);
    
    //============//
    // Inspectors //
    //============//
    // Given the current resampler parameters and the specified phase increment
    // per sample, returns an estimate of the number of output samples that
    // will be produced by the next call to resample() if len input samples are
    // provided.  Due to floating point rounding errors with certain compiler
    // optimizations, this should be taken as an estimate and my differ by a
    // sample from the number actually produced.
    //    phaseIncr = phase increment in cycles of the interpolator, this
    //        should be (input sampling rate)/(output sampling rate)
    //    len = number of input samples that will be provided to next call of
    //        resample()
    //    Resampler::getNumOutputSamps = number of resampled points that will
    //        be produced by next call to resample()
    inline int getNumOutputSamps(const double phaseIncr, const int len) const;
    // Wrapper for getNumOutputSamps which returns the number of output
    // samples from a flush() call rather than resample()
    inline int getNumOutputSamps(const double phaseIncr);
    
    // Returns the whole and fractional sample index of the next output sample
    // in cycles referenced from the first input sample since the last call
    // to reset() or clear().  This method returns the values that the
    // interpolator is attempting to achieve which may differ from the values
    // that are actually achieved.  It is guaranteed that frac will be in the
    // range [0:1].
    //    whole = whole sample index of next interpolated output point,
    //        referenced from first input sample
    //    frac = fractional sample index (phase in cycles) of next interpolated
    //        output point
    inline void getDesiredSamp(long long &whole, double &frac) const;
    
    // Returns the whole and fractional sample index of the next output sample
    // in cycles referenced from the first input sample since the last call
    // to reset() or clear().  This method returns the values that the
    // interpolator is actually able to achieve given the fact that there are
    // a discrete number of interpolation filters in the bank. It is
    // guaranteed that frac will be in the range [0:1].
    //    whole = whole sample index of next interpolated output point,
    //        referenced from first input sample
    //    frac = fractional sample index (phase in cycles) of next interpolated
    //        output point
    inline void getActualSamp(long long &whole, double &frac) const;
    
    // Returns the current phase accumulator value in cycles, should be
    // identical to the frac value returned by getDesiredPhase, this method
    // just provided for compatibility with older apps
    inline double getPhase() const {return _phaseAccum;};
    
    // Returns whether the resampler is currently configured for real or
    // complex input data
    //    Resampler::isComplex = true if resampler is configured for Complex8
    //        input data, false if it is configured for float input data
    inline bool isComplex() const {return _isCx;};
    
  protected:
    double _phaseAccum; // phase accumulator in cycles
    double _phaseAdjust;// phase accumulator adjustment when converting
                        // interpolation phase to index of filter in bank
    int _savedElems;    // number of elements saved from last input
    int _numFilters;    // number of filters (legs) in filter bank
    int _numTaps;       // number of taps per leg of filter bank
    int _initPad;       // initial zero-padding to produce transient response
    double _bw;         // bandwidth of interpolation filter as a fraction of
                        // input sampling rate
    bool _isCx;         // whether or not filter is currently in complex mode
    long long _sampIdx; // current input sample that _phaseaccum corresponds to
    double _rolloff;    // rolloff factor for SRRC
    int _elemSize;      // input element size in bytes depending on whether
                        // processing mode is set to real or complex
    int _dotpLen;       // length of the filtering dot-product in samples
    int _threads;       // number of OpenMP threads class is allowed to use
    SincBank _filtBank; // bank of sinc or SRRC interpolation filters
    GenFilter::WinType _win;    // tapering window for sinc filter
    vector<char> _vsave;// saved data from previous processing call
    vector<int> _vsaveIdx;  // index in saved data buffer to start
                            // dot-product on for each output sample
    vector<int> _vfiltIdx;  // filter index of each dot-product for
                            // current output buffer
    
    // Reset internal variables and clear the saved data buffer, but does
    // not modify the filter bank.
    //    numTaps = taps per leg of the polyphase interpolation filter bank
    //    numFilters = number of filters (legs) in the polyphase interpolation
    //        filter bank, the phase precision of the resampler in cycles will
    //        be the inverse of this number
    //    bw = bandwidth of the resampling filter as a fraction of the input
    //        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
    //    complex = true if input data will be complex float, false if it will
    //        be scalar float
    //    startPhase = starting phase of interpolator in cycles of the input
    //        rate, will skip or zero pad data as needed
    inline void _reset(const int numTaps, const int numFilters, const double bw,
        const bool complex, const double startPhase);
    
    // Returns the index of the interpolation filter that will be used to
    // interpolate the input argument based on the rounding scheme and
    // compiler optimization level
    //    phase = phase in cycles to interpolate
    //    Resampler::_getFiltIdx = index of interpolation filter that will
    //        be used to interpolate phase
    inline int _getFiltIdx(const double phase) const;
    
    // Given the phase increment per sample, returns an estimate of the number
    // of output samples, computes the starting sample in vsave and the filter
    // index in the filter bank of each dot-product and stores the values in
    // _vsaveIdx and _vfiltIdx
    //    phaseIncr = phase increment in cycles of the interpolator, this
    //        should be (input sampling rate)/(output sampling rate)
    //    Resampler::_computeInterpPhases = number of output samples that can
    //        be computed
    inline int _computeInterpPhases(const double phaseIncr);
};

//==============//
// Constructors //
//==============//
// Default Constructor, use reset() method to setup resampler
inline Resampler::Resampler()
{
  setThreads(1);
  reset(3,10,1,GenFilter::HANN,false,0.0);
}

// Windowed Sinc constructor
//    numTaps = taps per leg of the polyphase interpolation filter bank
//    numFilters = number of filters (legs) in the polyphase interpolation
//        filter bank, the phase precision of the resampler in cycles will
//        be the inverse of this number
//    bw = bandwidth of the resampling filter as a fraction of the input
//        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
//    win = tapering window for sinc interpolation filter, options are
//        NONE, HANN, HAMM, BH61, BH67, BH74, BH92
//    complex = true if input data will be complex, false if it will
//        be real
//    startPhase = starting phase of interpolator in cycles of the input
//        rate, will skip or zero pad data as needed
inline Resampler::Resampler(const int numTaps, const int numFilters,
    const double bw, const GenFilter::WinType win, const bool complex,
    const double startPhase)
{
  setThreads(1);
  reset(numTaps, numFilters, bw, win, complex, startPhase);
}

// SRRC constructor
//    numTaps = taps per leg of the polyphase interpolation filter bank
//    numFilters = number of filters (legs) in the polyphase interpolation
//        filter bank, the phase precision of the resampler in cycles will
//        be the inverse of this number
//    bw = bandwidth of the resampling filter as a fraction of the input
//        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
//    rolloff = rolloff factor for SRRC
//    complex = true if input data will be complex, false if it will
//        be real
//    startPhase = starting phase of interpolator in cycles of the input
//        rate, will skip or zero pad data as needed
inline Resampler::Resampler(const int numTaps, const int numFilters,
    const double bw, const double rolloff, const bool complex,
    const double startPhase)
{
  setThreads(1);
  reset(numTaps, numFilters, bw, rolloff, complex, startPhase);
}

//============//
// Processors //
//============//
// Resamples the input data to the specified rate, the output buffer
// must be sized large enough to hold the resampled output data (use
// the getNumOutputSamps method for assistance).  The data size and
// is determined by the "complex" flag which was set on either reset
// or constructor.  It can be checked with isComplex().  If the
// complex flag is set to true, input and output data are expected to
// be complex float data type (4 bytes real, 4 bytes imaginary, i.e.
// Complex8).  If the complex flag is false, input and output data are
// expected to be float data type.  Individual wrappers are also provided
// which throw an exception if the input and output data types do not
// match the current complex flag.
//    in = input data buffer with either Complex8 or float values
//    inLen = number of input elements to use (samples, not bytes)
//    phaseIncr = phase increment in cycles of the interpolator, this
//        should be (input sampling rate)/(output sampling rate)
//    out = output data buffer, most be sized large enough to hold all
//        available output samples prior to this method call
//    outLen = number of resampled output samples that were produced
//        (samples, not bytes)
inline void Resampler::resample(const char *in, const size_t inLen,
    const double phaseIncr, char *out, size_t &outLen)
{
  int tmp = _savedElems+inLen; // new size of saved data buffer
  // copy new samples to vsave, taking into account the fact that _savedElems
  // might be negative indicating that we need to skip over some elements
  int pts2copy = tmp - max(_savedElems,0);
  if(pts2copy > 0)
  {
    if(int(_vsave.size()) < tmp*_elemSize+EVM_ALIGNMENT)
    {
      _vsave.resize(tmp*_elemSize+EVM_ALIGNMENT);
    }
    memcpy(&_vsave[max(_savedElems,0)*_elemSize],
        &in[(inLen-pts2copy)*_elemSize], pts2copy*_elemSize);
  }
  // compute number of output samples and filter indices for dot-product
  _savedElems = tmp;
  tmp = _computeInterpPhases(phaseIncr); // store number of output samps as int
  if(_isCx)
  {
    #ifdef HAVE_OPENMP
    #pragma omp parallel for
    #endif
    // resample the saved data buffer, storing interpolated samples in out
    for(int i=0; i<tmp; ++i)
    {
      EVM::int_dotpsUA(
          reinterpret_cast<Complex8*>(&_vsave[_vsaveIdx[i]*_elemSize]),
          _filtBank.getFiltPointer(_vfiltIdx[i]), _dotpLen,
          reinterpret_cast<Complex8&>(out[i*_elemSize]));
    }
  }
  else
  {
    #ifdef HAVE_OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i<tmp; ++i)
    {
      EVM::int_dotpUA(
          reinterpret_cast<float*>(&_vsave[_vsaveIdx[i]*_elemSize]),
          _filtBank.getFiltPointer(_vfiltIdx[i]), _dotpLen,
          reinterpret_cast<float&>(out[i*_elemSize]));
    }
  }
  // move any remaining elements to the beginning of the save buffer for next
  // pass
  outLen = size_t(tmp);
  if(outLen > 0)
  {
    // determine the next saved data index that is required and how far it is
    // from the end of _vsave
    tmp = int(trunc(_phaseAccum+_phaseAdjust));
    _phaseAccum -= tmp;
    tmp += _vsaveIdx[outLen-1];
    pts2copy = _savedElems-tmp;
    if(pts2copy > 0)
    {
      // use overlap-safe memmove here just in case
      memmove(&_vsave[0],&_vsave[tmp*_elemSize],pts2copy*_elemSize);
    }
    _savedElems = pts2copy;
    _sampIdx += tmp;
  }
  return;
}

// Wrappers for resample method that throw an exception if input/output
// data types don't match the complex flag
inline void Resampler::resample(const Complex8 *in, const size_t inLen,
    const double phaseIncr, Complex8 *out, size_t &outLen)
{
  if(!_isCx)
  {
    throw logic_error("[Resampler::resample] Complex input data provided "
        "but current configuration requires real (float) data.");
  }
  resample(reinterpret_cast<const char*>(in), inLen, phaseIncr,
      reinterpret_cast<char*>(out), outLen);
  return;
}

inline void Resampler::resample(const float *in, const size_t inLen,
    const double phaseIncr, float *out, size_t &outLen)
{
  if(_isCx)
  {
    throw logic_error("[Resampler::resample] Real (float) input data provided "
        "but current configuration requires complex data.");
  }
  resample(reinterpret_cast<const char*>(in), inLen, phaseIncr,
      reinterpret_cast<char*>(out), outLen);
  return;
}

// Wrappers for resample method that uses STL vectors
template <typename T, typename U>
inline void Resampler::resample(const vector<Complex8,T> &vin,
    const double phaseIncr, vector<Complex8,U> &vout)
{
  size_t maxSize = getNumOutputSamps(phaseIncr, vin.size())+1;
  vout.resize(maxSize);
  resample(&vin[0],vin.size(),phaseIncr,&vout[0],maxSize);
  if(maxSize != vout.size()) vout.resize(maxSize);
  return;
}
template <typename T, typename U>
inline void Resampler::resample(const vector<float,T> &vin,
    const double phaseIncr, vector<float,U> &vout)
{
  size_t maxSize = getNumOutputSamps(phaseIncr, vin.size())+1;
  vout.resize(maxSize);
  resample(&vin[0],vin.size(),phaseIncr,&vout[0],maxSize);
  if(maxSize != vout.size()) vout.resize(maxSize);
  return;
}

// Flushes the internal saved data buffer and computes final output samples
// by zero-padding.  This may be useful when the transient filter response
// is desired.
//    phaseIncr = phase increment in cycles of the interpolator, this
//        should be (input sampling rate)/(output sampling rate)
//    out = output data buffer, most be sized large enough to hold all
//        available output samples prior to this method call
//    outLen = number of resampled output points that were produced
inline void Resampler::flush(const double phaseIncr, char *out, size_t &outLen)
{
  // create a final input buffer of 0-padding to account for 1/2 filter length
  // then just resample as normal
  size_t zpad = _numTaps-1-_initPad;
  vector<char> in(zpad*_elemSize,0);
  resample(&in[0],zpad,phaseIncr,out,outLen);
  return;
}

// Wrappers for flush method that throw an exception if input/output
// data types don't match the complex flag
inline void Resampler::flush(const double phaseIncr, Complex8 *out,
    size_t &outLen)
{
  if(!_isCx)
  {
    throw logic_error("[Resampler::flush] Output buffer format is complex "
        "but current resampler configuration is for real (float) data.");
  }
  flush(phaseIncr, reinterpret_cast<char*>(out), outLen);
  return;
}
inline void Resampler::flush(const double phaseIncr, float *out, size_t &outLen)
{
  if(_isCx)
  {
    throw logic_error("[Resampler::flush] Output buffer format is real (float) "
        "but current resampler configuration is for complex data.");
  }
  flush(phaseIncr, reinterpret_cast<char*>(out), outLen);
  return;
}

// Wrappers for flush method that uses STL vectors
template <typename T>
inline void Resampler::flush(const double phaseIncr, vector<Complex8,T> &vout)
{
  size_t maxSize = getNumOutputSamps(phaseIncr)+1;
  vout.resize(maxSize);
  flush(phaseIncr,&vout[0],maxSize);
  if(maxSize != vout.size()) vout.resize(maxSize);
  return;
}
template <typename T>
inline void Resampler::flush(const double phaseIncr, vector<float,T> &vout)
{
  size_t maxSize = getNumOutputSamps(phaseIncr)+1;
  vout.resize(maxSize);
  flush(phaseIncr,&vout[0],maxSize);
  if(maxSize != vout.size()) vout.resize(maxSize);
  return;
}

//===========//
// Modifiers //
//===========//
// Updates the resampler with a new windowed sinc filter bank and resets
// all internal buffers and parameteres (including sample count)
//    numTaps = taps per leg of the polyphase interpolation filter bank
//    numFilters = number of filters (legs) in the polyphase interpolation
//        filter bank, the phase precision of the resampler in cycles will
//        be the inverse of this number
//    bw = bandwidth of the resampling filter as a fraction of the input
//        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
//    win = tapering window for sinc interpolation filter, options are
//        NONE, HANN, HAMM, BH61, BH67, BH74, BH92
//    complex = true if input data will be complex, false if it will
//        be real
//    startPhase = starting phase of interpolator in cycles of the input
//        rate, will skip or zero pad data as needed
inline void Resampler::reset(const int numTaps, const int numFilters,
    const double bw, const GenFilter::WinType win, const bool complex,
    const double startPhase)
{
  // reset buffers and common parameters
  _reset(numTaps,numFilters,bw,complex,startPhase);
  // construct the windowed sinc interpolation filter bank
  // for complex data, repeat each tap (mode=2) so we can use the more
  // efficient complex-real dot-product
  if(_isCx)
  {
    _filtBank.update(_numFilters,_numTaps,_bw,win,2,true,true);
  }
  else
  {
    _filtBank.update(_numFilters,_numTaps,_bw,win,0,true,true);
  }
  return;
}

// Updates the resampler with a new SRRC filter bank and resets
// all internal buffers and parameteres (including sample count)
//    numTaps = taps per leg of the polyphase interpolation filter bank
//    numFilters = number of filters (legs) in the polyphase interpolation
//        filter bank, the phase precision of the resampler in cycles will
//        be the inverse of this number
//    bw = bandwidth of the resampling filter as a fraction of the input
//        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
//    rolloff = rolloff factor for SRRC
//    complex = true if input data will be complex, false if it will
//        be real
//    startPhase = starting phase of interpolator in cycles of the input
//        rate, will skip or zero pad data as needed
inline void Resampler::reset(const int numTaps, const int numFilters,
    const double bw, const double rolloff, const bool complex,
    const double startPhase)
{
  // reset buffers and common parameters
  _reset(numTaps,numFilters,bw,complex,startPhase);
  // construct the windowed sinc interpolation filter bank
  // for complex data, repeat each tap (mode=2) so we can use the more
  // efficient complex-real dot-product
  if(_isCx)
  {
    _filtBank.update(_numFilters,_numTaps,_bw,rolloff,2,true,true);
  }
  else
  {
    _filtBank.update(_numFilters,_numTaps,_bw,rolloff,0,true,true);
  }
  return;
}

// Clears the internal save buffer and resets the resampler phase to
// startPhase without reconstructing the filter bank.  If startPhase is
// <0 the data will be zero-padded to accomodate it.  If startPhase is
// >1, input data will be consumed to accomodate it.
//    startPhase = desired initial interpolation phase in cycles
inline void Resampler::clear(const double startPhase)
{
  _sampIdx = 0;
  // convert the phase argument to be in the range [0:1), if it was negative
  // need to add more zero-padding to account for samples we don't have, if
  // it was positive may need to skip forward (discard) some samples
  int zpad = _initPad-int(floor(startPhase));
  _sampIdx += int(floor(startPhase));
  _phaseAccum = startPhase-floor(startPhase); // now should be in [0:1) range
  
  // if the bank has an odd number of taps per filter, the interpolation phases
  // will actually go rom [-.5 : .5] rather than [0 : 1]
  if(_phaseAdjust > 0 && _phaseAccum > .5)
  {
    // convert phase to [-.5 : 0) range and move forward by a sample
    _phaseAccum -= 1;
    ++_sampIdx;
    --zpad;
  }
  // reset the saved data buffer with the transient response padding
  _savedElems = zpad;
  _vsave.resize(max(zpad,0)*_elemSize);
  memset(&_vsave[0],0,_vsave.size()); //_vsave is char so size is in bytes
  _vsaveIdx.resize(0);
  _vfiltIdx.resize(0);
  return;
}

// Sets the number of OpenMP threads and the threads blocktime (currently
// hardcoded depending on number of threads).  If the HAVE_OPENMP macro
// is not defined, threads will be set to 1.
//    threads = number of OpenMP threads to allow the resampler to use
inline void Resampler::setThreads(const int threads)
{
  #ifdef HAVE_OPENMP
  _threads = max(threads,1);
  //if(_threads > 1) kmp_set_blocktime(10);
  //else kmp_set_blocktime(0);  // runs faster on 1 thread if OMP isn't blocking
  omp_set_num_threads(_threads);
  #else
  _threads = 1;
  #endif
  return;
}

//============//
// Inspectors //
//============//
// Given the current resampler parameters and the specified phase increment
// per sample, returns an estimate of the number of output samples that
// will be produced by the next call to resample() if len input samples are
// provided.  Due to floating point rounding errors with certain compiler
// optimizations, this should be taken as an estimate and my differ by a
// sample from the number actually produced.
//    phaseIncr = phase increment in cycles of the interpolator, this
//        should be (input sampling rate)/(output sampling rate)
//    len = number of input samples that will be provided to next call of
//        resample()
//    Resampler::getNumOutputSamps = number of resampled points that will
//        be produced by next call to resample()
inline int Resampler::getNumOutputSamps(const double phaseIncr, const int len)
    const
{
  // compute last sample we can interpolate from the saved data buffer once
  // the new elements are added
  double lastInterp = len+_savedElems-_numTaps+1-_phaseAdjust;
  // compute how many times we can increment the current phase accumulator
  // before we exceed the last start sample, note that the number of output
  // samples will be the number of increments we can do plus 1
  return max(int(floor((lastInterp-_phaseAccum)/phaseIncr))+1,0);
}

// Wrapper for getNumOutputSamps which returns the number of output
// samples from a flush() call rather than resample()
inline int Resampler::getNumOutputSamps(const double phaseIncr)
{
  return getNumOutputSamps(phaseIncr,_numTaps-1-_initPad);
}

// Returns the whole and fractional sample index of the next output sample
// in cycles referenced from the first input sample since the last call
// to reset() or clear().  This method returns the values that the
// interpolator is attempting to achieve which may differ from the values
// that are actually achieved.  It is guaranteed that frac will be in the
// range [0:1].
//    whole = whole sample index of next interpolated output point,
//        referenced from first input sample
//    frac = fractional sample index (phase in cycles) of next interpolated
//        output point
inline void Resampler::getDesiredSamp(long long &whole, double &frac) const
{
  whole = _sampIdx;
  frac = _phaseAccum;
  if(frac < 0)
  {
    frac += 1;
    whole -= 1;
  }
  return;
}

// Returns the whole and fractional sample index of the next output sample
// in cycles referenced from the first input sample since the last call
// to reset() or clear().  This method returns the values that the
// interpolator is actually able to achieve given the fact that there are
// a discrete number of interpolation filters in the bank. It is
// guaranteed that frac will be in the range [0:1].
//    whole = whole sample index of next interpolated output point,
//        referenced from first input sample
//    frac = fractional sample index (phase in cycles) of next interpolated
//        output point
inline void Resampler::getActualSamp(long long &whole, double &frac) const
{
  whole = _sampIdx;
  // determine the actual phase achieved by the resampler by computing the
  // index of the filter it will choose and converting this index back to
  // the actual interpolation phase
  frac = _getFiltIdx(_phaseAccum)/double(_numFilters)-_phaseAdjust;
  if(frac < 0)
  {
    frac += 1;
    whole -= 1;
  }
  return;
}

//===================//
// Protected Methods //
//===================//
// Reset internal variables and clear the saved data buffer
//    numTaps = taps per leg of the polyphase interpolation filter bank
//    numFilters = number of filters (legs) in the polyphase interpolation
//        filter bank, the phase precision of the resampler in cycles will
//        be the inverse of this number
//    bw = bandwidth of the resampling filter as a fraction of the input
//        signal sampling rate (i.e. bw=.95 is 95% of the input rate)
//    complex = true if input data will be complex float, false if it will
//        be scalar float
//    startPhase = starting phase of interpolator in cycles of the input
//        rate, will skip or zero pad data as needed
inline void Resampler::_reset(const int numTaps, const int numFilters,
    const double bw, const bool complex, const double startPhase)
{
  // set imternal params
  _numTaps = numTaps;
  _numFilters = numFilters;
  _bw = bw;
  _isCx = complex;
  // do some bounds checking
  if(_numTaps <= 0)
  {
    throw logic_error("[Resampler::reset] Invalid number of taps specified, "
          "taps must be a positive non-zero value");
  }
  if(_numFilters <= 0)
  {
    throw logic_error("[Resampler::reset] Invalid number of filters specified, "
          "filters must be a positive non-zero value");
  }
  if(_bw <= 0)
  {
    throw logic_error("[Resampler::reset] Invalid resampling filter bandwidth "
          "specified, bandwidth must be a positive non-zero value");
  }
  // set dependent parameters
  if(_isCx) _elemSize = sizeof(Complex8);
  else _elemSize = sizeof(float);
  // for filter banks with odd number of taps per leg, the phase=0 filter will
  // be in the center of the bank and the index=0 filter will be phase=-.5,
  // for this reason we need to apply an adjustment factor of .5 when doing
  // the phase to index conversion for banks with odd number of taps
  if(_numTaps%2==0)
  {
    _phaseAdjust = 0;
    _initPad = _numTaps/2-1;
  }
  else
  {
    _phaseAdjust = .5;
    _initPad = (_numTaps-1)/2;
  }
  // ensure dotproduct length is a multiple of the required byte alignment
  // for EVM internal methods
  if(complex) _dotpLen = EVM::getSize<Complex8>(_numTaps);
  else _dotpLen = EVM::getSize<float>(_numTaps);
  clear(startPhase);
  return;
}

// Returns the index of the interpolation filter that will be used to
// interpolate the input argument based on the rounding scheme and
// compiler optimization level
//    phase = phase in cycles to interpolate
//    Resampler::_getFiltIdx = index of interpolation filter that will
//        be used to interpolate phase
inline int Resampler::_getFiltIdx(const double phase) const
{
  #ifdef HAVE_SSE2
  // even though SSE2 rounds .5 toward the nearest even value, this is
  // ok here since rounding up or down will not affect the alignment of
  // our filter with the data, just bear in mind that when there are an
  // odd number of filters in the bank, this may result in very slight
  // differences in output when compiled with SSE2 vs. without SSE2 due
  // to choice of a different filter when _phaseAccum+_phaseAdjust is exactly .5
  return int(_mm_cvtsd_si32(_mm_set_sd((phase+_phaseAdjust)*
      _numFilters)));
  #else
  return int(round((phase+_phaseAdjust)*_numFilters));
  #endif
}

// Given the phase increment per sample, returns an estimate of the number
// of output samples, computes the starting sample in vsave and the filter
// index in the filter bank of each dot-product and stores the values in
// _vsaveIdx and _vfiltIdx
//    phaseIncr = phase increment in cycles of the interpolator, this
//        should be (input sampling rate)/(output sampling rate)
//    Resampler::_computeInterpPhases = number of output samples that can
//        be computed
inline int Resampler::_computeInterpPhases(const double phaseIncr)
{
  // need to ensure that all buffers are sized appropriately to avoid
  // over-indexing when doing SSE store operations, for this reason compute
  // how many extra integers are needed to produce a buffer size that is a
  // multiple of the alignment
  int extraInts = int(ceil(float(EVM_ALIGNMENT)/sizeof(int)));
  int nout = getNumOutputSamps(phaseIncr,0);
  // resize the save index and filter index buffers if they are too small
  // to accomodate all the output samples.  Add an extra sample of buffer
  // room since we may be using SSE2 instructions to compute the filter index
  // and thus processing 2 samples at a time
  if(_vsaveIdx.size() < size_t(nout+extraInts))
  {
    _vsaveIdx.resize(nout+extraInts);
    _vfiltIdx.resize(_vsaveIdx.size());
  }
  int curSave = 0; // starting index in _vsave for current _phaseAccum 
  #ifdef HAVE_SSE2
  __m128d incr = _mm_set1_pd(2*phaseIncr);
  __m128d adj = _mm_set1_pd(_phaseAdjust);
  __m128d nfilt = _mm_set1_pd(double(_numFilters));
  __m128i sidx = _mm_set1_epi32(curSave); // initial saved data index
  __m128d ph = _mm_setr_pd(_phaseAccum,_phaseAccum+phaseIncr);// starting phase
  __m128d trncd;  // double precision truncated phase
  __m128i trnci;  // integer truncated phase
  for(int i=0; i<nout; i+=2)
  {
    // add the phase adjustment offset and truncate the phase value to determine
    // how many whole samples we must skip forward
    trncd = _mm_add_pd(ph,adj);      // add offset (.5 or 0)
    trnci = _mm_cvttpd_epi32(trncd); // truncate
    sidx = _mm_add_epi32(sidx,trnci);// update save buffer indices for phase ph
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&_vsaveIdx[i]),sidx);
    // convert the truncated phase values (whole index updates) to dp float
    // and subtract from the phase + adjustment to get a value between [0:1)
    // for computing the filter index
    trncd = _mm_cvtepi32_pd(trnci);
    ph = _mm_sub_pd(ph,trncd);
    trncd = _mm_add_pd(ph,adj);
    // convert the value between [0:1) to the filter index
    trncd = _mm_mul_pd(trncd,nfilt); // multiply by number of filters
    trnci = _mm_cvtpd_epi32(trncd);  // round to nearest filter index
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&_vfiltIdx[i]),trnci);
    // increment the phase for the next loop
    ph = _mm_add_pd(ph,incr);
  }
  // Since the resample method assumes the phase accumulator will be one
  // increment past the last sample and will be referenced from the last
  // element in _vsaveIdx, need to adjust the final increment to ph to be
  // be a single increment instead of a double and need to make sure we
  // do this to the specific phase that was applied to the final output
  // sample
  incr = _mm_set1_pd(phaseIncr);
  ph = _mm_sub_pd(ph, incr);
  if(nout%2==0) _mm_storeh_pd(&_phaseAccum,ph);
  else _mm_storel_pd(&_phaseAccum,ph);
  #else
  int sampIncr;
  for(int i=0; i<nout; ++i)
  {
    // check if we need to skip forward by a sample to accomodate this phase
    sampIncr = int(trunc(_phaseAccum+_phaseAdjust));
    _phaseAccum -= sampIncr;
    curSave += sampIncr;
    // save current parameters and increment phase accumulator
    _vsaveIdx[i] = curSave;
    _vfiltIdx[i] = _getFiltIdx(_phaseAccum);
    _phaseAccum += phaseIncr;
  }
  #endif
  return nout;
}

#endif // RESAMPLERPS_H
