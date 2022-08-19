/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*- vi:set ts=8 sts=4 sw=4: */
/*
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "resampler.h"
// #include "../base///Profiler.h"

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <algorithm>

#include "../system/memallocators.h"

#define USE_SPEEX // temp fixed def

#ifdef HAVE_IPP
#include <ippversion.h>
#if (IPP_VERSION_MAJOR < 7)
#include <ipps.h>
#include <ippsr.h>
#include <ippac.h>
#else
#include <ipps.h>
#endif
#endif

#ifdef HAVE_LIBSAMPLERATE
#include <samplerate.h>
#endif

#ifdef HAVE_LIBRESAMPLE
#include <libresample.h>
#endif

#ifdef USE_SPEEX
#include "../speex/speex_resampler.h"
#endif

#ifndef HAVE_IPP
#ifndef HAVE_LIBSAMPLERATE
#ifndef HAVE_LIBRESAMPLE
#ifndef USE_SPEEX
#error No resampler implementation selected!
#endif
#endif
#endif
#endif

namespace audiomod {

class resamplerImpl
{
public:
    virtual ~resamplerImpl() { }
    
    virtual int doresample(const float *const R__ *const R__ in, 
                         float *const R__ *const R__ out,
                         int incount,
                         float ratio,
                         bool final) = 0;
    
    // virtual int resampleInterleaved(const float *const R__ in, 
    //                                 float *const R__ out,
    //                                 int incount,
    //                                 float ratio,
    //                                 bool final) = 0;

    virtual int numchannels() const = 0;

    virtual void reset() = 0;
};

namespace resamplers {

#ifdef HAVE_IPP

class RS_IPP : public resamplerImpl
{
public:
    RS_IPP(resampler::RSMode quality, int chans, int maxBufferSize);
    ~RS_IPP();

    int doresample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final);

    int resampleInterleaved(const float *const R__ in,
                            float *const R__ out,
                            int incount,
                            float ratio,
                            bool final = false);

    int numchannels() const { return chans_; }

    void reset();

protected:
    // to m_outbuf
    int doResample(int outcount, double ratio, bool final);

    IppsResamplingPolyphase_32f **m_state;
    float **m_inbuf;
    size_t m_inbufsz;
    float **m_outbuf;
    size_t m_outbufsz;
    int bufsize_;
    int chans_;
    int m_window;
    float m_factor;
    int m_history;
    int *m_lastread;
    double *m_time;
    
    void setBufSize(int);
};

RS_IPP::RS_IPP(resampler::RSMode quality, int chans, int maxBufferSize) :
    m_state(0),
    chans_(chans)
{

    int nStep = 16;
    IppHintAlgorithm hint = ippAlgHintFast;

    switch (quality) {

    case resampler::Best:
        m_window = 64;
        nStep = 80;
        hint = ippAlgHintAccurate;
        break;

    case resampler::FastestTolerable:
        nStep = 16;
        m_window = 16;
        hint = ippAlgHintFast;
        break;

    case resampler::Fastest:
        m_window = 24;
        nStep = 64;
        hint = ippAlgHintFast;
        break;
    }

    m_factor = 8; // initial upper bound on m_ratio, may be amended later
    m_history = int(m_window * 0.5 * std::max(1.0, 1.0 / m_factor)) + 1;

    m_state = new IppsResamplingPolyphase_32f *[chans_];

    m_lastread = new int[chans_];
    m_time = new double[chans_];

    m_inbufsz = 0;
    m_outbufsz = 0;
    m_inbuf = 0;
    m_outbuf = 0;
    bufsize_ = 0;

    setBufSize(maxBufferSize + m_history);

#if (IPP_VERSION_MAJOR >= 7)
    int specSize = 0;
    ippsResamplePolyphaseGetSize_32f(float(m_window),
                                     nStep,
                                     &specSize,
                                     hint);
    if (specSize == 0) {
#ifndef NO_EXCEPTIONS
        throw resampler::ImplementationError;
#else        
        abort();
#endif
    }
#endif

    for (int c = 0; c < chans_; ++c) {
#if (IPP_VERSION_MAJOR < 7)
        ippsResamplePolyphaseInitAlloc_32f(&m_state[c],
                                           float(m_window),
                                           nStep,
                                           0.95f,
                                           9.0f,
                                           hint);
#else
        m_state[c] = (IppsResamplingPolyphase_32f *)ippsMalloc_8u(specSize);
        ippsResamplePolyphaseInit_32f(float(m_window),
                                      nStep,
                                      0.95f,
                                      9.0f,
                                      m_state[c],
                                      hint);
#endif
        
        m_lastread[c] = m_history;
        m_time[c] = m_history;
    }

}

RS_IPP::~RS_IPP()
{
#if (IPP_VERSION_MAJOR < 7)
    for (int c = 0; c < chans_; ++c) {
        ippsResamplePolyphaseFree_32f(m_state[c]);
    }
#else
    for (int c = 0; c < chans_; ++c) {
        ippsFree(m_state[c]);
    }
#endif

    mem_deallocate_channels(m_inbuf, chans_);
    mem_deallocate_channels(m_outbuf, chans_);

    delete[] m_lastread;
    delete[] m_time;
    delete[] m_state;
}

void
RS_IPP::setBufSize(int sz)
{

    bufsize_ = sz;

    int n1 = bufsize_ + m_history + 2;
    int n2 = lrintf(ceil((bufsize_ - m_history) * m_factor + 2));

    m_inbuf = mem_reallocate_and_zeroext_channels
        (m_inbuf, chans_, m_inbufsz, chans_, n1);

    m_outbuf = mem_reallocate_and_zeroext_channels
        (m_outbuf, chans_, m_outbufsz, chans_, n2);
            
    m_inbufsz = n1;
    m_outbufsz = n2;
}

int
RS_IPP::doresample(const float *const R__ *const R__ in,
                float *const R__ *const R__ out,
                int incount,
                float ratio,
                bool final)
{
    if (ratio > m_factor) {
        m_factor = ratio;
        m_history = int(m_window * 0.5 * std::max(1.0, 1.0 / m_factor)) + 1;
    }

    for (int c = 0; c < chans_; ++c) {
        if (m_lastread[c] + incount + m_history > bufsize_) {
            setBufSize(m_lastread[c] + incount + m_history);
        }
    }

    for (int c = 0; c < chans_; ++c) {
        for (int i = 0; i < incount; ++i) {
            m_inbuf[c][m_lastread[c] + i] = in[c][i];
        }
        m_lastread[c] += incount;
    }

    int got = doResample(int(round(incount * ratio)), ratio, final);
    
    for (int c = 0; c < chans_; ++c) {
        vector_copy(out[c], m_outbuf[c], got);
    }

    return got;
}

int
RS_IPP::doResample(int outspace, double ratio, bool final)
{
    int outcount = 0;
    
    for (int c = 0; c < chans_; ++c) {

        int n = m_lastread[c] - m_history - int(m_time[c]);

        // We're committed to not overrunning outspace, so we need to
        // offer the resampler only enough samples to ensure it won't

        int limit = int(floor(outspace / ratio));
        if (n > limit) {

            n = limit;
        }
        
#if (IPP_VERSION_MAJOR < 7)
        ippsResamplePolyphase_32f(m_state[c],
                                  m_inbuf[c],
                                  n,
                                  m_outbuf[c],
                                  ratio,
                                  1.0f,
                                  &m_time[c],
                                  &outcount);
#else
        ippsResamplePolyphase_32f(m_inbuf[c],
                                  n,
                                  m_outbuf[c],
                                  ratio,
                                  1.0f,
                                  &m_time[c],
                                  &outcount,
                                  m_state[c]);
#endif

        int t = int(round(m_time[c]));

        vector_move(m_inbuf[c],
               m_inbuf[c] + t - m_history,
               m_lastread[c] + m_history - t);

        m_lastread[c] -= t - m_history;
        m_time[c] -= t - m_history;

        
        if (final && n < limit) {

            // Looks like this actually produces too many samples
            // (additionalcount is a few samples too large).

            // Also, we aren't likely to have enough space in the
            // output buffer as the caller won't have allowed for
            // all the samples we're retrieving here.

            // What to do?

            int additionalcount = 0;
            
            for (int i = 0; i < m_history; ++i) {
                m_inbuf[c][m_lastread[c] + i] = 0.f;
            }

            int nAdditional = m_lastread[c] - int(m_time[c]);

            if (n + nAdditional > limit) {
                nAdditional = limit - n;
            }
            
#if (IPP_VERSION_MAJOR < 7)
            ippsResamplePolyphase_32f(m_state[c],
                                      m_inbuf[c],
                                      nAdditional,
                                      m_outbuf[c],
                                      ratio,
                                      1.0f,
                                      &m_time[c],
                                      &additionalcount);
#else
            ippsResamplePolyphase_32f(m_inbuf[c],
                                      nAdditional,
                                      m_outbuf[c],
                                      ratio,
                                      1.0f,
                                      &m_time[c],
                                      &additionalcount,
                                      m_state[c]);
#endif

            if (c == 0) {
                outcount += additionalcount;
            }
        }
    }
    
    return outcount;
}

void
RS_IPP::reset()
{
    //!!!
}

#endif /* HAVE_IPP */

#ifdef HAVE_LIBSAMPLERATE

class RS_SRC : public resamplerImpl
{
public:
    RS_SRC(resampler::RSMode quality, int chans, int maxBufferSize);
    ~RS_SRC();

    int doresample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final);


    int numchannels() const { return chans_; }

    void reset();

protected:
    SRC_STATE *m_src;
    float *iin_;
    float *iout_;
    float m_lastRatio;
    int chans_;
    int iin_size;
    int iout_size;
};

RS_SRC::RS_SRC(resampler::RSMode quality, int chans, int maxBufferSize) :
    m_src(0),
    iin_(0),
    iout_(0),
    m_lastRatio(1.f),
    chans_(chans),
    iin_size(0),
    iout_size(0)
{

    int err = 0;
    m_src = src_new(quality == resampler::Best ? SRC_SINC_BEST_QUALITY :
                    quality == resampler::Fastest ? SRC_LINEAR :
                    SRC_SINC_FASTEST,
                    chans, &err);

    if (err) {
        std::cerr << "resampler::resampler: failed to create libsamplerate resampler: " 
                  << src_strerror(err) << std::endl;
#ifndef NO_EXCEPTIONS
        throw resampler::ImplementationError;
#endif
    }

    if (maxBufferSize > 0 && chans_ > 1) {
        iin_size = maxBufferSize * chans_;
        iout_size = maxBufferSize * chans_ * 2;
        iin_ = mem_allocate<float>(iin_size);
        iout_ = mem_allocate<float>(iout_size);
    }

    reset();
}

RS_SRC::~RS_SRC()
{
    src_delete(m_src);
    mem_deallocate(iin_);
    mem_deallocate(iout_);
}

int
RS_SRC::doresample(const float *const R__ *const R__ in,
                float *const R__ *const R__ out,
                int incount,
                float ratio,
                bool final)
{
    SRC_DATA data;

    int outcount = lrintf(ceilf(incount * ratio));

    if (chans_ == 1) {
        data.data_in = const_cast<float *>(*in); //!!!???
        data.data_out = *out;
    } else {
        if (incount * chans_ > iin_size) {
            iin_ = mem_reallocate<float>(iin_, iin_size, incount * chans_);
            iin_size = incount * chans_;
        }
        if (outcount * chans_ > iout_size) {
            iout_ = mem_reallocate<float>(iout_, iout_size, outcount * chans_);
            iout_size = outcount * chans_;
        }
        vector_interleave(iin_, in, chans_, incount);
        data.data_in = iin_;
        data.data_out = iout_;
    }

    data.input_frames = incount;
    data.output_frames = outcount;
    data.src_ratio = ratio;
    data.end_of_input = (final ? 1 : 0);

    int err = src_process(m_src, &data);

    if (err) {
        std::cerr << "resampler::process: libsamplerate error: "
                  << src_strerror(err) << std::endl;
#ifndef NO_EXCEPTIONS
        throw resampler::ImplementationError;
#endif
    }

    if (chans_ > 1) {
        vector_deinterleave(out, iout_, chans_, data.output_frames_gen);
    }

    m_lastRatio = ratio;

    return data.output_frames_gen;
}

void
RS_SRC::reset()
{
    src_reset(m_src);
}

#endif /* HAVE_LIBSAMPLERATE */

#ifdef HAVE_LIBRESAMPLE

class RS_Resample : public resamplerImpl
{
public:
    RS_Resample(resampler::RSMode quality, int chans, int maxBufferSize);
    ~RS_Resample();

    int doresample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final);

    int numchannels() const { return chans_; }

    void reset();

protected:
    void *m_src;
    float *iin_;
    float *iout_;
    float m_lastRatio;
    int chans_;
    int iin_size;
    int iout_size;
};

RS_Resample::RS_Resample(resampler::RSMode quality, int chans, int maxBufferSize) :
    m_src(0),
    iin_(0),
    iout_(0),
    m_lastRatio(1.f),
    chans_(chans),
    iin_size(0),
    iout_size(0)
{

    float min_factor = 0.125f;
    float max_factor = 8.0f;

    m_src = resample_open(quality == resampler::Best ? 1 : 0, min_factor, max_factor);

    if (!m_src) {
        std::cerr << "resampler::resampler: failed to create libresample resampler: " 
                  << std::endl;
        throw resampler::ImplementationError; //!!! of course, need to catch this!
    }

    if (maxBufferSize > 0 && chans_ > 1) {
        iin_size = maxBufferSize * chans_;
        iout_size = maxBufferSize * chans_ * 2;
        iin_ = mem_allocate<float>(iin_size);
        iout_ = mem_allocate<float>(iout_size);
    }

    reset();
}

RS_Resample::~RS_Resample()
{
    resample_close(m_src);
    if (iin_size > 0) {
        mem_deallocate(iin_);
    }
    if (iout_size > 0) {
        mem_deallocate(iout_);
    }
}

int
RS_Resample::doresample(const float *const R__ *const R__ in,
                     float *const R__ *const R__ out,
                     int incount,
                     float ratio,
                     bool final)
{
    float *data_in;
    float *data_out;
    int input_frames, output_frames, end_of_input, source_used;
    float src_ratio;

    int outcount = lrintf(ceilf(incount * ratio));

    if (chans_ == 1) {
        data_in = const_cast<float *>(*in); //!!!???
        data_out = *out;
    } else {
        if (incount * chans_ > iin_size) {
            iin_ = mem_reallocate<float>(iin_, iin_size, incount * chans_);
            iin_size = incount * chans_;
        }
        if (outcount * chans_ > iout_size) {
            iout_ = mem_reallocate<float>(iout_, iout_size, outcount * chans_);
            iout_size = outcount * chans_;
        }
        vector_interleave(iin_, in, chans_, incount);
        data_in = iin_;
        data_out = iout_;
    }

    input_frames = incount;
    output_frames = outcount;
    src_ratio = ratio;
    end_of_input = (final ? 1 : 0);

    int output_frames_gen = resample_process(m_src,
                                             src_ratio,
                                             data_in,
                                             input_frames,
                                             end_of_input,
                                             &source_used,
                                             data_out,
                                             output_frames);

    if (output_frames_gen < 0) {
        std::cerr << "resampler::process: libresample error: "
                  << std::endl;
        throw resampler::ImplementationError; //!!! of course, need to catch this!
    }

    if (chans_ > 1) {
        vector_deinterleave(out, iout_, chans_, output_frames_gen);
    }

    m_lastRatio = ratio;

    return output_frames_gen;
}

void
RS_Resample::reset()
{
}

#endif /* HAVE_LIBRESAMPLE */

#ifdef USE_SPEEX
    
class RS_Speex : public resamplerImpl
{
public:
    RS_Speex(resampler::RSMode quality, int chans, int maxBufferSize);
    ~RS_Speex();

    int doresample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final);

    int numchannels() const { return chans_; }

    void reset();

protected:
    SpeexResamplerState *resampler_;
    float *iin_; // internal inbuffer, the size could be changed with every resample operation
    float *iout_; // internal outbuffer, the size could be changed with every resample operation
    int chans_;
    unsigned int iin_size;
    unsigned int iout_size;
    float lastratio_;
    bool initial_;

    void setratio(float);
};

RS_Speex::RS_Speex(resampler::RSMode quality, int chans, int maxBufferSize) :
    resampler_(0),
    iin_(0),
    iout_(0),
    chans_(chans),
    iin_size(0),
    iout_size(0),
    lastratio_(1),
    initial_(true)
{
    int q = (quality == resampler::Best ? 10 :
             quality == resampler::Fastest ? 0 : 4);

    int err = 0;
    resampler_ = speex_resampler_init_frac(chans_,
                                            1, 1,
                                            48000, 48000, // irrelevant
                                            q,
                                            &err);
    

    if (err) {
        std::cerr << "resampler::resampler: failed to create Speex resampler" 
                  << std::endl;
#ifndef NO_EXCEPTIONS
        throw resampler::ImplementationError;
#endif
    }

    if (maxBufferSize > 0 && chans_ > 1) {
        iin_size = maxBufferSize * chans_;
        iout_size = maxBufferSize * chans_ * 2;
        iin_ = mem_allocate<float>(iin_size);
        iout_ = mem_allocate<float>(iout_size);
    }
}

RS_Speex::~RS_Speex()
{
    speex_resampler_destroy(resampler_);
    mem_deallocate<float>(iin_);
    mem_deallocate<float>(iout_);
}

void
RS_Speex::setratio(float ratio)
{
    // Speex wants a ratio of two unsigned integers, not a single
    // float.  Let's do that.

    unsigned int big = 272408136U; 
    unsigned int denom = 1, num = 1;

    if (ratio < 1.f) {
        denom = big;
        double dnum = double(big) * double(ratio);
        num = (unsigned int)dnum;
    } else if (ratio > 1.f) {
        num = big;
        double ddenom = double(big) / double(ratio);
        denom = (unsigned int)ddenom;
    }
    
    speex_resampler_set_rate_frac
        (resampler_, denom, num, 48000, 48000);
    
    speex_resampler_get_ratio(resampler_, &denom, &num);
    
    lastratio_ = ratio;

    if (initial_) {
        speex_resampler_skip_zeros(resampler_);
        initial_ = false;
    }
}

int
RS_Speex::doresample(const float *const R__ *const R__ in,
                  float *const R__ *const R__ out,
                  int incount,
                  float ratio,
                  bool final)
{
    if (ratio != lastratio_) {
        setratio(ratio);
    }

    unsigned int uincount = incount;
    unsigned int outcount = lrintf(ceilf(incount * ratio)); //!!! inexact now

    float *data_in, *data_out;

    if (chans_ == 1) {
        data_in = const_cast<float *>(*in);
        data_out = *out;
    } else {
        if (uincount * chans_ > iin_size) {
            iin_ = mem_reallocate<float>(iin_, iin_size, uincount * chans_);
            iin_size = uincount * chans_;
        }
        if (outcount * chans_ > iout_size) {
            iout_ = mem_reallocate<float>(iout_, iout_size, outcount * chans_);
            iout_size = outcount * chans_;
        }
        vector_interleave(iin_, in, chans_, incount);
        data_in = iin_;
        data_out = iout_;
    }

    speex_resampler_process_interleaved_float(resampler_,
                                                        data_in,
                                                        &uincount,
                                                        data_out,
                                                        &outcount);


    if (chans_ > 1) {
        vector_deinterleave(out, iout_, chans_, outcount);
    }

    return outcount;
}

void
RS_Speex::reset()
{
    lastratio_ = -1.0; // force reset of ratio
    initial_ = true;
    speex_resampler_reset_mem(resampler_);
}

#endif

} /* end namespace resamplers */

resampler::resampler(resampler::RSMode quality, int chans, int maxBufferSize)
{
    method_ = -1;
    
    switch (quality) {

    case resampler::Best:
#ifdef HAVE_IPP
        method_ = 0;
#endif
#ifdef USE_SPEEX
        method_ = 2;
#endif
#ifdef HAVE_LIBRESAMPLE
        method_ = 3;
#endif
#ifdef HAVE_LIBSAMPLERATE
        method_ = 1;
#endif
        break;

    case resampler::FastestTolerable:
#ifdef HAVE_IPP
        method_ = 0;
#endif
#ifdef HAVE_LIBRESAMPLE
        method_ = 3;
#endif
#ifdef HAVE_LIBSAMPLERATE
        method_ = 1;
#endif
#ifdef USE_SPEEX
        method_ = 2;
#endif
        break;

    case resampler::Fastest:
#ifdef HAVE_IPP
        method_ = 0;
#endif
#ifdef HAVE_LIBRESAMPLE
        method_ = 3;
#endif
#ifdef USE_SPEEX
        method_ = 2;
#endif
#ifdef HAVE_LIBSAMPLERATE
        method_ = 1;
#endif
        break;
    }

    if (method_ == -1) {
        std::cerr << "resampler::resampler(" << quality << ", " << chans
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
    }

    switch (method_) {
    case 0:
#ifdef HAVE_IPP
        rs = new resamplers::RS_IPP(quality, chans, maxBufferSize);
#else
        std::cerr << "resampler::resampler(" << quality << ", " << chans
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
#endif
        break;

    case 1:
#ifdef HAVE_LIBSAMPLERATE
        rs = new resamplers::RS_SRC(quality, chans, maxBufferSize);
#else
        std::cerr << "resampler::resampler(" << quality << ", " << chans
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
#endif
        break;

    case 2:
#ifdef USE_SPEEX
        rs = new resamplers::RS_Speex(quality, chans, maxBufferSize);
#else
        std::cerr << "resampler::resampler(" << quality << ", " << chans
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
#endif
        break;

    case 3:
#ifdef HAVE_LIBRESAMPLE
        rs = new resamplers::RS_Resample(quality, chans, maxBufferSize);
#else
        std::cerr << "resampler::resampler(" << quality << ", " << chans
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
#endif
        break;
    }

    if (!rs) {
        std::cerr << "resampler::resampler(" << quality << ", " << chans
                  << ", " << maxBufferSize
                  << "): Internal error: No implementation selected"
                  << std::endl;
        abort();
    }
}

resampler::~resampler()
{
    delete rs;
}

int 
resampler::doresample(const float *const R__ *const R__ in,
                    float *const R__ *const R__ out,
                    int incount, float ratio, bool final)
{
    //Profiler profiler("resampler::doresample");
    return rs->doresample(in, out, incount, ratio, final);
}

int
resampler::numchannels() const
{
    return rs->numchannels();
}

void
resampler::reset()
{
    rs->reset();
}

}
