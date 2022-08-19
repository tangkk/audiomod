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

#include "FFT.h"
#include "../system/memallocators.h"
#include "../system/vectorization.h"


#define USE_KISSFFT // temp fixed def

#ifdef HAVE_IPP
#include <ippversion.h>
#include <ipps.h>
#endif

#ifdef HAVE_FFTW3
#include <fftw3.h>
#endif

#ifdef HAVE_VDSP
#include <Accelerate/Accelerate.h>
#endif

#ifdef HAVE_MEDIALIB
#include <mlib_signal.h>
#endif

#ifdef HAVE_OPENMAX
#include <omxSP.h>
#endif

#ifdef HAVE_SFFT
extern "C" {
#include <sfft.h>
}
#endif

#ifdef USE_KISSFFT
#include "../kissfft/kiss_fftr.h"
//#include "../../../ali_denoise/kiss_fftr.h" // FIXME: try move kissfft to a separate folder
#endif

#ifndef HAVE_IPP
#ifndef HAVE_FFTW3
#ifndef USE_KISSFFT
#ifndef USE_BUILTIN_FFT
#ifndef HAVE_VDSP
#ifndef HAVE_MEDIALIB
#ifndef HAVE_OPENMAX
#ifndef HAVE_SFFT
#error No FFT implementation selected!
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#include <cmath>
#include <iostream>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef FFT_MEASUREMENT
#ifndef _WIN32
#include <unistd.h>
#endif
#endif

namespace audiomod {

class FFTImpl
{
public:
    virtual ~FFTImpl() { }

    virtual FFT::Precisions getSupportedPrecisions() const = 0;

    virtual void initFloat() = 0;
    virtual void initDouble() = 0;

    virtual void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) = 0;
    virtual void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) = 0;
    virtual void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) = 0;
    virtual void forwardMagnitude(const double *R__ realIn, double *R__ magOut) = 0;

    virtual void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) = 0;
    virtual void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) = 0;
    virtual void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) = 0;
    virtual void forwardMagnitude(const float *R__ realIn, float *R__ magOut) = 0;

    virtual void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) = 0;
    virtual void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) = 0;
    virtual void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) = 0;
    virtual void inverseCepstral(const double *R__ magIn, double *R__ cepOut) = 0;

    virtual void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) = 0;
    virtual void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) = 0;
    virtual void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) = 0;
    virtual void inverseCepstral(const float *R__ magIn, float *R__ cepOut) = 0;
};    

namespace FFTs {

#ifdef HAVE_IPP

class D_IPP : public FFTImpl
{
public:
    D_IPP(int size) :
        size_(size), m_fspec(0), m_dspec(0)
    { 
        for (int i = 0; ; ++i) {
            if (size_ & (1 << i)) {
                m_order = i;
                break;
            }
        }
    }

    ~D_IPP() {
        if (m_fspec) {
#if (IPP_VERSION_MAJOR >= 9)
            ippsFree(m_fspecbuf);
#else
            ippsFFTFree_R_32f(m_fspec);
#endif
            ippsFree(fbuf_);
            ippsFree(fpacked_);
            ippsFree(m_fspare);
        }
        if (m_dspec) {
#if (IPP_VERSION_MAJOR >= 9)
            ippsFree(m_dspecbuf);
#else
            ippsFFTFree_R_64f(m_dspec);
#endif
            ippsFree(m_dbuf);
            ippsFree(m_dpacked);
            ippsFree(m_dspare);
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision | FFT::DoublePrecision;
    }

    //!!! rv check

    void initFloat() {
        if (m_fspec) return;
#if (IPP_VERSION_MAJOR >= 9)
        int specSize, specBufferSize, bufferSize;
        ippsFFTGetSize_R_32f(m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                             &specSize, &specBufferSize, &bufferSize);
        m_fspecbuf = ippsMalloc_8u(specSize);
        Ipp8u *tmp = ippsMalloc_8u(specBufferSize);
        fbuf_ = ippsMalloc_8u(bufferSize);
        fpacked_ = ippsMalloc_32f(size_ + 2);
        m_fspare = ippsMalloc_32f(size_ / 2 + 1);
        ippsFFTInit_R_32f(&m_fspec,
                          m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                          m_fspecbuf, tmp);
        ippsFree(tmp);
#else
        int specSize, specBufferSize, bufferSize;
        ippsFFTGetSize_R_32f(m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                             &specSize, &specBufferSize, &bufferSize);
        fbuf_ = ippsMalloc_8u(bufferSize);
        fpacked_ = ippsMalloc_32f(size_ + 2);
        m_fspare = ippsMalloc_32f(size_ / 2 + 1);
        ippsFFTInitAlloc_R_32f(&m_fspec, m_order, IPP_FFT_NODIV_BY_ANY, 
                               ippAlgHintFast);
#endif
    }

    void initDouble() {
        if (m_dspec) return;
#if (IPP_VERSION_MAJOR >= 9)
        int specSize, specBufferSize, bufferSize;
        ippsFFTGetSize_R_64f(m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                             &specSize, &specBufferSize, &bufferSize);
        m_dspecbuf = ippsMalloc_8u(specSize);
        Ipp8u *tmp = ippsMalloc_8u(specBufferSize);
        m_dbuf = ippsMalloc_8u(bufferSize);
        m_dpacked = ippsMalloc_64f(size_ + 2);
        m_dspare = ippsMalloc_64f(size_ / 2 + 1);
        ippsFFTInit_R_64f(&m_dspec,
                          m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                          m_dspecbuf, tmp);
        ippsFree(tmp);
#else
        int specSize, specBufferSize, bufferSize;
        ippsFFTGetSize_R_64f(m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                             &specSize, &specBufferSize, &bufferSize);
        m_dbuf = ippsMalloc_8u(bufferSize);
        m_dpacked = ippsMalloc_64f(size_ + 2);
        m_dspare = ippsMalloc_64f(size_ / 2 + 1);
        ippsFFTInitAlloc_R_64f(&m_dspec, m_order, IPP_FFT_NODIV_BY_ANY, 
                               ippAlgHintFast);
#endif
    }

    void packFloat(const float *R__ real, const float *R__ imag) {
        //Profiler profiler("D_IPP::packFloat");
        int index = 0;
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            fpacked_[index++] = real[i];
            index++;
        }
        index = 0;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                fpacked_[index++] = imag[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                fpacked_[index++] = 0.f;
            }
        }
    }

    void packDouble(const double *R__ real, const double *R__ imag) {
        //Profiler profiler("D_IPP::packDouble");
        int index = 0;
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            m_dpacked[index++] = real[i];
            index++;
        }
        index = 0;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_dpacked[index++] = imag[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_dpacked[index++] = 0.0;
            }
        }
    }

    void unpackFloat(float *real, float *R__ imag) { // real may be equal to fpacked_
        //Profiler profiler("D_IPP::unpackFloat");
        int index = 0;
        const int hs = size_/2;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                imag[i] = fpacked_[index++];
            }
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            real[i] = fpacked_[index++];
            index++;
        }
    }        

    void unpackDouble(double *real, double *R__ imag) { // real may be equal to m_dpacked
        //Profiler profiler("D_IPP::unpackDouble");
        int index = 0;
        const int hs = size_/2;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                imag[i] = m_dpacked[index++];
            }
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            real[i] = m_dpacked[index++];
            index++;
        }
    }        

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        //Profiler profiler("D_IPP::forward [d]");
        if (!m_dspec) initDouble();
        ippsFFTFwd_RToCCS_64f(realIn, m_dpacked, m_dspec, m_dbuf);
        unpackDouble(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        //Profiler profiler("D_IPP::forwardInterleaved [d]");
        if (!m_dspec) initDouble();
        ippsFFTFwd_RToCCS_64f(realIn, complexOut, m_dspec, m_dbuf);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        //Profiler profiler("D_IPP::forwardPolar [d]");
        if (!m_dspec) initDouble();
        ippsFFTFwd_RToCCS_64f(realIn, m_dpacked, m_dspec, m_dbuf);
        unpackDouble(m_dpacked, m_dspare);
        //Profiler profiler2("D_IPP::forwardPolar [d] conv");
        ippsCartToPolar_64f(m_dpacked, m_dspare, magOut, phaseOut, size_/2+1);
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        //Profiler profiler("D_IPP::forwardMagnitude [d]");
        if (!m_dspec) initDouble();
        ippsFFTFwd_RToCCS_64f(realIn, m_dpacked, m_dspec, m_dbuf);
        unpackDouble(m_dpacked, m_dspare);
        ippsMagnitude_64f(m_dpacked, m_dspare, magOut, size_/2+1);
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        //Profiler profiler("D_IPP::forward [f]");
        if (!m_fspec) initFloat();
        ippsFFTFwd_RToCCS_32f(realIn, fpacked_, m_fspec, fbuf_);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        //Profiler profiler("D_IPP::forwardInterleaved [f]");
        if (!m_fspec) initFloat();
        ippsFFTFwd_RToCCS_32f(realIn, complexOut, m_fspec, fbuf_);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        //Profiler profiler("D_IPP::forwardPolar [f]");
        if (!m_fspec) initFloat();
        ippsFFTFwd_RToCCS_32f(realIn, fpacked_, m_fspec, fbuf_);
        unpackFloat(fpacked_, m_fspare);
        //Profiler profiler2("D_IPP::forwardPolar [f] conv");
        ippsCartToPolar_32f(fpacked_, m_fspare, magOut, phaseOut, size_/2+1);
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        //Profiler profiler("D_IPP::forwardMagnitude [f]");
        if (!m_fspec) initFloat();
        ippsFFTFwd_RToCCS_32f(realIn, fpacked_, m_fspec, fbuf_);
        unpackFloat(fpacked_, m_fspare);
        ippsMagnitude_32f(fpacked_, m_fspare, magOut, size_/2+1);
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        //Profiler profiler("D_IPP::inverse [d]");
        if (!m_dspec) initDouble();
        packDouble(realIn, imagIn);
        ippsFFTInv_CCSToR_64f(m_dpacked, realOut, m_dspec, m_dbuf);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        //Profiler profiler("D_IPP::inverse [d]");
        if (!m_dspec) initDouble();
        ippsFFTInv_CCSToR_64f(complexIn, realOut, m_dspec, m_dbuf);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        //Profiler profiler("D_IPP::inversePolar [d]");
        if (!m_dspec) initDouble();
        ippsPolarToCart_64f(magIn, phaseIn, realOut, m_dspare, size_/2+1);
        //Profiler profiler2("D_IPP::inversePolar [d] postconv");
        packDouble(realOut, m_dspare); // to m_dpacked
        ippsFFTInv_CCSToR_64f(m_dpacked, realOut, m_dspec, m_dbuf);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        //Profiler profiler("D_IPP::inverseCepstral [d]");
        if (!m_dspec) initDouble();
        const int hs1 = size_/2 + 1;
        ippsCopy_64f(magIn, m_dspare, hs1);
        ippsAddC_64f_I(0.000001, m_dspare, hs1);
        ippsLn_64f_I(m_dspare, hs1);
        packDouble(m_dspare, 0);
        ippsFFTInv_CCSToR_64f(m_dpacked, cepOut, m_dspec, m_dbuf);
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        //Profiler profiler("D_IPP::inverse [f]");
        if (!m_fspec) initFloat();
        packFloat(realIn, imagIn);
        ippsFFTInv_CCSToR_32f(fpacked_, realOut, m_fspec, fbuf_);
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        //Profiler profiler("D_IPP::inverse [f]");
        if (!m_fspec) initFloat();
        ippsFFTInv_CCSToR_32f(complexIn, realOut, m_fspec, fbuf_);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        //Profiler profiler("D_IPP::inversePolar [f]");
        if (!m_fspec) initFloat();
        ippsPolarToCart_32f(magIn, phaseIn, realOut, m_fspare, size_/2+1);
        //Profiler profiler2("D_IPP::inversePolar [f] postconv");
        packFloat(realOut, m_fspare); // to fpacked_
        ippsFFTInv_CCSToR_32f(fpacked_, realOut, m_fspec, fbuf_);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        //Profiler profiler("D_IPP::inverseCepstral [f]");
        if (!m_fspec) initFloat();
        const int hs1 = size_/2 + 1;
        ippsCopy_32f(magIn, m_fspare, hs1);
        ippsAddC_32f_I(0.000001f, m_fspare, hs1);
        ippsLn_32f_I(m_fspare, hs1);
        packFloat(m_fspare, 0);
        ippsFFTInv_CCSToR_32f(fpacked_, cepOut, m_fspec, fbuf_);
    }

private:
    const int size_;
    int m_order;
    IppsFFTSpec_R_32f *m_fspec;
    IppsFFTSpec_R_64f *m_dspec;
    Ipp8u *m_fspecbuf;
    Ipp8u *m_dspecbuf;
    Ipp8u *fbuf_;
    Ipp8u *m_dbuf;
    float *fpacked_;
    float *m_fspare;
    double *m_dpacked;
    double *m_dspare;
};

#endif /* HAVE_IPP */

#ifdef HAVE_VDSP

class D_VDSP : public FFTImpl
{
public:
    D_VDSP(int size) :
        size_(size), m_fspec(0), m_dspec(0),
        fpacked_(0), m_fspare(0),
        m_dpacked(0), m_dspare(0)
    { 
        for (int i = 0; ; ++i) {
            if (size_ & (1 << i)) {
                m_order = i;
                break;
            }
        }
    }

    ~D_VDSP() {
        if (m_fspec) {
            vDSP_destroy_fftsetup(m_fspec);
            mem_deallocate(m_fspare);
            mem_deallocate(m_fspare2);
            mem_deallocate(fbuf_->realp);
            mem_deallocate(fbuf_->imagp);
            delete fbuf_;
            mem_deallocate(fpacked_->realp);
            mem_deallocate(fpacked_->imagp);
            delete fpacked_;
        }
        if (m_dspec) {
            vDSP_destroy_fftsetupD(m_dspec);
            mem_deallocate(m_dspare);
            mem_deallocate(m_dspare2);
            mem_deallocate(m_dbuf->realp);
            mem_deallocate(m_dbuf->imagp);
            delete m_dbuf;
            mem_deallocate(m_dpacked->realp);
            mem_deallocate(m_dpacked->imagp);
            delete m_dpacked;
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision | FFT::DoublePrecision;
    }

    //!!! rv check

    void initFloat() {
        if (m_fspec) return;
        m_fspec = vDSP_create_fftsetup(m_order, FFT_RADIX2);
        fbuf_ = new DSPSplitComplex;
        //!!! "If possible, tempBuffer->realp and tempBuffer->imagp should be 32-byte aligned for best performance."
        fbuf_->realp = mem_allocate<float>(size_);
        fbuf_->imagp = mem_allocate<float>(size_);
        fpacked_ = new DSPSplitComplex;
        fpacked_->realp = mem_allocate<float>(size_ / 2 + 1);
        fpacked_->imagp = mem_allocate<float>(size_ / 2 + 1);
        m_fspare = mem_allocate<float>(size_ + 2);
        m_fspare2 = mem_allocate<float>(size_ + 2);
    }

    void initDouble() {
        if (m_dspec) return;
        m_dspec = vDSP_create_fftsetupD(m_order, FFT_RADIX2);
        m_dbuf = new DSPDoubleSplitComplex;
        //!!! "If possible, tempBuffer->realp and tempBuffer->imagp should be 32-byte aligned for best performance."
        m_dbuf->realp = mem_allocate<double>(size_);
        m_dbuf->imagp = mem_allocate<double>(size_);
        m_dpacked = new DSPDoubleSplitComplex;
        m_dpacked->realp = mem_allocate<double>(size_ / 2 + 1);
        m_dpacked->imagp = mem_allocate<double>(size_ / 2 + 1);
        m_dspare = mem_allocate<double>(size_ + 2);
        m_dspare2 = mem_allocate<double>(size_ + 2);
    }

    void packReal(const float *R__ const real) {
        // Pack input for forward transform 
        vDSP_ctoz((DSPComplex *)real, 2, fpacked_, 1, size_/2);
    }
    void packComplex(const float *R__ const real, const float *R__ const imag) {
        // Pack input for inverse transform 
        if (real) vector_copy(fpacked_->realp, real, size_/2 + 1);
        else vector_zeros(fpacked_->realp, size_/2 + 1);
        if (imag) vector_copy(fpacked_->imagp, imag, size_/2 + 1);
        else vector_zeros(fpacked_->imagp, size_/2 + 1);
        fnyq();
    }

    void unpackReal(float *R__ const real) {
        // Unpack output for inverse transform
        vDSP_ztoc(fpacked_, 1, (DSPComplex *)real, 2, size_/2);
    }
    void unpackComplex(float *R__ const real, float *R__ const imag) {
        // Unpack output for forward transform
        // vDSP forward FFTs are scaled 2x (for some reason)
        float two = 2.f;
        vDSP_vsdiv(fpacked_->realp, 1, &two, real, 1, size_/2 + 1);
        vDSP_vsdiv(fpacked_->imagp, 1, &two, imag, 1, size_/2 + 1);
    }
    void unpackComplex(float *R__ const cplx) {
        // Unpack output for forward transform
        // vDSP forward FFTs are scaled 2x (for some reason)
        const int hs1 = size_/2 + 1;
        for (int i = 0; i < hs1; ++i) {
            cplx[i*2] = fpacked_->realp[i] / 2.f;
            cplx[i*2+1] = fpacked_->imagp[i] / 2.f;
        }
    }

    void packReal(const double *R__ const real) {
        // Pack input for forward transform
        vDSP_ctozD((DSPDoubleComplex *)real, 2, m_dpacked, 1, size_/2);
    }
    void packComplex(const double *R__ const real, const double *R__ const imag) {
        // Pack input for inverse transform
        if (real) vector_copy(m_dpacked->realp, real, size_/2 + 1);
        else vector_zeros(m_dpacked->realp, size_/2 + 1);
        if (imag) vector_copy(m_dpacked->imagp, imag, size_/2 + 1);
        else vector_zeros(m_dpacked->imagp, size_/2 + 1);
        dnyq();
    }

    void unpackReal(double *R__ const real) {
        // Unpack output for inverse transform
        vDSP_ztocD(m_dpacked, 1, (DSPDoubleComplex *)real, 2, size_/2);
    }
    void unpackComplex(double *R__ const real, double *R__ const imag) {
        // Unpack output for forward transform
        // vDSP forward FFTs are scaled 2x (for some reason)
        double two = 2.0;
        vDSP_vsdivD(m_dpacked->realp, 1, &two, real, 1, size_/2 + 1);
        vDSP_vsdivD(m_dpacked->imagp, 1, &two, imag, 1, size_/2 + 1);
    }
    void unpackComplex(double *R__ const cplx) {
        // Unpack output for forward transform
        // vDSP forward FFTs are scaled 2x (for some reason)
        const int hs1 = size_/2 + 1;
        for (int i = 0; i < hs1; ++i) {
            cplx[i*2] = m_dpacked->realp[i] / 2.0;
            cplx[i*2+1] = m_dpacked->imagp[i] / 2.0;
        }
    }

    void fdenyq() {
        // for fft result in packed form, unpack the DC and Nyquist bins
        const int hs = size_/2;
        fpacked_->realp[hs] = fpacked_->imagp[0];
        fpacked_->imagp[hs] = 0.f;
        fpacked_->imagp[0] = 0.f;
    }
    void ddenyq() {
        // for fft result in packed form, unpack the DC and Nyquist bins
        const int hs = size_/2;
        m_dpacked->realp[hs] = m_dpacked->imagp[0];
        m_dpacked->imagp[hs] = 0.;
        m_dpacked->imagp[0] = 0.;
    }

    void fnyq() {
        // for ifft input in packed form, pack the DC and Nyquist bins
        const int hs = size_/2;
        fpacked_->imagp[0] = fpacked_->realp[hs];
        fpacked_->realp[hs] = 0.f;
        fpacked_->imagp[hs] = 0.f;
    }
    void dnyq() {
        // for ifft input in packed form, pack the DC and Nyquist bins
        const int hs = size_/2;
        m_dpacked->imagp[0] = m_dpacked->realp[hs];
        m_dpacked->realp[hs] = 0.;
        m_dpacked->imagp[hs] = 0.;
    }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        //Profiler profiler("D_VDSP::forward [d]");
        if (!m_dspec) initDouble();
        packReal(realIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_FORWARD);
        ddenyq();
        unpackComplex(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        //Profiler profiler("D_VDSP::forward [d]");
        if (!m_dspec) initDouble();
        packReal(realIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_FORWARD);
        ddenyq();
        unpackComplex(complexOut);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        //Profiler profiler("D_VDSP::forwardPolar [d]");
        if (!m_dspec) initDouble();
        const int hs1 = size_/2+1;
        packReal(realIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_FORWARD);
        ddenyq();
        // vDSP forward FFTs are scaled 2x (for some reason)
        for (int i = 0; i < hs1; ++i) m_dpacked->realp[i] /= 2.0;
        for (int i = 0; i < hs1; ++i) m_dpacked->imagp[i] /= 2.0;
        vector_cartesian_to_polar(magOut, phaseOut,
                             m_dpacked->realp, m_dpacked->imagp, hs1);
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        //Profiler profiler("D_VDSP::forwardMagnitude [d]");
        if (!m_dspec) initDouble();
        packReal(realIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_FORWARD);
        ddenyq();
        const int hs1 = size_/2+1;
        vDSP_zvmagsD(m_dpacked, 1, m_dspare, 1, hs1);
        vvsqrt(m_dspare2, m_dspare, &hs1);
        // vDSP forward FFTs are scaled 2x (for some reason)
        double two = 2.0;
        vDSP_vsdivD(m_dspare2, 1, &two, magOut, 1, hs1);
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        //Profiler profiler("D_VDSP::forward [f]");
        if (!m_fspec) initFloat();
        packReal(realIn);
        vDSP_fft_zript(m_fspec, fpacked_, 1, fbuf_, m_order, FFT_FORWARD);
        fdenyq();
        unpackComplex(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        //Profiler profiler("D_VDSP::forward [f]");
        if (!m_fspec) initFloat();
        packReal(realIn);
        vDSP_fft_zript(m_fspec, fpacked_, 1, fbuf_, m_order, FFT_FORWARD);
        fdenyq();
        unpackComplex(complexOut);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        //Profiler profiler("D_VDSP::forwardPolar [f]");
        if (!m_fspec) initFloat();
        const int hs1 = size_/2+1;
        packReal(realIn);
        vDSP_fft_zript(m_fspec, fpacked_, 1, fbuf_, m_order, FFT_FORWARD);
        fdenyq();
        // vDSP forward FFTs are scaled 2x (for some reason)
        for (int i = 0; i < hs1; ++i) fpacked_->realp[i] /= 2.f;
        for (int i = 0; i < hs1; ++i) fpacked_->imagp[i] /= 2.f;
        vector_cartesian_to_polar(magOut, phaseOut,
                             fpacked_->realp, fpacked_->imagp, hs1);
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        //Profiler profiler("D_VDSP::forwardMagnitude [f]");
        if (!m_fspec) initFloat();
        packReal(realIn);
        vDSP_fft_zript(m_fspec, fpacked_, 1, fbuf_, m_order, FFT_FORWARD);
        fdenyq();
        const int hs1 = size_/2 + 1;
        vDSP_zvmags(fpacked_, 1, m_fspare, 1, hs1);
        vvsqrtf(m_fspare2, m_fspare, &hs1);
        // vDSP forward FFTs are scaled 2x (for some reason)
        float two = 2.f;
        vDSP_vsdiv(m_fspare2, 1, &two, magOut, 1, hs1);
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        //Profiler profiler("D_VDSP::inverse [d]");
        if (!m_dspec) initDouble();
        packComplex(realIn, imagIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        //Profiler profiler("D_VDSP::inverseInterleaved [d]");
        if (!m_dspec) initDouble();
        double *d[2] = { m_dpacked->realp, m_dpacked->imagp };
        vector_deinterleave(d, complexIn, 2, size_/2 + 1);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        //Profiler profiler("D_VDSP::inversePolar [d]");
        if (!m_dspec) initDouble();
        const int hs1 = size_/2+1;
        vvsincos(m_dpacked->imagp, m_dpacked->realp, phaseIn, &hs1);
        double *const rp = m_dpacked->realp;
        double *const ip = m_dpacked->imagp;
        for (int i = 0; i < hs1; ++i) rp[i] *= magIn[i];
        for (int i = 0; i < hs1; ++i) ip[i] *= magIn[i];
        dnyq();
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        //Profiler profiler("D_VDSP::inverseCepstral [d]");
        if (!m_dspec) initDouble();
        const int hs1 = size_/2 + 1;
        vector_copy(m_dspare, magIn, hs1);
        for (int i = 0; i < hs1; ++i) m_dspare[i] += 0.000001;
        vvlog(m_dspare2, m_dspare, &hs1);
        inverse(m_dspare2, 0, cepOut);
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        //Profiler profiler("D_VDSP::inverse [f]");
        if (!m_fspec) initFloat();
        packComplex(realIn, imagIn);
        vDSP_fft_zript(m_fspec, fpacked_, 1, fbuf_, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        //Profiler profiler("D_VDSP::inverseInterleaved [f]");
        if (!m_fspec) initFloat();
        float *f[2] = { fpacked_->realp, fpacked_->imagp };
        vector_deinterleave(f, complexIn, 2, size_/2 + 1);
        vDSP_fft_zript(m_fspec, fpacked_, 1, fbuf_, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        //Profiler profiler("D_VDSP::inversePolar [f]");
        if (!m_fspec) initFloat();

        const int hs1 = size_/2+1;
        vvsincosf(fpacked_->imagp, fpacked_->realp, phaseIn, &hs1);
        float *const rp = fpacked_->realp;
        float *const ip = fpacked_->imagp;
        for (int i = 0; i < hs1; ++i) rp[i] *= magIn[i];
        for (int i = 0; i < hs1; ++i) ip[i] *= magIn[i];
        fnyq();
        vDSP_fft_zript(m_fspec, fpacked_, 1, fbuf_, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        //Profiler profiler("D_VDSP::inverseCepstral [f]");
        if (!m_fspec) initFloat();
        const int hs1 = size_/2 + 1;
        vector_copy(m_fspare, magIn, hs1);
        for (int i = 0; i < hs1; ++i) m_fspare[i] += 0.000001f;
        vvlogf(m_fspare2, m_fspare, &hs1);
        inverse(m_fspare2, 0, cepOut);
    }

private:
    const int size_;
    int m_order;
    FFTSetup m_fspec;
    FFTSetupD m_dspec;
    DSPSplitComplex *fbuf_;
    DSPDoubleSplitComplex *m_dbuf;
    DSPSplitComplex *fpacked_;
    float *m_fspare;
    float *m_fspare2;
    DSPDoubleSplitComplex *m_dpacked;
    double *m_dspare;
    double *m_dspare2;
};

#endif /* HAVE_VDSP */

#ifdef HAVE_MEDIALIB

class D_MEDIALIB : public FFTImpl
{
public:
    D_MEDIALIB(int size) :
        size_(size),
        m_dpacked(0), fpacked_(0)
    { 
        for (int i = 0; ; ++i) {
            if (size_ & (1 << i)) {
                m_order = i;
                break;
            }
        }
    }

    ~D_MEDIALIB() {
        if (m_dpacked) {
            mem_deallocate(m_dpacked);
        }
        if (fpacked_) {
            mem_deallocate(fpacked_);
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision | FFT::DoublePrecision;
    }

    //!!! rv check

    void initFloat() {
        fpacked_ = mem_allocate<float>(size_*2);
    }

    void initDouble() {
        m_dpacked = mem_allocate<double>(size_*2);
    }

    void packFloatConjugates() {
        const int hs = size_ / 2;
        for (int i = 1; i <= hs; ++i) {
            fpacked_[(size_-i)*2] = fpacked_[2*i];
            fpacked_[(size_-i)*2 + 1] = -fpacked_[2*i + 1];
        }
    }

    void packDoubleConjugates() {
        const int hs = size_ / 2;
        for (int i = 1; i <= hs; ++i) {
            m_dpacked[(size_-i)*2] = m_dpacked[2*i];
            m_dpacked[(size_-i)*2 + 1] = -m_dpacked[2*i + 1];
        }
    }

    void packFloat(const float *R__ real, const float *R__ imag) {
        int index = 0;
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            fpacked_[index++] = real[i];
            index++;
        }
        index = 0;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                fpacked_[index++] = imag[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                fpacked_[index++] = 0.f;
            }
        }
        packFloatConjugates();
    }

    void packDouble(const double *R__ real, const double *R__ imag) {
        int index = 0;
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            m_dpacked[index++] = real[i];
            index++;
        }
        index = 0;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_dpacked[index++] = imag[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_dpacked[index++] = 0.0;
            }
        }
        packDoubleConjugates();
    }

    void unpackFloat(float *real, float *R__ imag) { // real may be equal to fpacked_
        int index = 0;
        const int hs = size_/2;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                imag[i] = fpacked_[index++];
            }
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            real[i] = fpacked_[index++];
            index++;
        }
    }        

    void unpackDouble(double *real, double *R__ imag) { // real may be equal to m_dpacked
        int index = 0;
        const int hs = size_/2;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                imag[i] = m_dpacked[index++];
            }
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            real[i] = m_dpacked[index++];
            index++;
        }
    }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        //Profiler profiler("D_MEDIALIB::forward [d]");
        if (!m_dpacked) initDouble();
        mlib_SignalFFT_1_D64C_D64(m_dpacked, realIn, m_order);
        unpackDouble(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        //Profiler profiler("D_MEDIALIB::forwardInterleaved [d]");
        if (!m_dpacked) initDouble();
        // mlib FFT gives the whole redundant complex result
        mlib_SignalFFT_1_D64C_D64(m_dpacked, realIn, m_order);
        vector_copy(complexOut, m_dpacked, size_ + 2);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        //Profiler profiler("D_MEDIALIB::forwardPolar [d]");
        if (!m_dpacked) initDouble();
        mlib_SignalFFT_1_D64C_D64(m_dpacked, realIn, m_order);
        const int hs = size_/2;
        int index = 0;
        for (int i = 0; i <= hs; ++i) {
            int reali = index;
            ++index;
            magOut[i] = sqrt(m_dpacked[reali] * m_dpacked[reali] +
                             m_dpacked[index] * m_dpacked[index]);
            phaseOut[i] = atan2(m_dpacked[index], m_dpacked[reali]) ;
            ++index;
        }
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        //Profiler profiler("D_MEDIALIB::forwardMagnitude [d]");
        if (!m_dpacked) initDouble();
        mlib_SignalFFT_1_D64C_D64(m_dpacked, realIn, m_order);
        const int hs = size_/2;
        int index = 0;
        for (int i = 0; i <= hs; ++i) {
            int reali = index;
            ++index;
            magOut[i] = sqrt(m_dpacked[reali] * m_dpacked[reali] +
                             m_dpacked[index] * m_dpacked[index]);
            ++index;
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        //Profiler profiler("D_MEDIALIB::forward [f]");
        if (!fpacked_) initFloat();
        mlib_SignalFFT_1_F32C_F32(fpacked_, realIn, m_order);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        //Profiler profiler("D_MEDIALIB::forwardInterleaved [f]");
        if (!fpacked_) initFloat();
        // mlib FFT gives the whole redundant complex result
        mlib_SignalFFT_1_F32C_F32(fpacked_, realIn, m_order);
        vector_copy(complexOut, fpacked_, size_ + 2);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        //Profiler profiler("D_MEDIALIB::forwardPolar [f]");
        if (!fpacked_) initFloat();
        mlib_SignalFFT_1_F32C_F32(fpacked_, realIn, m_order);
        const int hs = size_/2;
        int index = 0;
        for (int i = 0; i <= hs; ++i) {
            int reali = index;
            ++index;
            magOut[i] = sqrtf(fpacked_[reali] * fpacked_[reali] +
                              fpacked_[index] * fpacked_[index]);
            phaseOut[i] = atan2f(fpacked_[index], fpacked_[reali]);
            ++index;
        }
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        //Profiler profiler("D_MEDIALIB::forwardMagnitude [f]");
        if (!fpacked_) initFloat();
        mlib_SignalFFT_1_F32C_F32(fpacked_, realIn, m_order);
        const int hs = size_/2;
        int index = 0;
        for (int i = 0; i <= hs; ++i) {
            int reali = index;
            ++index;
            magOut[i] = sqrtf(fpacked_[reali] * fpacked_[reali] +
                              fpacked_[index] * fpacked_[index]);
            ++index;
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        //Profiler profiler("D_MEDIALIB::inverse [d]");
        if (!m_dpacked) initDouble();
        packDouble(realIn, imagIn);
        mlib_SignalIFFT_2_D64_D64C(realOut, m_dpacked, m_order);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        //Profiler profiler("D_MEDIALIB::inverseInterleaved [d]");
        if (!m_dpacked) initDouble();
        vector_copy(m_dpacked, complexIn, size_ + 2);
        packDoubleConjugates();
        mlib_SignalIFFT_2_D64_D64C(realOut, m_dpacked, m_order);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        //Profiler profiler("D_MEDIALIB::inversePolar [d]");
        if (!m_dpacked) initDouble();
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            double real = magIn[i] * cos(phaseIn[i]);
            double imag = magIn[i] * sin(phaseIn[i]);
            m_dpacked[i*2] = real;
            m_dpacked[i*2 + 1] = imag;
        }
        packDoubleConjugates();
        mlib_SignalIFFT_2_D64_D64C(realOut, m_dpacked, m_order);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        //Profiler profiler("D_MEDIALIB::inverseCepstral [d]");
        if (!m_dpacked) initDouble();
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            m_dpacked[i*2] = log(magIn[i] + 0.000001);
            m_dpacked[i*2 + 1] = 0.0;
        }
        packDoubleConjugates();
        mlib_SignalIFFT_2_D64_D64C(cepOut, m_dpacked, m_order);
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        //Profiler profiler("D_MEDIALIB::inverse [f]");
        if (!fpacked_) initFloat();
        packFloat(realIn, imagIn);
        mlib_SignalIFFT_2_F32_F32C(realOut, fpacked_, m_order);
    }
    
    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        //Profiler profiler("D_MEDIALIB::inverseInterleaved [f]");
        if (!fpacked_) initFloat();
        vector_assign(fpacked_, complexIn, size_ + 2);
        packFloatConjugates();
        mlib_SignalIFFT_2_F32_F32C(realOut, fpacked_, m_order);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        //Profiler profiler("D_MEDIALIB::inversePolar [f]");
        if (!fpacked_) initFloat();
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            double real = magIn[i] * cos(phaseIn[i]);
            double imag = magIn[i] * sin(phaseIn[i]);
            fpacked_[i*2] = real;
            fpacked_[i*2 + 1] = imag;
        }
        packFloatConjugates();
        mlib_SignalIFFT_2_F32_F32C(realOut, fpacked_, m_order);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        //Profiler profiler("D_MEDIALIB::inverseCepstral [f]");
        if (!fpacked_) initFloat();
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            fpacked_[i*2] = logf(magIn[i] + 0.000001);
            fpacked_[i*2 + 1] = 0.f;
        }
        packFloatConjugates();
        mlib_SignalIFFT_2_F32_F32C(cepOut, fpacked_, m_order);
    }

private:
    const int size_;
    int m_order;
    double *m_dpacked;
    float *fpacked_;
};

#endif /* HAVE_MEDIALIB */

#ifdef HAVE_OPENMAX

class D_OPENMAX : public FFTImpl
{
    // Convert a signed 32-bit integer to a float in the range [-1,1)
    static inline float i2f(OMX_S32 i)
    {
        return float(i) / float(OMX_MAX_S32);
    }

    // Convert a signed 32-bit integer to a double in the range [-1,1)
    static inline double i2d(OMX_S32 i)
    {
        return double(i) / double(OMX_MAX_S32);
    }

    // Convert a float in the range [-1,1) to a signed 32-bit integer
    static inline OMX_S32 f2i(float f)
    {
        return OMX_S32(f * OMX_MAX_S32);
    }

    // Convert a double in the range [-1,1) to a signed 32-bit integer
    static inline OMX_S32 d2i(double d)
    {
        return OMX_S32(d * OMX_MAX_S32);
    }

public:
    D_OPENMAX(int size) :
        size_(size),
        m_packed(0)
    { 
        for (int i = 0; ; ++i) {
            if (size_ & (1 << i)) {
                m_order = i;
                break;
            }
        }
    }

    ~D_OPENMAX() {
        if (m_packed) {
            mem_deallocate(m_packed);
            mem_deallocate(m_buf);
            mem_deallocate(fbuf_);
            mem_deallocate(m_spec);
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision;
    }

    //!!! rv check

    // The OpenMAX implementation uses a fixed-point representation in
    // 32-bit signed integers, with a downward scaling factor (0-32
    // bits) supplied as an argument to the FFT function.

    void initFloat() {
        initDouble();
    }

    void initDouble() {
        if (!m_packed) {
            m_buf = mem_allocate<OMX_S32>(size_);
            m_packed = mem_allocate<OMX_S32>(size_*2 + 2);
            fbuf_ = mem_allocate<float>(size_*2 + 2);
            OMX_INT sz = 0;
            omxSP_FFTGetBufSize_R_S32(m_order, &sz);
            m_spec = (OMXFFTSpec_R_S32 *)mem_allocate<char>(sz);
            omxSP_FFTInit_R_S32(m_spec, m_order);
        }
    }

    void packFloat(const float *R__ real) {
        // prepare fixed point input for forward transform
        for (int i = 0; i < size_; ++i) {
            m_buf[i] = f2i(real[i]);
        }
    }

    void packDouble(const double *R__ real) {
        // prepare fixed point input for forward transform
        for (int i = 0; i < size_; ++i) {
            m_buf[i] = d2i(real[i]);
        }
    }

    void unpackFloat(float *R__ real, float *R__ imag) {
        // convert fixed point output for forward transform
        int index = 0;
        const int hs = size_/2;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                imag[i] = i2f(m_packed[index++]);
            }
            vector_gain(imag, size_, hs + 1);
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            real[i] = i2f(m_packed[index++]);
            index++;
        }
        vector_gain(real, size_, hs + 1);
    }        

    void unpackDouble(double *R__ real, double *R__ imag) {
        // convert fixed point output for forward transform
        int index = 0;
        const int hs = size_/2;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                imag[i] = i2d(m_packed[index++]);
            }
            vector_gain(imag, size_, hs + 1);
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            real[i] = i2d(m_packed[index++]);
            index++;
        }
        vector_gain(real, size_, hs + 1);
    }

    void unpackFloatInterleaved(float *R__ cplx) {
        // convert fixed point output for forward transform
        for (int i = 0; i < size_ + 2; ++i) {
            cplx[i] = i2f(m_packed[i]);
        }            
        vector_gain(cplx, size_, size_ + 2);
    }

    void unpackDoubleInterleaved(double *R__ cplx) {
        // convert fixed point output for forward transform
        for (int i = 0; i < size_ + 2; ++i) {
            cplx[i] = i2d(m_packed[i]);
        }            
        vector_gain(cplx, size_, size_ + 2);
    }

    void packFloat(const float *R__ real, const float *R__ imag) {
        // prepare fixed point input for inverse transform
        int index = 0;
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            m_packed[index++] = f2i(real[i]);
            index++;
        }
        index = 0;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_packed[index++] = f2i(imag[i]);
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_packed[index++] = 0;
            }
        }
    }

    void packDouble(const double *R__ real, const double *R__ imag) {
        // prepare fixed point input for inverse transform
        int index = 0;
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            m_packed[index++] = d2i(real[i]);
            index++;
        }
        index = 0;
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_packed[index++] = d2i(imag[i]);
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_packed[index++] = 0;
            }
        }
    }

    void convertFloat(const float *R__ f) {
        // convert interleaved input for inverse interleaved transform
        const int n = size_ + 2;
        for (int i = 0; i < n; ++i) {
            m_packed[i] = f2i(f[i]);
        }
    }        

    void convertDouble(const double *R__ d) {
        // convert interleaved input for inverse interleaved transform
        const int n = size_ + 2;
        for (int i = 0; i < n; ++i) {
            m_packed[i] = d2i(d[i]);
        }
    }        

    void unpackFloat(float *R__ real) {
        // convert fixed point output for inverse transform
        for (int i = 0; i < size_; ++i) {
            real[i] = i2f(m_buf[i]) * size_;
        }
    }

    void unpackDouble(double *R__ real) {
        // convert fixed point output for inverse transform
        for (int i = 0; i < size_; ++i) {
            real[i] = i2d(m_buf[i]) * size_;
        }
    }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        //Profiler profiler("D_OPENMAX::forward [d]");
        if (!m_packed) initDouble();
        packDouble(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackDouble(realOut, imagOut);
    }
    
    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        //Profiler profiler("D_OPENMAX::forwardInterleaved [d]");
        if (!m_packed) initDouble();
        packDouble(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackDoubleInterleaved(complexOut);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        //Profiler profiler("D_OPENMAX::forwardPolar [d]");
        if (!m_packed) initDouble();
        packDouble(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackDouble(magOut, phaseOut); // temporarily
        // at this point we actually have real/imag in the mag/phase arrays
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            double real = magOut[i];
            double imag = phaseOut[i];
            complex_magphase(magOut + i, phaseOut + i, real, imag);
        }
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        //Profiler profiler("D_OPENMAX::forwardMagnitude [d]");
        if (!m_packed) initDouble();
        packDouble(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            int reali = i * 2;
            int imagi = reali + 1;
            double real = i2d(m_packed[reali]) * size_;
            double imag = i2d(m_packed[imagi]) * size_;
            magOut[i] = sqrt(real * real + imag * imag);
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        //Profiler profiler("D_OPENMAX::forward [f]");
        if (!m_packed) initFloat();
        packFloat(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        //Profiler profiler("D_OPENMAX::forwardInterleaved [f]");
        if (!m_packed) initFloat();
        packFloat(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackFloatInterleaved(complexOut);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        //Profiler profiler("D_OPENMAX::forwardPolar [f]");
        if (!m_packed) initFloat();

        packFloat(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackFloat(magOut, phaseOut); // temporarily
        // at this point we actually have real/imag in the mag/phase arrays
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            float real = magOut[i];
            float imag = phaseOut[i];
            complex_magphase(magOut + i, phaseOut + i, real, imag);
        }
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        //Profiler profiler("D_OPENMAX::forwardMagnitude [f]");
        if (!m_packed) initFloat();
        packFloat(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            int reali = i * 2;
            int imagi = reali + 1;
            float real = i2f(m_packed[reali]) * size_;
            float imag = i2f(m_packed[imagi]) * size_;
            magOut[i] = sqrtf(real * real + imag * imag);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        //Profiler profiler("D_OPENMAX::inverse [d]");
        if (!m_packed) initDouble();
        packDouble(realIn, imagIn);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackDouble(realOut);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        //Profiler profiler("D_OPENMAX::inverseInterleaved [d]");
        if (!m_packed) initDouble();
        convertDouble(complexIn);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackDouble(realOut);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        //Profiler profiler("D_OPENMAX::inversePolar [d]");
        if (!m_packed) initDouble();
        int index = 0;
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            double real, imag;
            complex_phasor(&real, &imag, phaseIn[i]);
            fbuf_[index++] = float(real);
            fbuf_[index++] = float(imag);
        }
        convertFloat(fbuf_);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackDouble(realOut);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        //Profiler profiler("D_OPENMAX::inverseCepstral [d]");
        if (!m_packed) initDouble();
        //!!! implement
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        //Profiler profiler("D_OPENMAX::inverse [f]");
        if (!m_packed) initFloat();
        packFloat(realIn, imagIn);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackFloat(realOut);
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        //Profiler profiler("D_OPENMAX::inverse [f]");
        if (!m_packed) initFloat();
        convertFloat(complexIn);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackFloat(realOut);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        //Profiler profiler("D_OPENMAX::inversePolar [f]");
        if (!m_packed) initFloat();
        const int hs = size_/2;
        vector_polar_to_cartesian_interleaved(fbuf_, magIn, phaseIn, hs+1);
        convertFloat(fbuf_);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackFloat(realOut);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        //Profiler profiler("D_OPENMAX::inverseCepstral [f]");
        if (!m_packed) initFloat();
        //!!! implement
    }

private:
    const int size_;
    int m_order;
    OMX_S32 *m_packed;
    OMX_S32 *m_buf;
    float *fbuf_;
    OMXFFTSpec_R_S32 *m_spec;

};

#endif /* HAVE_OPENMAX */

#ifdef HAVE_FFTW3

/*
 Define FFTW_DOUBLE_ONLY to make all uses of FFTW functions be
 double-precision (so "float" FFTs are calculated by casting to
 doubles and using the double-precision FFTW function).

 Define FFTW_SINGLE_ONLY to make all uses of FFTW functions be
 single-precision (so "double" FFTs are calculated by casting to
 floats and using the single-precision FFTW function).

 Neither of these flags is desirable for either performance or
 precision. The main reason to define either flag is to avoid linking
 against both fftw3 and fftw3f libraries.
*/

//#define FFTW_DOUBLE_ONLY 1
//#define FFTW_SINGLE_ONLY 1

#if defined(FFTW_DOUBLE_ONLY) && defined(FFTW_SINGLE_ONLY)
// Can't meaningfully define both
#error Can only define one of FFTW_DOUBLE_ONLY and FFTW_SINGLE_ONLY
#endif

#if defined(FFTW_FLOAT_ONLY)
#warning FFTW_FLOAT_ONLY is deprecated, use FFTW_SINGLE_ONLY instead
#define FFTW_SINGLE_ONLY 1
#endif

#ifdef FFTW_DOUBLE_ONLY
#define fft_float_type double
#define fftwf_complex fftw_complex
#define fftwf_plan fftw_plan
#define fftwf_plan_dft_r2c_1d fftw_plan_dft_r2c_1d
#define fftwf_plan_dft_c2r_1d fftw_plan_dft_c2r_1d
#define fftwf_destroy_plan fftw_destroy_plan
#define fftwf_malloc fftw_malloc
#define fftwf_free fftw_free
#define fftwf_cleanup fftw_cleanup
#define fftwf_execute fftw_execute
#define atan2f atan2
#define sqrtf sqrt
#define cosf cos
#define sinf sin
#else
#define fft_float_type float
#endif /* FFTW_DOUBLE_ONLY */

#ifdef FFTW_SINGLE_ONLY
#define fft_double_type float
#define fftw_complex fftwf_complex
#define fftw_plan fftwf_plan
#define fftw_plan_dft_r2c_1d fftwf_plan_dft_r2c_1d
#define fftw_plan_dft_c2r_1d fftwf_plan_dft_c2r_1d
#define fftw_destroy_plan fftwf_destroy_plan
#define fftw_malloc fftwf_malloc
#define fftw_free fftwf_free
#define fftw_cleanup fftwf_cleanup
#define fftw_execute fftwf_execute
#define atan2 atan2f
#define sqrt sqrtf
#define cos cosf
#define sin sinf
#else
#define fft_double_type double
#endif /* FFTW_SINGLE_ONLY */

class D_FFTW : public FFTImpl
{
public:
    D_FFTW(int size) :
        fplanf_(0), m_dplanf(0), size_(size)
    {
    }

    ~D_FFTW() {
        if (fplanf_) {
#ifndef NO_THREADING
            m_commonMutex.lock();
#endif
            bool save = false;
            if (m_extantf > 0 && --m_extantf == 0) save = true;
            (void)save;
#ifndef FFTW_DOUBLE_ONLY
            if (save) saveWisdom('f');
#endif
            fftwf_destroy_plan(fplanf_);
            fftwf_destroy_plan(fplani_);
            fftwf_free(fbuf_);
            fftwf_free(fpacked_);
#ifndef NO_THREADING
            m_commonMutex.unlock();
#endif
        }
        if (m_dplanf) {
#ifndef NO_THREADING
            m_commonMutex.lock();
#endif
            bool save = false;
            if (m_extantd > 0 && --m_extantd == 0) save = true;
            (void)save;
#ifndef FFTW_SINGLE_ONLY
            if (save) saveWisdom('d');
#endif
            fftw_destroy_plan(m_dplanf);
            fftw_destroy_plan(m_dplani);
            fftw_free(m_dbuf);
            fftw_free(m_dpacked);
#ifndef NO_THREADING
            m_commonMutex.unlock();
#endif
        }
#ifndef NO_THREADING
        m_commonMutex.lock();
#endif
        if (m_extantf <= 0 && m_extantd <= 0) {
            fftw_cleanup();
        }
#ifndef NO_THREADING
        m_commonMutex.unlock();
#endif
    }

    FFT::Precisions
    getSupportedPrecisions() const {
#ifdef FFTW_SINGLE_ONLY
        return FFT::SinglePrecision;
#else
#ifdef FFTW_DOUBLE_ONLY
        return FFT::DoublePrecision;
#else
        return FFT::SinglePrecision | FFT::DoublePrecision;
#endif
#endif
    }

    void initFloat() {
        if (fplanf_) return;
        bool load = false;
#ifndef NO_THREADING
        m_commonMutex.lock();
#endif
        if (m_extantf++ == 0) load = true;
#ifdef FFTW_DOUBLE_ONLY
        if (load) loadWisdom('d');
#else
        if (load) loadWisdom('f');
#endif
        fbuf_ = (fft_float_type *)fftw_malloc(size_ * sizeof(fft_float_type));
        fpacked_ = (fftwf_complex *)fftw_malloc
            ((size_/2 + 1) * sizeof(fftwf_complex));
        fplanf_ = fftwf_plan_dft_r2c_1d
            (size_, fbuf_, fpacked_, FFTW_MEASURE);
        fplani_ = fftwf_plan_dft_c2r_1d
            (size_, fpacked_, fbuf_, FFTW_MEASURE);
#ifndef NO_THREADING
        m_commonMutex.unlock();
#endif
    }

    void initDouble() {
        if (m_dplanf) return;
        bool load = false;
#ifndef NO_THREADING
        m_commonMutex.lock();
#endif
        if (m_extantd++ == 0) load = true;
#ifdef FFTW_SINGLE_ONLY
        if (load) loadWisdom('f');
#else
        if (load) loadWisdom('d');
#endif
        m_dbuf = (fft_double_type *)fftw_malloc(size_ * sizeof(fft_double_type));
        m_dpacked = (fftw_complex *)fftw_malloc
            ((size_/2 + 1) * sizeof(fftw_complex));
        m_dplanf = fftw_plan_dft_r2c_1d
            (size_, m_dbuf, m_dpacked, FFTW_MEASURE);
        m_dplani = fftw_plan_dft_c2r_1d
            (size_, m_dpacked, m_dbuf, FFTW_MEASURE);
#ifndef NO_THREADING
        m_commonMutex.unlock();
#endif
    }

    void loadWisdom(char type) { wisdom(false, type); }
    void saveWisdom(char type) { wisdom(true, type); }

    void wisdom(bool save, char type) {

#ifdef FFTW_DOUBLE_ONLY
        if (type == 'f') return;
#endif
#ifdef FFTW_SINGLE_ONLY
        if (type == 'd') return;
#endif

        const char *home = getenv("HOME");
        if (!home) return;

        char fn[256];
        snprintf(fn, 256, "%s/%s.%c", home, ".common.wisdom", type);

        FILE *f = fopen(fn, save ? "wb" : "rb");
        if (!f) return;

        if (save) {
            switch (type) {
#ifdef FFTW_DOUBLE_ONLY
            case 'f': break;
#else
            case 'f': fftwf_export_wisdom_to_file(f); break;
#endif
#ifdef FFTW_SINGLE_ONLY
            case 'd': break;
#else
            case 'd': fftw_export_wisdom_to_file(f); break;
#endif
            default: break;
            }
        } else {
            switch (type) {
#ifdef FFTW_DOUBLE_ONLY
            case 'f': break;
#else
            case 'f': fftwf_import_wisdom_from_file(f); break;
#endif
#ifdef FFTW_SINGLE_ONLY
            case 'd': break;
#else
            case 'd': fftw_import_wisdom_from_file(f); break;
#endif
            default: break;
            }
        }

        fclose(f);
    }

    void packFloat(const float *R__ real, const float *R__ imag) {
        const int hs = size_/2;
        fftwf_complex *const R__ fpacked = fpacked_; 
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][0] = real[i];
        }
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                fpacked[i][1] = imag[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                fpacked[i][1] = 0.f;
            }
        }                
    }

    void packDouble(const double *R__ real, const double *R__ imag) {
        const int hs = size_/2;
        fftw_complex *const R__ dpacked = m_dpacked; 
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][0] = real[i];
        }
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                dpacked[i][1] = imag[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                dpacked[i][1] = 0.0;
            }
        }
    }

    void unpackFloat(float *R__ real, float *R__ imag) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            real[i] = fpacked_[i][0];
        }
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                imag[i] = fpacked_[i][1];
            }
        }
    }        

    void unpackDouble(double *R__ real, double *R__ imag) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            real[i] = m_dpacked[i][0];
        }
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                imag[i] = m_dpacked[i][1];
            }
        }
    }        

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        if (!m_dplanf) initDouble();
        const int sz = size_;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realIn != dbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                dbuf[i] = realIn[i];
            }
        fftw_execute(m_dplanf);
        unpackDouble(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        if (!m_dplanf) initDouble();
        const int sz = size_;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realIn != dbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                dbuf[i] = realIn[i];
            }
        fftw_execute(m_dplanf);
        vector_assign(complexOut, (fft_double_type *)m_dpacked, sz + 2);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        if (!m_dplanf) initDouble();
        fft_double_type *const R__ dbuf = m_dbuf;
        const int sz = size_;
#ifndef FFTW_SINGLE_ONLY
        if (realIn != dbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                dbuf[i] = realIn[i];
            }
        fftw_execute(m_dplanf);
        vector_cartesian_interleaved_to_polar(magOut, phaseOut,
                                         (double *)m_dpacked, size_/2+1);
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        if (!m_dplanf) initDouble();
        fft_double_type *const R__ dbuf = m_dbuf;
        const int sz = size_;
#ifndef FFTW_SINGLE_ONLY
        if (realIn != m_dbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                dbuf[i] = realIn[i];
            }
        fftw_execute(m_dplanf);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_dpacked[i][0] * m_dpacked[i][0] +
                             m_dpacked[i][1] * m_dpacked[i][1]);
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        if (!fplanf_) initFloat();
        fft_float_type *const R__ fbuf = fbuf_;
        const int sz = size_;
#ifndef FFTW_DOUBLE_ONLY
        if (realIn != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                fbuf[i] = realIn[i];
            }
        fftwf_execute(fplanf_);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        if (!fplanf_) initFloat();
        fft_float_type *const R__ fbuf = fbuf_;
        const int sz = size_;
#ifndef FFTW_DOUBLE_ONLY
        if (realIn != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                fbuf[i] = realIn[i];
            }
        fftwf_execute(fplanf_);
        vector_assign(complexOut, (fft_float_type *)fpacked_, sz + 2);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        if (!fplanf_) initFloat();
        fft_float_type *const R__ fbuf = fbuf_;
        const int sz = size_;
#ifndef FFTW_DOUBLE_ONLY
        if (realIn != fbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                fbuf[i] = realIn[i];
            }
        fftwf_execute(fplanf_);
        vector_cartesian_interleaved_to_polar(magOut, phaseOut,
                                         (float *)fpacked_, size_/2+1);
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        if (!fplanf_) initFloat();
        fft_float_type *const R__ fbuf = fbuf_;
        const int sz = size_;
#ifndef FFTW_DOUBLE_ONLY
        if (realIn != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                fbuf[i] = realIn[i];
            }
        fftwf_execute(fplanf_);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrtf(fpacked_[i][0] * fpacked_[i][0] +
                              fpacked_[i][1] * fpacked_[i][1]);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, imagIn);
        fftw_execute(m_dplani);
        const int sz = size_;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realOut != dbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = dbuf[i];
            }
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        vector_assign((double *)m_dpacked, complexIn, size_ + 2);
        fftw_execute(m_dplani);
        const int sz = size_;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realOut != dbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = dbuf[i];
            }
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        const int hs = size_/2;
        fftw_complex *const R__ dpacked = m_dpacked;
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][0] = magIn[i] * cos(phaseIn[i]);
        }
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][1] = magIn[i] * sin(phaseIn[i]);
        }
        fftw_execute(m_dplani);
        const int sz = size_;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realOut != dbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = dbuf[i];
            }
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        if (!m_dplanf) initDouble();
        fft_double_type *const R__ dbuf = m_dbuf;
        fftw_complex *const R__ dpacked = m_dpacked;
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][0] = log(magIn[i] + 0.000001);
        }
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][1] = 0.0;
        }
        fftw_execute(m_dplani);
        const int sz = size_;
#ifndef FFTW_SINGLE_ONLY
        if (cepOut != dbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                cepOut[i] = dbuf[i];
            }
    }

    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        if (!fplanf_) initFloat();
        packFloat(realIn, imagIn);
        fftwf_execute(fplani_);
        const int sz = size_;
        fft_float_type *const R__ fbuf = fbuf_;
#ifndef FFTW_DOUBLE_ONLY
        if (realOut != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = fbuf[i];
            }
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        if (!fplanf_) initFloat();
        vector_copy((float *)fpacked_, complexIn, size_ + 2);
        fftwf_execute(fplani_);
        const int sz = size_;
        fft_float_type *const R__ fbuf = fbuf_;
#ifndef FFTW_DOUBLE_ONLY
        if (realOut != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = fbuf[i];
            }
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        if (!fplanf_) initFloat();
        const int hs = size_/2;
        fftwf_complex *const R__ fpacked = fpacked_;
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][0] = magIn[i] * cosf(phaseIn[i]);
        }
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][1] = magIn[i] * sinf(phaseIn[i]);
        }
        fftwf_execute(fplani_);
        const int sz = size_;
        fft_float_type *const R__ fbuf = fbuf_;
#ifndef FFTW_DOUBLE_ONLY
        if (realOut != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = fbuf[i];
            }
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        if (!fplanf_) initFloat();
        const int hs = size_/2;
        fftwf_complex *const R__ fpacked = fpacked_;
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][0] = logf(magIn[i] + 0.000001f);
        }
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][1] = 0.f;
        }
        fftwf_execute(fplani_);
        const int sz = size_;
        fft_float_type *const R__ fbuf = fbuf_;
#ifndef FFTW_DOUBLE_ONLY
        if (cepOut != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                cepOut[i] = fbuf[i];
            }
    }

private:
    fftwf_plan fplanf_;
    fftwf_plan fplani_;
#ifdef FFTW_DOUBLE_ONLY
    double *fbuf_;
#else
    float *fbuf_;
#endif
    fftwf_complex *fpacked_;
    fftw_plan m_dplanf;
    fftw_plan m_dplani;
#ifdef FFTW_SINGLE_ONLY
    float *m_dbuf;
#else
    double *m_dbuf;
#endif
    fftw_complex *m_dpacked;
    const int size_;
    static int m_extantf;
    static int m_extantd;
#ifndef NO_THREADING
    static Mutex m_commonMutex;
#endif
};

int
D_FFTW::m_extantf = 0;

int
D_FFTW::m_extantd = 0;

#ifndef NO_THREADING
Mutex
D_FFTW::m_commonMutex;
#endif

#endif /* HAVE_FFTW3 */

#ifdef HAVE_SFFT

/*
 Define SFFT_DOUBLE_ONLY to make all uses of SFFT functions be
 double-precision (so "float" FFTs are calculated by casting to
 doubles and using the double-precision SFFT function).

 Define SFFT_SINGLE_ONLY to make all uses of SFFT functions be
 single-precision (so "double" FFTs are calculated by casting to
 floats and using the single-precision SFFT function).

 Neither of these flags is desirable for either performance or
 precision.
*/

//#define SFFT_DOUBLE_ONLY 1
//#define SFFT_SINGLE_ONLY 1

#if defined(SFFT_DOUBLE_ONLY) && defined(SFFT_SINGLE_ONLY)
// Can't meaningfully define both
#error Can only define one of SFFT_DOUBLE_ONLY and SFFT_SINGLE_ONLY
#endif

#ifdef SFFT_DOUBLE_ONLY
#define fft_float_type double
#define FLAG_SFFT_FLOAT SFFT_DOUBLE
#define atan2f atan2
#define sqrtf sqrt
#define cosf cos
#define sinf sin
#define logf log
#else
#define FLAG_SFFT_FLOAT SFFT_FLOAT
#define fft_float_type float
#endif /* SFFT_DOUBLE_ONLY */

#ifdef SFFT_SINGLE_ONLY
#define fft_double_type float
#define FLAG_SFFT_DOUBLE SFFT_FLOAT
#define atan2 atan2f
#define sqrt sqrtf
#define cos cosf
#define sin sinf
#define log logf
#else
#define FLAG_SFFT_DOUBLE SFFT_DOUBLE
#define fft_double_type double
#endif /* SFFT_SINGLE_ONLY */

class D_SFFT : public FFTImpl
{
public:
    D_SFFT(int size) :
        fplanf_(0), fplani_(0), m_dplanf(0), m_dplani(0), size_(size)
    {
    }

    ~D_SFFT() {
        if (fplanf_) {
            sfft_free(fplanf_);
            sfft_free(fplani_);
            mem_deallocate(fbuf_);
            mem_deallocate(m_fresult);
        }
        if (m_dplanf) {
            sfft_free(m_dplanf);
            sfft_free(m_dplani);
            mem_deallocate(m_dbuf);
            mem_deallocate(m_dresult);
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
#ifdef SFFT_SINGLE_ONLY
        return FFT::SinglePrecision;
#else
#ifdef SFFT_DOUBLE_ONLY
        return FFT::DoublePrecision;
#else
        return FFT::SinglePrecision | FFT::DoublePrecision;
#endif
#endif
    }

    void initFloat() {
        if (fplanf_) return;
        fbuf_ = mem_allocate<fft_float_type>(2 * size_);
        m_fresult = mem_allocate<fft_float_type>(2 * size_);
        fplanf_ = sfft_init(size_, SFFT_FORWARD | FLAG_SFFT_FLOAT);
        fplani_ = sfft_init(size_, SFFT_BACKWARD | FLAG_SFFT_FLOAT);
        if (!fplanf_ || !fplani_) {
            if (!fplanf_) {
                std::cerr << "D_SFFT: Failed to construct forward float transform for size " << size_ << " (check SFFT library's target configuration)" << std::endl;
            } else {
                std::cerr << "D_SFFT: Failed to construct inverse float transform for size " << size_ << " (check SFFT library's target configuration)" << std::endl;
            }
#ifndef NO_EXCEPTIONS
            throw FFT::InternalError;
#else
            abort();
#endif
        }
    }

    void initDouble() {
        if (m_dplanf) return;
        m_dbuf = mem_allocate<fft_double_type>(2 * size_);
        m_dresult = mem_allocate<fft_double_type>(2 * size_);
        m_dplanf = sfft_init(size_, SFFT_FORWARD | FLAG_SFFT_DOUBLE);
        m_dplani = sfft_init(size_, SFFT_BACKWARD | FLAG_SFFT_DOUBLE);
        if (!m_dplanf || !m_dplani) {
            if (!m_dplanf) {
                std::cerr << "D_SFFT: Failed to construct forward double transform for size " << size_ << " (check SFFT library's target configuration)" << std::endl;
            } else {
                std::cerr << "D_SFFT: Failed to construct inverse double transform for size " << size_ << " (check SFFT library's target configuration)" << std::endl;
            }
#ifndef NO_EXCEPTIONS
            throw FFT::InternalError;
#else
            abort();
#endif
        }
    }

    void packFloat(const float *R__ real, const float *R__ imag, fft_float_type *target, int n) {
        for (int i = 0; i < n; ++i) target[i*2] = real[i];
        if (imag) {
            for (int i = 0; i < n; ++i) target[i*2+1] = imag[i]; 
        } else {
            for (int i = 0; i < n; ++i) target[i*2+1] = 0.f;
        }                
    }

    void packDouble(const double *R__ real, const double *R__ imag, fft_double_type *target, int n) {
        for (int i = 0; i < n; ++i) target[i*2] = real[i];
        if (imag) {
            for (int i = 0; i < n; ++i) target[i*2+1] = imag[i];
        } else {
            for (int i = 0; i < n; ++i) target[i*2+1] = 0.0;
        }                
    }

    void unpackFloat(const fft_float_type *source, float *R__ real, float *R__ imag, int n) {
        for (int i = 0; i < n; ++i) real[i] = source[i*2];
        if (imag) {
            for (int i = 0; i < n; ++i) imag[i] = source[i*2+1];
        }
    }        

    void unpackDouble(const fft_double_type *source, double *R__ real, double *R__ imag, int n) {
        for (int i = 0; i < n; ++i) real[i] = source[i*2];
        if (imag) {
            for (int i = 0; i < n; ++i) imag[i] = source[i*2+1];
        }
    }        

    template<typename T>
    void mirror(T *R__ cplx, int n) {
        for (int i = 1; i <= n/2; ++i) {
            int j = n-i;
            cplx[j*2] = cplx[i*2];
            cplx[j*2+1] = -cplx[i*2+1];
        }
    }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, 0, m_dbuf, size_);
        sfft_execute(m_dplanf, m_dbuf, m_dresult);
        unpackDouble(m_dresult, realOut, imagOut, size_/2+1);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, 0, m_dbuf, size_);
        sfft_execute(m_dplanf, m_dbuf, m_dresult);
        vector_assign(complexOut, m_dresult, size_+2); // i.e. size_/2+1 complex
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, 0, m_dbuf, size_);
        sfft_execute(m_dplanf, m_dbuf, m_dresult);
        vector_cartesian_interleaved_to_polar(magOut, phaseOut,
                                         m_dresult, size_/2+1);
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, 0, m_dbuf, size_);
        sfft_execute(m_dplanf, m_dbuf, m_dresult);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_dresult[i*2] * m_dresult[i*2] +
                             m_dresult[i*2+1] * m_dresult[i*2+1]);
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        if (!fplanf_) initFloat();
        packFloat(realIn, 0, fbuf_, size_);
        sfft_execute(fplanf_, fbuf_, m_fresult);
        unpackFloat(m_fresult, realOut, imagOut, size_/2+1);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        if (!fplanf_) initFloat();
        packFloat(realIn, 0, fbuf_, size_);
        sfft_execute(fplanf_, fbuf_, m_fresult);
        vector_assign(complexOut, m_fresult, size_+2); // i.e. size_/2+1 complex
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        if (!fplanf_) initFloat();
        packFloat(realIn, 0, fbuf_, size_);
        sfft_execute(fplanf_, fbuf_, m_fresult);
        vector_cartesian_interleaved_to_polar(magOut, phaseOut,
                                         m_fresult, size_/2+1);
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        if (!fplanf_) initFloat();
        packFloat(realIn, 0, fbuf_, size_);
        sfft_execute(fplanf_, fbuf_, m_fresult);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrtf(m_fresult[i*2] * m_fresult[i*2] +
                              m_fresult[i*2+1] * m_fresult[i*2+1]);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, imagIn, m_dbuf, size_/2+1);
        mirror(m_dbuf, size_);
        sfft_execute(m_dplani, m_dbuf, m_dresult);
        for (int i = 0; i < size_; ++i) {
            realOut[i] = m_dresult[i*2];
        }
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        vector_assign((double *)m_dbuf, complexIn, size_ + 2);
        mirror(m_dbuf, size_);
        sfft_execute(m_dplani, m_dbuf, m_dresult);
        for (int i = 0; i < size_; ++i) {
            realOut[i] = m_dresult[i*2];
        }
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            m_dbuf[i*2] = magIn[i] * cos(phaseIn[i]);
            m_dbuf[i*2+1] = magIn[i] * sin(phaseIn[i]);
        }
        mirror(m_dbuf, size_);
        sfft_execute(m_dplani, m_dbuf, m_dresult);
        for (int i = 0; i < size_; ++i) {
            realOut[i] = m_dresult[i*2];
        }
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        if (!m_dplanf) initDouble();
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            m_dbuf[i*2] = log(magIn[i] + 0.000001);
            m_dbuf[i*2+1] = 0.0;
        }
        mirror(m_dbuf, size_);
        sfft_execute(m_dplani, m_dbuf, m_dresult);
        for (int i = 0; i < size_; ++i) {
            cepOut[i] = m_dresult[i*2];
        }
    }

    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        if (!fplanf_) initFloat();
        packFloat(realIn, imagIn, fbuf_, size_/2+1);
        mirror(fbuf_, size_);
        sfft_execute(fplani_, fbuf_, m_fresult);
        for (int i = 0; i < size_; ++i) {
            realOut[i] = m_fresult[i*2];
        }
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        if (!fplanf_) initFloat();
        vector_assign((float *)fbuf_, complexIn, size_ + 2);
        mirror(fbuf_, size_);
        sfft_execute(fplani_, fbuf_, m_fresult);
        for (int i = 0; i < size_; ++i) {
            realOut[i] = m_fresult[i*2];
        }
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        if (!fplanf_) initFloat();
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            fbuf_[i*2] = magIn[i] * cosf(phaseIn[i]);
            fbuf_[i*2+1] = magIn[i] * sinf(phaseIn[i]);
        }
        mirror(fbuf_, size_);
        sfft_execute(fplani_, fbuf_, m_fresult);
        for (int i = 0; i < size_; ++i) {
            realOut[i] = m_fresult[i*2];
        }
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        if (!fplanf_) initFloat();
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            fbuf_[i*2] = logf(magIn[i] + 0.00001);
            fbuf_[i*2+1] = 0.0f;
        }
        sfft_execute(fplani_, fbuf_, m_fresult);
        for (int i = 0; i < size_; ++i) {
            cepOut[i] = m_fresult[i*2];
        }
    }

private:
    sfft_plan_t *fplanf_;
    sfft_plan_t *fplani_;
    fft_float_type *fbuf_;
    fft_float_type *m_fresult;

    sfft_plan_t *m_dplanf;
    sfft_plan_t *m_dplani;
    fft_double_type *m_dbuf;
    fft_double_type *m_dresult;

    const int size_;
};

#endif /* HAVE_SFFT */

#ifdef USE_KISSFFT

class D_KISSFFT : public FFTImpl
{
public:
    D_KISSFFT(int size) :
        size_(size),
        fplanf_(0),  
        fplani_(0)
    {
#ifdef FIXED_POINT
#error KISSFFT is not configured for float values
#endif
        if (sizeof(kiss_fft_scalar) != sizeof(float)) {
            std::cerr << "ERROR: KISSFFT is not configured for float values"
                      << std::endl;
        }

        fbuf_ = new kiss_fft_scalar[size_ + 2];
        fpacked_ = new kiss_fft_cpx[size_ + 2];
        fplanf_ = kiss_fftr_alloc(size_, 0, NULL, NULL);
        fplani_ = kiss_fftr_alloc(size_, 1, NULL, NULL);
    }

    ~D_KISSFFT() {
        kiss_fftr_free(fplanf_);
        kiss_fftr_free(fplani_);
        kiss_fft_cleanup();

        delete[] fbuf_;
        delete[] fpacked_;
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision;
    }

    void initFloat() { }
    void initDouble() { }

    void packFloat(const float *R__ real, const float *R__ imag) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            fpacked_[i].r = real[i];
        }
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                fpacked_[i].i = imag[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                fpacked_[i].i = 0.f;
            }
        }
    }

    void unpackFloat(float *R__ real, float *R__ imag) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            real[i] = fpacked_[i].r;
        }
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                imag[i] = fpacked_[i].i;
            }
        }
    }        

    void packDouble(const double *R__ real, const double *R__ imag) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            fpacked_[i].r = float(real[i]);
        }
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                fpacked_[i].i = float(imag[i]);
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                fpacked_[i].i = 0.f;
            }
        }
    }

    void unpackDouble(double *R__ real, double *R__ imag) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            real[i] = double(fpacked_[i].r);
        }
        if (imag) {
            for (int i = 0; i <= hs; ++i) {
                imag[i] = double(fpacked_[i].i);
            }
        }
    }        

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {

        vector_assign(fbuf_, realIn, size_);
        kiss_fftr(fplanf_, fbuf_, fpacked_);
        unpackDouble(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {

        vector_assign(fbuf_, realIn, size_);
        kiss_fftr(fplanf_, fbuf_, fpacked_);
        vector_assign(complexOut, (float *)fpacked_, size_ + 2);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {

        for (int i = 0; i < size_; ++i) {
            fbuf_[i] = float(realIn[i]);
        }

        kiss_fftr(fplanf_, fbuf_, fpacked_);

        const int hs = size_/2;

        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(double(fpacked_[i].r) * double(fpacked_[i].r) +
                             double(fpacked_[i].i) * double(fpacked_[i].i));
        }

        for (int i = 0; i <= hs; ++i) {
            phaseOut[i] = atan2(double(fpacked_[i].i), double(fpacked_[i].r));
        }
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {

        for (int i = 0; i < size_; ++i) {
            fbuf_[i] = float(realIn[i]);
        }

        kiss_fftr(fplanf_, fbuf_, fpacked_);

        const int hs = size_/2;

        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(double(fpacked_[i].r) * double(fpacked_[i].r) +
                             double(fpacked_[i].i) * double(fpacked_[i].i));
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {

        kiss_fftr(fplanf_, realIn, fpacked_);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {

        kiss_fftr(fplanf_, realIn, (kiss_fft_cpx *)complexOut);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {

        kiss_fftr(fplanf_, realIn, fpacked_);

        const int hs = size_/2;

        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrtf(fpacked_[i].r * fpacked_[i].r +
                              fpacked_[i].i * fpacked_[i].i);
        }

        for (int i = 0; i <= hs; ++i) {
            phaseOut[i] = atan2f(fpacked_[i].i, fpacked_[i].r);
        }
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {

        kiss_fftr(fplanf_, realIn, fpacked_);

        const int hs = size_/2;

        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrtf(fpacked_[i].r * fpacked_[i].r +
                              fpacked_[i].i * fpacked_[i].i);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {

        packDouble(realIn, imagIn);

        kiss_fftri(fplani_, fpacked_, fbuf_);

        for (int i = 0; i < size_; ++i) {
            realOut[i] = fbuf_[i];
        }
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {

        vector_assign((float *)fpacked_, complexIn, size_ + 2);

        kiss_fftri(fplani_, fpacked_, fbuf_);

        for (int i = 0; i < size_; ++i) {
            realOut[i] = fbuf_[i];
        }
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {

        const int hs = size_/2;

        for (int i = 0; i <= hs; ++i) {
            fpacked_[i].r = float(magIn[i] * cos(phaseIn[i]));
            fpacked_[i].i = float(magIn[i] * sin(phaseIn[i]));
        }

        kiss_fftri(fplani_, fpacked_, fbuf_);

        for (int i = 0; i < size_; ++i) {
            realOut[i] = fbuf_[i];
        }
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {

        const int hs = size_/2;

        for (int i = 0; i <= hs; ++i) {
            fpacked_[i].r = float(log(magIn[i] + 0.000001));
            fpacked_[i].i = 0.0f;
        }

        kiss_fftri(fplani_, fpacked_, fbuf_);

        for (int i = 0; i < size_; ++i) {
            cepOut[i] = fbuf_[i];
        }
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {

        packFloat(realIn, imagIn);
        kiss_fftri(fplani_, fpacked_, realOut);
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {

        vector_copy((float *)fpacked_, complexIn, size_ + 2);
        kiss_fftri(fplani_, fpacked_, realOut);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {

        const int hs = size_/2;

        for (int i = 0; i <= hs; ++i) {
            fpacked_[i].r = magIn[i] * cosf(phaseIn[i]);
            fpacked_[i].i = magIn[i] * sinf(phaseIn[i]);
        }

        kiss_fftri(fplani_, fpacked_, realOut);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {

        const int hs = size_/2;

        for (int i = 0; i <= hs; ++i) {
            fpacked_[i].r = logf(magIn[i] + 0.000001f);
            fpacked_[i].i = 0.0f;
        }

        kiss_fftri(fplani_, fpacked_, cepOut);
    }

private:
    const int size_;
    kiss_fftr_cfg fplanf_;
    kiss_fftr_cfg fplani_;
    kiss_fft_scalar *fbuf_;
    kiss_fft_cpx *fpacked_;
};

#endif /* USE_KISSFFT */

#ifdef USE_BUILTIN_FFT

class D_Cross : public FFTImpl
{
public:
    D_Cross(int size) : size_(size), m_table(0) {
        
        m_a = new double[size];
        m_b = new double[size];
        m_c = new double[size];
        m_d = new double[size];

        m_table = new int[size_];
    
        int bits;
        int i, j, k, m;

        for (i = 0; ; ++i) {
            if (size_ & (1 << i)) {
                bits = i;
                break;
            }
        }
        
        for (i = 0; i < size_; ++i) {
            
            m = i;
            
            for (j = k = 0; j < bits; ++j) {
                k = (k << 1) | (m & 1);
                m >>= 1;
            }
            
            m_table[i] = k;
        }
    }

    ~D_Cross() {
        delete[] m_table;
        delete[] m_a;
        delete[] m_b;
        delete[] m_c;
        delete[] m_d;
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::DoublePrecision;
    }

    void initFloat() { }
    void initDouble() { }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        basefft(false, realIn, 0, m_c, m_d);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) realOut[i] = m_c[i];
        if (imagOut) {
            for (int i = 0; i <= hs; ++i) imagOut[i] = m_d[i];
        }
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        basefft(false, realIn, 0, m_c, m_d);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) complexOut[i*2] = m_c[i];
        for (int i = 0; i <= hs; ++i) complexOut[i*2+1] = m_d[i];
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        basefft(false, realIn, 0, m_c, m_d);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_c[i] * m_c[i] + m_d[i] * m_d[i]);
            phaseOut[i] = atan2(m_d[i], m_c[i]) ;
        }
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        basefft(false, realIn, 0, m_c, m_d);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_c[i] * m_c[i] + m_d[i] * m_d[i]);
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        for (int i = 0; i < size_; ++i) m_a[i] = realIn[i];
        basefft(false, m_a, 0, m_c, m_d);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) realOut[i] = m_c[i];
        if (imagOut) {
            for (int i = 0; i <= hs; ++i) imagOut[i] = m_d[i];
        }
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        for (int i = 0; i < size_; ++i) m_a[i] = realIn[i];
        basefft(false, m_a, 0, m_c, m_d);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) complexOut[i*2] = m_c[i];
        for (int i = 0; i <= hs; ++i) complexOut[i*2+1] = m_d[i];
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        for (int i = 0; i < size_; ++i) m_a[i] = realIn[i];
        basefft(false, m_a, 0, m_c, m_d);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_c[i] * m_c[i] + m_d[i] * m_d[i]);
            phaseOut[i] = atan2(m_d[i], m_c[i]) ;
        }
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        for (int i = 0; i < size_; ++i) m_a[i] = realIn[i];
        basefft(false, m_a, 0, m_c, m_d);
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_c[i] * m_c[i] + m_d[i] * m_d[i]);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            double real = realIn[i];
            double imag = imagIn[i];
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[size_-i] = real;
                m_b[size_-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, realOut, m_d);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            double real = complexIn[i*2];
            double imag = complexIn[i*2+1];
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[size_-i] = real;
                m_b[size_-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, realOut, m_d);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            double real = magIn[i] * cos(phaseIn[i]);
            double imag = magIn[i] * sin(phaseIn[i]);
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[size_-i] = real;
                m_b[size_-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, realOut, m_d);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            double real = log(magIn[i] + 0.000001);
            m_a[i] = real;
            m_b[i] = 0.0;
            if (i > 0) {
                m_a[size_-i] = real;
                m_b[size_-i] = 0.0;
            }
        }
        basefft(true, m_a, m_b, cepOut, m_d);
    }

    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            float real = realIn[i];
            float imag = imagIn[i];
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[size_-i] = real;
                m_b[size_-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, m_c, m_d);
        for (int i = 0; i < size_; ++i) realOut[i] = m_c[i];
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            float real = complexIn[i*2];
            float imag = complexIn[i*2+1];
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[size_-i] = real;
                m_b[size_-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, m_c, m_d);
        for (int i = 0; i < size_; ++i) realOut[i] = m_c[i];
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            float real = magIn[i] * cosf(phaseIn[i]);
            float imag = magIn[i] * sinf(phaseIn[i]);
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[size_-i] = real;
                m_b[size_-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, m_c, m_d);
        for (int i = 0; i < size_; ++i) realOut[i] = m_c[i];
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        const int hs = size_/2;
        for (int i = 0; i <= hs; ++i) {
            float real = logf(magIn[i] + 0.000001);
            m_a[i] = real;
            m_b[i] = 0.0;
            if (i > 0) {
                m_a[size_-i] = real;
                m_b[size_-i] = 0.0;
            }
        }
        basefft(true, m_a, m_b, m_c, m_d);
        for (int i = 0; i < size_; ++i) cepOut[i] = m_c[i];
    }

private:
    const int size_;
    int *m_table;
    double *m_a;
    double *m_b;
    double *m_c;
    double *m_d;
    void basefft(bool inverse, const double *R__ ri, const double *R__ ii, double *R__ ro, double *R__ io);
};

void
D_Cross::basefft(bool inverse, const double *R__ ri, const double *R__ ii, double *R__ ro, double *R__ io)
{
    if (!ri || !ro || !io) return;

    int i, j, k, m;
    int blockSize, blockEnd;

    double tr, ti;

    double angle = 2.0 * M_PI;
    if (inverse) angle = -angle;

    const int n = size_;

    if (ii) {
	for (i = 0; i < n; ++i) {
	    ro[m_table[i]] = ri[i];
        }
	for (i = 0; i < n; ++i) {
	    io[m_table[i]] = ii[i];
	}
    } else {
	for (i = 0; i < n; ++i) {
	    ro[m_table[i]] = ri[i];
        }
	for (i = 0; i < n; ++i) {
	    io[m_table[i]] = 0.0;
	}
    }

    blockEnd = 1;

    for (blockSize = 2; blockSize <= n; blockSize <<= 1) {

	double delta = angle / (double)blockSize;
	double sm2 = -sin(-2 * delta);
	double sm1 = -sin(-delta);
	double cm2 = cos(-2 * delta);
	double cm1 = cos(-delta);
	double w = 2 * cm1;
	double ar[3], ai[3];

	for (i = 0; i < n; i += blockSize) {

	    ar[2] = cm2;
	    ar[1] = cm1;

	    ai[2] = sm2;
	    ai[1] = sm1;

	    for (j = i, m = 0; m < blockEnd; j++, m++) {

		ar[0] = w * ar[1] - ar[2];
		ar[2] = ar[1];
		ar[1] = ar[0];

		ai[0] = w * ai[1] - ai[2];
		ai[2] = ai[1];
		ai[1] = ai[0];

		k = j + blockEnd;
		tr = ar[0] * ro[k] - ai[0] * io[k];
		ti = ar[0] * io[k] + ai[0] * ro[k];

		ro[k] = ro[j] - tr;
		io[k] = io[j] - ti;

		ro[j] += tr;
		io[j] += ti;
	    }
	}

	blockEnd = blockSize;
    }

/* fftw doesn't rescale, so nor will we

    if (inverse) {

	double denom = (double)n;

	for (i = 0; i < n; i++) {
	    ro[i] /= denom;
	    io[i] /= denom;
	}
    }
*/
}

#endif /* USE_BUILTIN_FFT */

} /* end namespace FFTs */

std::string
FFT::implementation_;

std::set<std::string>
FFT::getImplementations()
{
    std::set<std::string> impls;
#ifdef HAVE_IPP
    impls.insert("ipp");
#endif
#ifdef HAVE_FFTW3
    impls.insert("fftw");
#endif
#ifdef USE_KISSFFT
    impls.insert("kissfft");
#endif
#ifdef HAVE_VDSP
    impls.insert("vdsp");
#endif
#ifdef HAVE_MEDIALIB
    impls.insert("medialib");
#endif
#ifdef HAVE_OPENMAX
    impls.insert("openmax");
#endif
#ifdef HAVE_SFFT
    impls.insert("sfft");
#endif
#ifdef USE_BUILTIN_FFT
    impls.insert("cross");
#endif
    return impls;
}

void
FFT::pickDefaultImplementation()
{
    if (implementation_ != "") return;

    std::set<std::string> impls = getImplementations();

    std::string best = "cross";
    if (impls.find("kissfft") != impls.end()) best = "kissfft";
    if (impls.find("medialib") != impls.end()) best = "medialib";
    if (impls.find("openmax") != impls.end()) best = "openmax";
    if (impls.find("sfft") != impls.end()) best = "sfft";
    if (impls.find("fftw") != impls.end()) best = "fftw";
    if (impls.find("vdsp") != impls.end()) best = "vdsp";
    if (impls.find("ipp") != impls.end()) best = "ipp";
    
    implementation_ = best;
}

std::string
FFT::getDefaultImplementation()
{
    return implementation_;
}

void
FFT::setDefaultImplementation(std::string i)
{
    implementation_ = i;
}

FFT::FFT(int size, int debugLevel) :
    d(0)
{
    if ((size < 2) ||
        (size & (size-1))) {
        std::cerr << "FFT::FFT(" << size << "): power-of-two sizes only supported, minimum size 2" << std::endl;
#ifndef NO_EXCEPTIONS
        throw InvalidSize;
#else
        abort();
#endif
    }

    if (implementation_ == "") pickDefaultImplementation();
    std::string impl = implementation_;

    if (debugLevel > 0) {
        std::cerr << "FFT::FFT(" << size << "): using implementation: "
                  << impl << std::endl;
    }

    if (impl == "ipp") {
#ifdef HAVE_IPP
        d = new FFTs::D_IPP(size);
#endif
    } else if (impl == "fftw") {
#ifdef HAVE_FFTW3
        d = new FFTs::D_FFTW(size);
#endif
    } else if (impl == "kissfft") {        
#ifdef USE_KISSFFT
        d = new FFTs::D_KISSFFT(size);
#endif
    } else if (impl == "vdsp") {
#ifdef HAVE_VDSP
        d = new FFTs::D_VDSP(size);
#endif
    } else if (impl == "medialib") {
#ifdef HAVE_MEDIALIB
        d = new FFTs::D_MEDIALIB(size);
#endif
    } else if (impl == "openmax") {
#ifdef HAVE_OPENMAX
        d = new FFTs::D_OPENMAX(size);
#endif
    } else if (impl == "sfft") {
#ifdef HAVE_SFFT
        d = new FFTs::D_SFFT(size);
#endif
    } else if (impl == "cross") {
#ifdef USE_BUILTIN_FFT
        d = new FFTs::D_Cross(size);
#endif
    }

    if (!d) {
        std::cerr << "FFT::FFT(" << size << "): ERROR: implementation "
                  << impl << " is not compiled in" << std::endl;
#ifndef NO_EXCEPTIONS
        throw InvalidImplementation;
#else
        abort();
#endif
    }
}

FFT::~FFT()
{
    delete d;
}

#ifndef NO_EXCEPTIONS
#define CHECK_NOT_NULL(x) \
    if (!(x)) { \
        std::cerr << "FFT: ERROR: Null argument " #x << std::endl;  \
        throw NullArgument; \
    }
#else
#define CHECK_NOT_NULL(x) \
    if (!(x)) { \
        std::cerr << "FFT: ERROR: Null argument " #x << std::endl;  \
        std::cerr << "FFT: Would be throwing NullArgument here, if exceptions were not disabled" << std::endl;  \
        return; \
    }
#endif

void
FFT::forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(realOut);
    CHECK_NOT_NULL(imagOut);
    d->forward(realIn, realOut, imagOut);
}

void
FFT::forwardInterleaved(const double *R__ realIn, double *R__ complexOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(complexOut);
    d->forwardInterleaved(realIn, complexOut);
}

void
FFT::forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(magOut);
    CHECK_NOT_NULL(phaseOut);
    d->forwardPolar(realIn, magOut, phaseOut);
}

void
FFT::forwardMagnitude(const double *R__ realIn, double *R__ magOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(magOut);
    d->forwardMagnitude(realIn, magOut);
}

void
FFT::forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(realOut);
    CHECK_NOT_NULL(imagOut);
    d->forward(realIn, realOut, imagOut);
}

void
FFT::forwardInterleaved(const float *R__ realIn, float *R__ complexOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(complexOut);
    d->forwardInterleaved(realIn, complexOut);
}

void
FFT::forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(magOut);
    CHECK_NOT_NULL(phaseOut);
    d->forwardPolar(realIn, magOut, phaseOut);
}

void
FFT::forwardMagnitude(const float *R__ realIn, float *R__ magOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(magOut);
    d->forwardMagnitude(realIn, magOut);
}

void
FFT::inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(imagIn);
    CHECK_NOT_NULL(realOut);
    d->inverse(realIn, imagIn, realOut);
}

void
FFT::inverseInterleaved(const double *R__ complexIn, double *R__ realOut)
{
    CHECK_NOT_NULL(complexIn);
    CHECK_NOT_NULL(realOut);
    d->inverseInterleaved(complexIn, realOut);
}

void
FFT::inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut)
{
    CHECK_NOT_NULL(magIn);
    CHECK_NOT_NULL(phaseIn);
    CHECK_NOT_NULL(realOut);
    d->inversePolar(magIn, phaseIn, realOut);
}

void
FFT::inverseCepstral(const double *R__ magIn, double *R__ cepOut)
{
    CHECK_NOT_NULL(magIn);
    CHECK_NOT_NULL(cepOut);
    d->inverseCepstral(magIn, cepOut);
}

void
FFT::inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(imagIn);
    CHECK_NOT_NULL(realOut);
    d->inverse(realIn, imagIn, realOut);
}

void
FFT::inverseInterleaved(const float *R__ complexIn, float *R__ realOut)
{
    CHECK_NOT_NULL(complexIn);
    CHECK_NOT_NULL(realOut);
    d->inverseInterleaved(complexIn, realOut);
}

void
FFT::inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut)
{
    CHECK_NOT_NULL(magIn);
    CHECK_NOT_NULL(phaseIn);
    CHECK_NOT_NULL(realOut);
    d->inversePolar(magIn, phaseIn, realOut);
}

void
FFT::inverseCepstral(const float *R__ magIn, float *R__ cepOut)
{
    CHECK_NOT_NULL(magIn);
    CHECK_NOT_NULL(cepOut);
    d->inverseCepstral(magIn, cepOut);
}

void
FFT::initFloat() 
{
    d->initFloat();
}

void
FFT::initDouble() 
{
    d->initDouble();
}

FFT::Precisions
FFT::getSupportedPrecisions() const
{
    return d->getSupportedPrecisions();
}

}
