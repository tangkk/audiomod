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

#pragma once

#ifdef HAVE_IPP
#ifndef _MSC_VER
#include <inttypes.h>
#endif
#include <ippversion.h>
#include <ipps.h>
#if (IPP_VERSION_MAJOR <= 7)
// Deprecated in v8, removed in v9
#include <ippac.h>
#endif
#endif

#ifdef HAVE_VDSP
#include <Accelerate/Accelerate.h>
#endif

#include <cstring>
#include "sys.h"

#if defined USE_POMMIER_MATHFUN
#if defined __ARMEL__ || defined __aarch64__
#include "../pommier/neon_mathfun.h"
#else
#include "../pommier/sse_mathfun.h"
#endif
#endif

namespace audiomod {

// Note that all functions with a "target" vector have their arguments
// in the same order as memcpy and friends, i.e. target vector first.
// This is the reverse order from the IPP functions.

// The ideal here is to write the basic loops in such a way as to be
// auto-vectorizable by a sensible compiler (definitely gcc-4.3 on
// Linux, ideally also gcc-4.0 on OS/X).

template<typename T>
inline void vector_zeros(T *const R__ ptr,
                   const int cnt)
{
    const T value = T(0);
    for (int i = 0; i < cnt; ++i) {
        ptr[i] = value;
    }
}

#if defined HAVE_IPP
template<> 
inline void vector_zeros(float *const R__ ptr, 
                   const int cnt)
{
    ippsZero_32f(ptr, cnt);
}
template<> 
inline void vector_zeros(double *const R__ ptr,
                   const int cnt)
{
    ippsZero_64f(ptr, cnt);
}
#elif defined HAVE_VDSP
template<> 
inline void vector_zeros(float *const R__ ptr, 
                   const int cnt)
{
    vDSP_vclr(ptr, 1, cnt);
}
template<> 
inline void vector_zeros(double *const R__ ptr,
                   const int cnt)
{
    vDSP_vclrD(ptr, 1, cnt);
}
#endif

template<typename T>
inline void vector_zeros_channels(T *const R__ *const R__ ptr,
                            const int chans,
                            const int cnt)
{
    for (int c = 0; c < chans; ++c) {
        vector_zeros(ptr[c], cnt);
    }
}

template<typename T>
inline void vector_set(T *const R__ ptr,
                  const T value,
                  const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        ptr[i] = value;
    }
}

template<typename T>
inline void vector_copy(T *const R__ dst,
                   const T *const R__ src,
                   const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] = src[i];
    }
}

#if defined HAVE_IPP
template<>
inline void vector_copy(float *const R__ dst,
                   const float *const R__ src,
                   const int cnt)
{
    ippsCopy_32f(src, dst, cnt);
}
template<>
inline void vector_copy(double *const R__ dst,
                   const double *const R__ src,
                   const int cnt)
{
    ippsCopy_64f(src, dst, cnt);
}
#endif

template<typename T>
inline void vector_copy_channels(T *const R__ *const R__ dst,
                            const T *const R__ *const R__ src,
                            const int chans,
                            const int cnt)
{
    for (int c = 0; c < chans; ++c) {
        vector_copy(dst[c], src[c], cnt);
    }
}

// src and dst alias by definition, so not restricted
template<typename T>
inline void vector_move(T *const dst,
                   const T *const src,
                   const int cnt)
{
    memmove(dst, src, cnt * sizeof(T));
}

#if defined HAVE_IPP
template<>
inline void vector_move(float *const dst,
                   const float *const src,
                   const int cnt)
{
    ippsMove_32f(src, dst, cnt);
}
template<>
inline void vector_move(double *const dst,
                   const double *const src,
                   const int cnt)
{
    ippsMove_64f(src, dst, cnt);
}
#endif

template<typename T, typename U>
inline void vector_assign(U *const R__ dst,
                      const T *const R__ src,
                      const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] = U(src[i]);
    }
}

template<>
inline void vector_assign(float *const R__ dst,
                      const float *const R__ src,
                      const int cnt)
{
    vector_copy(dst, src, cnt);
}
template<>
inline void vector_assign(double *const R__ dst,
                      const double *const R__ src,
                      const int cnt)
{
    vector_copy(dst, src, cnt);
}

#if defined HAVE_IPP
template<>
inline void vector_assign(double *const R__ dst,
                      const float *const R__ src,
                      const int cnt)
{
    ippsConvert_32f64f(src, dst, cnt);
}
template<>
inline void vector_assign(float *const R__ dst,
                      const double *const R__ src,
                      const int cnt)
{
    ippsConvert_64f32f(src, dst, cnt);
}
#elif defined HAVE_VDSP
template<>
inline void vector_assign(double *const R__ dst,
                      const float *const R__ src,
                      const int cnt)
{
    vDSP_vspdp((float *)src, 1, dst, 1, cnt);
}
template<>
inline void vector_assign(float *const R__ dst,
                      const double *const R__ src,
                      const int cnt)
{
    vDSP_vdpsp((double *)src, 1, dst, 1, cnt);
}
#endif

template<typename T, typename U>
inline void vector_assign_channels(U *const R__ *const R__ dst,
                               const T *const R__ *const R__ src,
                               const int chans,
                               const int cnt)
{
    for (int c = 0; c < chans; ++c) {
        vector_assign(dst[c], src[c], cnt);
    }
}

template<typename T>
inline void vector_add(T *const R__ dst,
                  const T *const R__ src,
                  const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] += src[i];
    }
}

template<typename T>
inline void vector_add(T *const R__ dst,
                  const T value,
                  const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] += value;
    }
}

#if defined HAVE_IPP
template<>
inline void vector_add(float *const R__ dst,
                  const float *const R__ src,
                  const int cnt)
{
    ippsAdd_32f_I(src, dst, cnt);
}    
inline void vector_add(double *const R__ dst,
                  const double *const R__ src,
                  const int cnt)
{
    ippsAdd_64f_I(src, dst, cnt);
}    
#endif

template<typename T>
inline void vector_add_channels(T *const R__ *const R__ dst,
                           const T *const R__ *const R__ src,
                           const int chans, const int cnt)
{
    for (int c = 0; c < chans; ++c) {
        vector_add(dst[c], src[c], cnt);
    }
}

template<typename T, typename G>
inline void vector_add_with_gain(T *const R__ dst,
                            const T *const R__ src,
                            const G gain,
                            const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] += src[i] * gain;
    }
}

template<typename T, typename G>
inline void vector_add_channels_with_gain(T *const R__ *const R__ dst,
                                     const T *const R__ *const R__ src,
                                     const G gain,
                                     const int chans,
                                     const int cnt)
{
    for (int c = 0; c < chans; ++c) {
        vector_add_with_gain(dst[c], src[c], gain, cnt);
    }
}

template<typename T>
inline void vector_subtract(T *const R__ dst,
                       const T *const R__ src,
                       const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] -= src[i];
    }
}

#if defined HAVE_IPP
template<>
inline void vector_subtract(float *const R__ dst,
                       const float *const R__ src,
                       const int cnt)
{
    ippsSub_32f_I(src, dst, cnt);
}    
inline void vector_subtract(double *const R__ dst,
                       const double *const R__ src,
                       const int cnt)
{
    ippsSub_64f_I(src, dst, cnt);
}    
#endif

template<typename T, typename G>
inline void vector_gain(T *const R__ dst,
                    const G gain,
                    const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] *= gain;
    }
}

#if defined HAVE_IPP 
template<>
inline void vector_gain(float *const R__ dst,
                    const float gain,
                    const int cnt)
{
    ippsMulC_32f_I(gain, dst, cnt);
}
template<>
inline void vector_gain(double *const R__ dst,
                    const double gain,
                    const int cnt)
{
    ippsMulC_64f_I(gain, dst, cnt);
}
#endif

template<typename T>
inline void vector_multiply(T *const R__ dst,
                       const T *const R__ src,
                       const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] *= src[i];
    }
}

#if defined HAVE_IPP 
template<>
inline void vector_multiply(float *const R__ dst,
                       const float *const R__ src,
                       const int cnt)
{
    ippsMul_32f_I(src, dst, cnt);
}
template<>
inline void vector_multiply(double *const R__ dst,
                       const double *const R__ src,
                       const int cnt)
{
    ippsMul_64f_I(src, dst, cnt);
}
#endif

template<typename T>
inline void vector_multiply(T *const R__ dst,
                       const T *const R__ src1,
                       const T *const R__ src2,
                       const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] = src1[i] * src2[i];
    }
}

template<typename T>
inline void vector_divide(T *const R__ dst,
                     const T *const R__ src,
                     const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] /= src[i];
    }
}

#if defined HAVE_IPP 
template<>
inline void vector_divide(float *const R__ dst,
                     const float *const R__ src,
                     const int cnt)
{
    ippsDiv_32f_I(src, dst, cnt);
}
template<>
inline void vector_divide(double *const R__ dst,
                     const double *const R__ src,
                     const int cnt)
{
    ippsDiv_64f_I(src, dst, cnt);
}
#endif

#if defined HAVE_IPP 
template<>
inline void vector_multiply(float *const R__ dst,
                       const float *const R__ src1,
                       const float *const R__ src2,
                       const int cnt)
{
    ippsMul_32f(src1, src2, dst, cnt);
}    
template<>
inline void vector_multiply(double *const R__ dst,
                       const double *const R__ src1,
                       const double *const R__ src2,
                       const int cnt)
{
    ippsMul_64f(src1, src2, dst, cnt);
}
#endif

template<typename T>
inline void vector_multiply_and_add(T *const R__ dst,
                               const T *const R__ src1,
                               const T *const R__ src2,
                               const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] += src1[i] * src2[i];
    }
}

#if defined HAVE_IPP
template<>
inline void vector_multiply_and_add(float *const R__ dst,
                               const float *const R__ src1,
                               const float *const R__ src2,
                               const int cnt)
{
    ippsAddProduct_32f(src1, src2, dst, cnt);
}
template<>
inline void vector_multiply_and_add(double *const R__ dst,
                               const double *const R__ src1,
                               const double *const R__ src2,
                               const int cnt)
{
    ippsAddProduct_64f(src1, src2, dst, cnt);
}
#endif

template<typename T>
inline T vector_sum(const T *const R__ src,
               const int cnt)
{
    T result = T();
    for (int i = 0; i < cnt; ++i) {
        result += src[i];
    }
    return result;
}

template<typename T>
inline void vector_log(T *const R__ dst,
                  const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] = log(dst[i]);
    }
}

#if defined HAVE_IPP
template<>
inline void vector_log(float *const R__ dst,
                  const int cnt)
{
    ippsLn_32f_I(dst, cnt);
}
template<>
inline void vector_log(double *const R__ dst,
                  const int cnt)
{
    ippsLn_64f_I(dst, cnt);
}
#elif defined HAVE_VDSP
// no in-place vForce functions for these -- can we use the
// out-of-place functions with equal input and output vectors? can we
// use an out-of-place one with temporary buffer and still be faster
// than doing it any other way?
template<>
inline void vector_log(float *const R__ dst,
                  const int cnt)
{
    float tmp[cnt];
    vvlogf(tmp, dst, &cnt);
    vector_copy(dst, tmp, cnt);
}
template<>
inline void vector_log(double *const R__ dst,
                  const int cnt)
{
    double tmp[cnt];
    vvlog(tmp, dst, &cnt);
    vector_copy(dst, tmp, cnt);
}
#endif

template<typename T>
inline void vector_exp(T *const R__ dst,
                  const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] = exp(dst[i]);
    }
}

#if defined HAVE_IPP
template<>
inline void vector_exp(float *const R__ dst,
                  const int cnt)
{
    ippsExp_32f_I(dst, cnt);
}
template<>
inline void vector_exp(double *const R__ dst,
                  const int cnt)
{
    ippsExp_64f_I(dst, cnt);
}
#elif defined HAVE_VDSP
// no in-place vForce functions for these -- can we use the
// out-of-place functions with equal input and output vectors? can we
// use an out-of-place one with temporary buffer and still be faster
// than doing it any other way?
template<>
inline void vector_exp(float *const R__ dst,
                  const int cnt)
{
    float tmp[cnt];
    vvexpf(tmp, dst, &cnt);
    vector_copy(dst, tmp, cnt);
}
template<>
inline void vector_exp(double *const R__ dst,
                  const int cnt)
{
    double tmp[cnt];
    vvexp(tmp, dst, &cnt);
    vector_copy(dst, tmp, cnt);
}
#endif

template<typename T>
inline void vector_sqrt(T *const R__ dst,
                   const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] = sqrt(dst[i]);
    }
}

#if defined HAVE_IPP
template<>
inline void vector_sqrt(float *const R__ dst,
                   const int cnt)
{
    ippsSqrt_32f_I(dst, cnt);
}
template<>
inline void vector_sqrt(double *const R__ dst,
                   const int cnt)
{
    ippsSqrt_64f_I(dst, cnt);
}
#elif defined HAVE_VDSP
// no in-place vForce functions for these -- can we use the
// out-of-place functions with equal input and output vectors? can we
// use an out-of-place one with temporary buffer and still be faster
// than doing it any other way?
template<>
inline void vector_sqrt(float *const R__ dst,
                   const int cnt)
{
    float tmp[cnt];
    vvsqrtf(tmp, dst, &cnt);
    vector_copy(dst, tmp, cnt);
}
template<>
inline void vector_sqrt(double *const R__ dst,
                   const int cnt)
{
    double tmp[cnt];
    vvsqrt(tmp, dst, &cnt);
    vector_copy(dst, tmp, cnt);
}
#endif

template<typename T>
inline void vector_square(T *const R__ dst,
                   const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] = dst[i] * dst[i];
    }
}

#if defined HAVE_IPP
template<>
inline void vector_square(float *const R__ dst,
                   const int cnt)
{
    ippsSqr_32f_I(dst, cnt);
}
template<>
inline void vector_square(double *const R__ dst,
                   const int cnt)
{
    ippsSqr_64f_I(dst, cnt);
}
#endif

template<typename T>
inline void vector_abs(T *const R__ dst,
                  const int cnt)
{
    for (int i = 0; i < cnt; ++i) {
        dst[i] = fabs(dst[i]);
    }
}

#if defined HAVE_IPP
template<>
inline void vector_abs(float *const R__ dst,
                  const int cnt)
{
    ippsAbs_32f_I(dst, cnt);
}
template<>
inline void vector_abs(double *const R__ dst,
                  const int cnt)
{
    ippsAbs_64f_I(dst, cnt);
}
#elif defined HAVE_VDSP
template<>
inline void vector_abs(float *const R__ dst,
                  const int cnt)
{
    float tmp[cnt];
#if (defined(MACOSX_DEPLOYMENT_TARGET) && MACOSX_DEPLOYMENT_TARGET <= 1070 && MAC_OS_X_VERSION_MIN_REQUIRED <= 1070)
    vvfabf(tmp, dst, &cnt);
#else
    vvfabsf(tmp, dst, &cnt);
#endif
    vector_copy(dst, tmp, cnt);
}
#endif

template<typename T>
inline void vector_interleave(T *const R__ dst,
                         const T *const R__ *const R__ src,
                         const int chans, 
                         const int cnt)
{
    int idx = 0;
    switch (chans) {
    case 2:
        // common case, may be vectorized by compiler if hardcoded
        for (int i = 0; i < cnt; ++i) {
            for (int j = 0; j < 2; ++j) {
                dst[idx++] = src[j][i];
            }
        }
        return;
    case 1:
        vector_copy(dst, src[0], cnt);
        return;
    default:
        for (int i = 0; i < cnt; ++i) {
            for (int j = 0; j < chans; ++j) {
                dst[idx++] = src[j][i];
            }
        }
    }
}

#if defined HAVE_IPP 
#if (IPP_VERSION_MAJOR <= 7)
// Deprecated in v8, removed in v9
template<>
inline void vector_interleave(float *const R__ dst,
                         const float *const R__ *const R__ src,
                         const int chans, 
                         const int cnt)
{
    ippsInterleave_32f((const Ipp32f **)src, chans, cnt, dst);
}
// IPP does not (currently?) provide double-precision interleave
#endif
#endif

template<typename T>
inline void vector_deinterleave(T *const R__ *const R__ dst,
                           const T *const R__ src,
                           const int chans, 
                           const int cnt)
{
    int idx = 0;
    switch (chans) {
    case 2:
        // common case, may be vectorized by compiler if hardcoded
        for (int i = 0; i < cnt; ++i) {
            for (int j = 0; j < 2; ++j) {
                dst[j][i] = src[idx++];
            }
        }
        return;
    case 1:
        vector_copy(dst[0], src, cnt);
        return;
    default:
        for (int i = 0; i < cnt; ++i) {
            for (int j = 0; j < chans; ++j) {
                dst[j][i] = src[idx++];
            }
        }
    }
}

#if defined HAVE_IPP
#if (IPP_VERSION_MAJOR <= 7)
// Deprecated in v8, removed in v9
template<>
inline void vector_deinterleave(float *const R__ *const R__ dst,
                           const float *const R__ src,
                           const int chans, 
                           const int cnt)
{
    ippsDeinterleave_32f((const Ipp32f *)src, chans, cnt, (Ipp32f **)dst);
}
// IPP does not (currently?) provide double-precision deinterleave
#endif
#endif


template<typename T>
inline void complex_phasor(T *real, T *imag, T phase)
{
    //!!! IPP contains ippsSinCos_xxx in ippvm.h -- these are
    //!!! fixed-accuracy, test and compare
#if defined HAVE_VDSP
    int one = 1;
    if (sizeof(T) == sizeof(float)) {
        vvsincosf((float *)imag, (float *)real, (const float *)&phase, &one);
    } else {
        vvsincos((double *)imag, (double *)real, (const double *)&phase, &one);
    }
#elif defined LACK_SINCOS
    if (sizeof(T) == sizeof(float)) {
        *real = cosf(phase);
        *imag = sinf(phase);
    } else {
        *real = cos(phase);
        *imag = sin(phase);
    }
#elif defined __GNUC__
    // if (sizeof(T) == sizeof(float)) {
    //     sincosf(phase, (float *)imag, (float *)real);
    // } else {
    //     sincos(phase, (double *)imag, (double *)real);
    // }
    if (sizeof(T) == sizeof(float)) {
        __sincosf(phase, (float *)imag, (float *)real);
    } else {
        __sincos(phase, (double *)imag, (double *)real);
    }
#else
    if (sizeof(T) == sizeof(float)) {
        *real = cosf(phase);
        *imag = sinf(phase);
    } else {
        *real = cos(phase);
        *imag = sin(phase);
    }
#endif
}

template<typename T>
inline void complex_magphase(T *mag, T *phase, T real, T imag)
{
    *mag = sqrt(real * real + imag * imag);
    *phase = atan2(imag, real);
}

#ifdef USE_APPROXIMATE_ATAN2
// NB arguments in opposite order from usual for atan2f
extern float approximate_atan2f(float real, float imag);
template<>
inline void complex_magphase(float *mag, float *phase, float real, float imag)
{
    float atan = approximate_atan2f(real, imag);
    *phase = atan;
    *mag = sqrtf(real * real + imag * imag);
}
#else
template<>
inline void complex_magphase(float *mag, float *phase, float real, float imag)
{
    *mag = sqrtf(real * real + imag * imag);
    *phase = atan2f(imag, real);
}
#endif


template<typename S, typename T> // S source, T target
void vector_polar_to_cartesian(T *const R__ real,
                          T *const R__ imag,
                          const S *const R__ mag,
                          const S *const R__ phase,
                          const int count)
{
    for (int i = 0; i < count; ++i) {
        complex_phasor<T>(real + i, imag + i, phase[i]);
    }
    vector_multiply(real, mag, count);
    vector_multiply(imag, mag, count);
}

template<typename T>
void vector_polar_interleaved_to_cartesian_inplace(T *const R__ srcdst,
                                              const int count)
{
    T real, imag;
    for (int i = 0; i < count*2; i += 2) {
        complex_phasor(&real, &imag, srcdst[i+1]);
        real *= srcdst[i];
        imag *= srcdst[i];
        srcdst[i] = real;
        srcdst[i+1] = imag;
    }
}

template<typename S, typename T> // S source, T target
void vector_polar_to_cartesian_interleaved(T *const R__ dst,
                                      const S *const R__ mag,
                                      const S *const R__ phase,
                                      const int count)
{
    T real, imag;
    for (int i = 0; i < count; ++i) {
        complex_phasor<T>(&real, &imag, phase[i]);
        real *= mag[i];
        imag *= mag[i];
        dst[i*2] = real;
        dst[i*2+1] = imag;
    }
}    

#if defined USE_POMMIER_MATHFUN
void vector_polar_to_cartesian_pommier(float *const R__ real,
                                  float *const R__ imag,
                                  const float *const R__ mag,
                                  const float *const R__ phase,
                                  const int count);
void vector_polar_interleaved_to_cartesian_inplace_pommier(float *const R__ srcdst,
                                                      const int count);
void vector_polar_to_cartesian_interleaved_pommier(float *const R__ dst,
                                              const float *const R__ mag,
                                              const float *const R__ phase,
                                              const int count);

template<>
inline void vector_polar_to_cartesian(float *const R__ real,
                                 float *const R__ imag,
                                 const float *const R__ mag,
                                 const float *const R__ phase,
                                 const int count)
{
    vector_polar_to_cartesian_pommier(real, imag, mag, phase, count);
}

template<>
inline void vector_polar_interleaved_to_cartesian_inplace(float *const R__ srcdst,
                                                     const int count)
{
    vector_polar_interleaved_to_cartesian_inplace_pommier(srcdst, count);
}

template<>
inline void vector_polar_to_cartesian_interleaved(float *const R__ dst,
                                             const float *const R__ mag,
                                             const float *const R__ phase,
                                             const int count)
{
    vector_polar_to_cartesian_interleaved_pommier(dst, mag, phase, count);
}

#endif

template<typename S, typename T> // S source, T target
void vector_cartesian_to_polar(T *const R__ mag,
                          T *const R__ phase,
                          const S *const R__ real,
                          const S *const R__ imag,
                          const int count)
{
    for (int i = 0; i < count; ++i) {
        complex_magphase<T>(mag + i, phase + i, real[i], imag[i]);
    }
}

template<typename S, typename T> // S source, T target
void vector_cartesian_interleaved_to_polar(T *const R__ mag,
                                      T *const R__ phase,
                                      const S *const R__ src,
                                      const int count)
{
    for (int i = 0; i < count; ++i) {
        complex_magphase<T>(mag + i, phase + i, src[i*2], src[i*2+1]);
    }
}

#ifdef HAVE_VDSP
template<>
inline void vector_cartesian_to_polar(float *const R__ mag,
                                 float *const R__ phase,
                                 const float *const R__ real,
                                 const float *const R__ imag,
                                 const int count)
{
    DSPSplitComplex c;
    c.realp = const_cast<float *>(real);
    c.imagp = const_cast<float *>(imag);
    vDSP_zvmags(&c, 1, phase, 1, count); // using phase as a temporary dest
    vvsqrtf(mag, phase, &count); // using phase as the source
    vvatan2f(phase, imag, real, &count);
}
template<>
inline void vector_cartesian_to_polar(double *const R__ mag,
                                 double *const R__ phase,
                                 const double *const R__ real,
                                 const double *const R__ imag,
                                 const int count)
{
    // double precision, this is significantly faster than using vDSP_polar
    DSPDoubleSplitComplex c;
    c.realp = const_cast<double *>(real);
    c.imagp = const_cast<double *>(imag);
    vDSP_zvmagsD(&c, 1, phase, 1, count); // using phase as a temporary dest
    vvsqrt(mag, phase, &count); // using phase as the source
    vvatan2(phase, imag, real, &count);
}
#endif

template<typename T>
void vector_cartesian_to_polar_interleaved_inplace(T *const R__ srcdst,
                                              const int count)
{
    T mag, phase;
    for (int i = 0; i < count * 2; i += 2) {
        complex_magphase(&mag, &phase, srcdst[i], srcdst[i+1]);
        srcdst[i] = mag;
        srcdst[i+1] = phase;
    }
}


#ifdef USE_APPROXIMATE_ATAN2
float approximate_atan2f(float real, float imag)
{
    static const float pi = M_PI;
    static const float pi2 = M_PI / 2;

    float atan;

    if (real == 0.f) {

        if (imag > 0.0f) atan = pi2;
        else if (imag == 0.0f) atan = 0.0f;
        else atan = -pi2;

    } else {

        float z = imag/real;

        if (fabsf(z) < 1.f) {
            atan = z / (1.f + 0.28f * z * z);
            if (real < 0.f) {
                if (imag < 0.f) atan -= pi;
                else atan += pi;
            }
        } else {
            atan = pi2 - z / (z * z + 0.28f);
            if (imag < 0.f) atan -= pi;
        }
    }

    return atan;
}
#endif

#if defined USE_POMMIER_MATHFUN

#if defined __ARMEL__ || defined __aarch64__
typedef union {
  float f[4];
  int i[4];
  v4sf  v;
} V4SF;
#else
typedef ALIGN16_BEG union {
  float f[4];
  int i[4];
  v4sf  v;
} ALIGN16_END V4SF;
#endif

void
v_polar_to_cartesian_pommier(float *const R__ real,
                             float *const R__ imag,
                             const float *const R__ mag,
                             const float *const R__ phase,
                             const int count)
{
    int idx = 0, tidx = 0;
    int i = 0;

    for (int i = 0; i + 4 < count; i += 4) {

	V4SF fmag, fphase, fre, fim;

        for (int j = 0; j < 3; ++j) {
            fmag.f[j] = mag[idx];
            fphase.f[j] = phase[idx++];
        }

	sincos_ps(fphase.v, &fim.v, &fre.v);

        for (int j = 0; j < 3; ++j) {
            real[tidx] = fre.f[j] * fmag.f[j];
            imag[tidx++] = fim.f[j] * fmag.f[j];
        }
    }

    while (i < count) {
        float re, im;
        complex_phasor(&re, &im, phase[i]);
        real[tidx] = re * mag[i];
        imag[tidx++] = im * mag[i];
        ++i;
    }
}    

void
v_polar_interleaved_to_cartesian_inplace_pommier(float *const R__ srcdst,
                                                 const int count)
{
    int i;
    int idx = 0, tidx = 0;

    for (i = 0; i + 4 < count; i += 4) {

	V4SF fmag, fphase, fre, fim;

        for (int j = 0; j < 3; ++j) {
            fmag.f[j] = srcdst[idx++];
            fphase.f[j] = srcdst[idx++];
        }

	sincos_ps(fphase.v, &fim.v, &fre.v);

        for (int j = 0; j < 3; ++j) {
            srcdst[tidx++] = fre.f[j] * fmag.f[j];
            srcdst[tidx++] = fim.f[j] * fmag.f[j];
        }
    }

    while (i < count) {
        float real, imag;
        float mag = srcdst[idx++];
        float phase = srcdst[idx++];
        complex_phasor(&real, &imag, phase);
        srcdst[tidx++] = real * mag;
        srcdst[tidx++] = imag * mag;
        ++i;
    }
}    

void
v_polar_to_cartesian_interleaved_pommier(float *const R__ dst,
                                         const float *const R__ mag,
                                         const float *const R__ phase,
                                         const int count)
{
    int i;
    int idx = 0, tidx = 0;

    for (i = 0; i + 4 <= count; i += 4) {

	V4SF fmag, fphase, fre, fim;

        for (int j = 0; j < 3; ++j) {
            fmag.f[j] = mag[idx];
            fphase.f[j] = phase[idx];
            ++idx;
        }

	sincos_ps(fphase.v, &fim.v, &fre.v);

        for (int j = 0; j < 3; ++j) {
            dst[tidx++] = fre.f[j] * fmag.f[j];
            dst[tidx++] = fim.f[j] * fmag.f[j];
        }
    }

    while (i < count) {
        float real, imag;
        complex_phasor(&real, &imag, phase[i]);
        dst[tidx++] = real * mag[i];
        dst[tidx++] = imag * mag[i];
        ++i;
    }
}    

#endif


} // end of namespace audiomod
