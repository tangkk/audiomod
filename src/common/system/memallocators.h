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

#include "vectorization.h"

#include <new> // for std::bad_alloc
#include <stdlib.h>

#define MALLOC_IS_ALIGNED // temp fixed def
#define NO_EXCEPTIONS // temp fixed def

#ifndef HAVE_POSIX_MEMALIGN
#ifndef _WIN32
#ifndef __APPLE__
#ifndef LACK_POSIX_MEMALIGN
#define HAVE_POSIX_MEMALIGN
#endif
#endif
#endif
#endif

#ifdef HAVE_POSIX_MEMALIGN
#include <sys/mman.h>
#endif

#ifdef LACK_BAD_ALLOC
namespace std { struct bad_alloc { }; }
#endif

namespace audiomod {

template <typename T>
T *mem_allocate(size_t cnt)
{
    void *ptr = 0;

    // 32-byte alignment is required for at least OpenMAX
    // static const int alignment = 32;
    
#ifdef USE_OWN_ALIGNED_MALLOC
    static const int alignment = 32;
    // Alignment must be a power of two, bigger than the pointer
    // size. Stuff the actual malloc'd pointer in just before the
    // returned value.  This is the least desirable way to do this --
    // the other options below are all better
    size_t allocd = cnt * sizeof(T) + alignment;
    void *buf = malloc(allocd);
    if (buf) {
        char *adj = (char *)buf;
        while ((unsigned long long)adj & (alignment-1)) --adj;
        ptr = ((char *)adj) + alignment;
        ((void **)ptr)[-1] = buf;
    }
#else /* !USE_OWN_ALIGNED_MALLOC */
#ifdef HAVE_POSIX_MEMALIGN
    static const int alignment = 32;
    if (posix_memalign(&ptr, alignment, cnt * sizeof(T))) {
        ptr = malloc(cnt * sizeof(T));
    }
#else /* !HAVE_POSIX_MEMALIGN */
#ifdef __MSVC__
    static const int alignment = 32;
    ptr = _aligned_malloc(cnt * sizeof(T), alignment);
#else /* !__MSVC__ */
#ifndef MALLOC_IS_ALIGNED
#warning "No aligned malloc available or defined"
#endif
    // Note that malloc always aligns to 16 byte boundaries on OS/X
    ptr = malloc(cnt * sizeof(T));
#endif /* !__MSVC__ */
#endif /* !HAVE_POSIX_MEMALIGN */
#endif /* !USE_OWN_ALIGNED_MALLOC */
    if (!ptr) {
#ifndef NO_EXCEPTIONS
        throw(std::bad_alloc());
#else
        abort();
#endif
    }
    return (T *)ptr;
}

#ifdef HAVE_IPP

template <>
float *mem_allocate(size_t cnt)
{
    float *ptr = ippsMalloc_32f(cnt);
    if (!ptr) throw (std::bad_alloc());
    return ptr;
}

template <>
double *mem_allocate(size_t cnt)
{
    double *ptr = ippsMalloc_64f(cnt);
    if (!ptr) throw (std::bad_alloc());
    return ptr;
}

#endif
	
template <typename T>
T *mem_allocate_and_zero(size_t cnt)
{
    T *ptr = mem_allocate<T>(cnt);
    vector_zeros(ptr, cnt);
    return ptr;
}

template <typename T>
void mem_deallocate(T *ptr)
{
#ifdef USE_OWN_ALIGNED_MALLOC
    if (ptr) free(((void **)ptr)[-1]);
#else /* !USE_OWN_ALIGNED_MALLOC */
#ifdef __MSVC__
    if (ptr) _aligned_free((void *)ptr);
#else /* !__MSVC__ */
    if (ptr) free((void *)ptr);
#endif /* !__MSVC__ */
#endif /* !USE_OWN_ALIGNED_MALLOC */
}

#ifdef HAVE_IPP

template <>
void mem_deallocate(float *ptr)
{
    if (ptr) ippsFree((void *)ptr);
}

template <>
void mem_deallocate(double *ptr)
{
    if (ptr) ippsFree((void *)ptr);
}

#endif

/// Reallocate preserving contents but leaving additional memory uninitialised	
template <typename T>
T *mem_reallocate(T *ptr, size_t oldcnt, size_t cnt)
{
    T *newptr = mem_allocate<T>(cnt);
    if (oldcnt && ptr) {
        vector_copy(newptr, ptr, oldcnt < cnt ? oldcnt : cnt);
    }
    if (ptr) mem_deallocate<T>(ptr);
    return newptr;
}

/// Reallocate, zeroing all contents
template <typename T>
T *mem_reallocate_zero(T *ptr, size_t oldcnt, size_t cnt)
{
    ptr = mem_reallocate(ptr, oldcnt, cnt);
    vector_zeros(ptr, cnt);
    return ptr;
}
	
/// Reallocate preserving contents and zeroing any additional memory	
template <typename T>
T *mem_reallocate_zeroext(T *ptr, size_t oldcnt, size_t cnt)
{
    ptr = mem_reallocate(ptr, oldcnt, cnt);
    if (cnt > oldcnt) vector_zeros(ptr + oldcnt, cnt - oldcnt);
    return ptr;
}

template <typename T>
T **mem_allocate_channels(size_t chans, size_t cnt)
{
    T **ptr = mem_allocate<T *>(chans);
    for (size_t c = 0; c < chans; ++c) {
        ptr[c] = mem_allocate<T>(cnt);
    }
    return ptr;
}
	
template <typename T>
T **mem_allocate_and_zero_channels(size_t chans, size_t cnt)
{
    T **ptr = mem_allocate<T *>(chans);
    for (size_t c = 0; c < chans; ++c) {
        ptr[c] = mem_allocate_and_zero<T>(cnt);
    }
    return ptr;
}

template <typename T>
void mem_deallocate_channels(T **ptr, size_t chans)
{
    if (!ptr) return;
    for (size_t c = 0; c < chans; ++c) {
        mem_deallocate<T>(ptr[c]);
    }
    mem_deallocate<T *>(ptr);
}
	
template <typename T>
T **mem_reallocate_channels(T **ptr,
                        size_t oldchans, size_t oldcnt,
                        size_t chans, size_t cnt)
{
    T **newptr = mem_allocate_channels<T>(chans, cnt);
    if (oldcnt && ptr) {
        vector_copy_channels(newptr, ptr, chans, oldcnt < cnt ? oldcnt : cnt);
    } 
    if (ptr) mem_deallocate_channels<T>(ptr, chans);
    return newptr;
}
	
template <typename T>
T **mem_reallocate_and_zeroext_channels(T **ptr,
                                        size_t oldchans, size_t oldcnt,
                                        size_t chans, size_t cnt)
{
    T **newptr = mem_allocate_and_zero_channels<T>(chans, cnt);
    if (oldcnt && ptr) {
        vector_copy_channels(newptr, ptr, chans, oldcnt < cnt ? oldcnt : cnt);
    } 
    if (ptr) mem_deallocate_channels<T>(ptr, chans);
    return newptr;
}

}

