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

#ifdef __MSVC__
#if _MSC_VER < 1800
#include "float_cast/float_cast.h"
#endif
#define R__ __restrict
#endif

#ifdef __clang__
#define R__ __restrict__
#else
#ifdef __GNUC__
#define R__ __restrict__
#endif
#endif

#ifndef R__
#define R__
#endif

#ifdef __MINGW32__
#include <malloc.h>
#else
#ifndef __MSVC__
#include <alloca.h>
#endif
#endif

#ifdef __MSVC__
#include <malloc.h>
#include <process.h>
#define alloca _alloca
#define getpid _getpid
#endif

#if defined(__MSVC__) && _MSC_VER < 1700
#define uint8_t unsigned __int8
#define uint16_t unsigned __int16
#define uint32_t unsigned __int32
#elif defined(__MSVC__)
#define ssize_t long
#include <stdint.h>
#else
#include <stdint.h>
#endif

#include <math.h>

namespace audiomod {

extern void system_specific_initialise();

enum ProcessStatus { ProcessRunning, ProcessNotRunning, UnknownProcessStatus };

#ifdef _WIN32

struct timeval { long tv_sec; long tv_usec; };
void gettimeofday(struct timeval *p, void *tz);

#endif

#ifdef __MSVC__

void usleep(unsigned long);

#endif

inline double mod(double x, double y) { return x - (y * floor(x / y)); }
inline float modf(float x, float y) { return x - (y * float(floor(x / y))); }

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline double princarg(double a) { return mod(a + M_PI, -2.0 * M_PI) + M_PI; }
inline float princargf(float a) { return modf(a + (float)M_PI, -2.f * (float)M_PI) + (float)M_PI; }

} // end namespace

// The following should be functions in the Common namespace, really

#ifdef _WIN32

#define MLOCK(a,b)   1
#define MUNLOCK(a,b) 1
#define MUNLOCK_SAMPLEBLOCK(a) 1

namespace audiomod {
extern void system_memorybarrier();
}
#define MBARRIER() audiomod::system_memorybarrier()

#define DLOPEN(a,b)  LoadLibrary((a).toStdWString().c_str())
#define DLSYM(a,b)   GetProcAddress((HINSTANCE)(a),(b))
#define DLCLOSE(a)   FreeLibrary((HINSTANCE)(a))
#define DLERROR()    ""

#else

#include <sys/mman.h>
#include <dlfcn.h>
#include <stdio.h>

#define MLOCK(a,b)   ::mlock((char *)(a),(b))
#define MUNLOCK(a,b) (::munlock((char *)(a),(b)) ? (::perror("munlock failed"), 0) : 0)
#define MUNLOCK_SAMPLEBLOCK(a) do { if (!(a).empty()) { const float &b = *(a).begin(); MUNLOCK(&b, (a).capacity() * sizeof(float)); } } while(0);

#ifdef __APPLE__
#include <libkern/OSAtomic.h>
// #define MBARRIER() OSMemoryBarrier()
#define MBARRIER() std::atomic_thread_fence(std::memory_order_relaxed)
#else
#if (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 1)
#define MBARRIER() __sync_synchronize()
#else
namespace audiomod {
extern void system_memorybarrier();
}
#define MBARRIER() ::audiomod::system_memorybarrier()
#endif
#endif

#define DLOPEN(a,b)  dlopen((a).toStdString().c_str(),(b))
#define DLSYM(a,b)   dlsym((a),(b))
#define DLCLOSE(a)   dlclose((a))
#define DLERROR()    dlerror()

#endif

#ifdef NO_THREADING
#undef MBARRIER
#define MBARRIER() 
#endif

