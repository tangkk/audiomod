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

#include "sys.h"

#ifdef _WIN32
#include <windows.h>
#include <fcntl.h>
#include <io.h>
#else /* !_WIN32 */
#include <signal.h>
#include <unistd.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#include <mach/mach.h>
#include <mach/mach_time.h>
#else /* !__APPLE__, !_WIN32 */
#include <stdio.h>
#include <string.h>
#endif /* !__APPLE__, !_WIN32 */
#endif /* !_WIN32 */

#ifdef __sun
#include <sys/processor.h>
#endif

#include <cstdlib>
#include <iostream>

#ifdef HAVE_IPP
#include <ippversion.h>
#include <ipp.h>
#endif

#ifdef HAVE_VDSP
#include <Accelerate/Accelerate.h>
#include <fenv.h>
#endif

#ifdef _WIN32
#include <fstream>
#endif


namespace audiomod {


#ifdef _WIN32

void gettimeofday(struct timeval *tv, void *tz)
{
    union { 
	long long ns100;  
	FILETIME ft; 
    } now; 
    
    ::GetSystemTimeAsFileTime(&now.ft); 
    tv->tv_usec = (long)((now.ns100 / 10LL) % 1000000LL); 
    tv->tv_sec = (long)((now.ns100 - 116444736000000000LL) / 10000000LL); 
}

void usleep(unsigned long usec)
{
    ::Sleep(usec == 0 ? 0 : usec < 1000 ? 1 : usec / 1000);
}

#endif

void system_specific_initialise()
{
#if defined HAVE_IPP
#ifndef USE_IPP_DYNAMIC_LIBS
#if (IPP_VERSION_MAJOR < 9)
    // This was removed in v9
    ippStaticInit();
#endif
#endif
    ippSetDenormAreZeros(1);
#elif defined HAVE_VDSP
#if defined __i386__ || defined __x86_64__ 
    fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);
#elif defined __arm64__
    fesetenv(FE_DFL_DISABLE_DENORMS_ENV);
#endif
#endif
#if defined __ARMEL__
    // ARM32
    static const unsigned int x = 0x04086060;
    static const unsigned int y = 0x03000000;
    int r;
    asm volatile (
        "fmrx	%0, fpscr   \n\t"
        "and	%0, %0, %1  \n\t"
        "orr	%0, %0, %2  \n\t"
        "fmxr	fpscr, %0   \n\t"
        : "=r"(r)
        : "r"(x), "r"(y)
	);
#endif
}




#ifdef _WIN32
void system_memorybarrier()
{
#ifdef __MSVC__
    MemoryBarrier();
#else /* (mingw) */
    LONG Barrier = 0;
    __asm__ __volatile__("xchgl %%eax,%0 "
                         : "=r" (Barrier));
#endif
}
#else /* !_WIN32 */
#if (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 1)
// Not required
#else
#include <pthread.h>
void system_memorybarrier()
{
    pthread_mutex_t dummy = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&dummy);
    pthread_mutex_unlock(&dummy);
}
#endif
#endif

}



