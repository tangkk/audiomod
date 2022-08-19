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

#include "channelinfo.h"

#include <algorithm>

#include "../common/dsp/resampler.h"
#include "../common/system/memallocators.h"

namespace audiomod 
{
      
phasevocodercore::Impl::channelinfo::channelinfo(size_t windowSize,
                                                    size_t fftSize,
                                                    size_t outbufSize)
{
    size_t bufferSize = windowSize; 
    if (fftSize > bufferSize) bufferSize = fftSize;
    size_t realSize = bufferSize / 2 + 1; // NOTE here realSize probably = fftSize / 2 + 1
    bufferSize *= 2; // enlarge buffersize after readSize calculation
    if (outbufSize < bufferSize) outbufSize = bufferSize;
    inbuf = new circularqueue<float>(bufferSize); // this is float type because it interfaces with input data
    outbuf = new circularqueue<float>(outbufSize); // this is float type because it interfaces with output data

    interfacebuffer = mem_allocate_and_zero<float>(bufferSize); // this is float because it interfaces with inbuf and outbuf
    internalbuffer = mem_allocate_and_zero<data_type>(bufferSize); // this is data_type because it is used internally

    outputAccumulator = mem_allocate_and_zero<float>(bufferSize); // used in the OLA process for the overlapped output information, interface with interfacebuffer
    windowAccumulator = mem_allocate_and_zero<float>(bufferSize); // used in the OLA process for the overlapped window information, interface with interfacebuffer

    mag = mem_allocate_and_zero<data_type>(realSize);
    phase = mem_allocate_and_zero<data_type>(realSize);
    real = mem_allocate_and_zero<data_type>(realSize);
    imag = mem_allocate_and_zero<data_type>(realSize);

    prev_phase = mem_allocate_and_zero<data_type>(realSize);
    prev_error = mem_allocate_and_zero<data_type>(realSize);
    prev_outphase = mem_allocate_and_zero<data_type>(realSize);
    locked_phase = mem_allocate_and_zero<data_type>(realSize);
    envelope = mem_allocate_and_zero<data_type>(realSize);

    fft = new FFT(fftSize);

    const int resampleBufferSize = 4096 * 16;
    res = new audiomod::resampler(resampler::FastestTolerable, 1, resampleBufferSize); // one channel because this resampler is only linked to one channel

    

    reset();

    realsize_ = realSize;
    bufsize_ = bufferSize;

}

phasevocodercore::Impl::channelinfo::~channelinfo()
{
    delete res;

    delete inbuf;
    delete outbuf;

    mem_deallocate(mag);
    mem_deallocate(phase);
    mem_deallocate(prev_phase);
    mem_deallocate(locked_phase);
    mem_deallocate(prev_error);
    mem_deallocate(prev_outphase);
    mem_deallocate(envelope);
    // mem_deallocate(interpolator);
    // mem_deallocate(ms);
    mem_deallocate(outputAccumulator);
    mem_deallocate(windowAccumulator);
    mem_deallocate(interfacebuffer);
    mem_deallocate(internalbuffer);

}

void
phasevocodercore::Impl::channelinfo::reset()
{
    inbuf->reset();
    outbuf->reset();

    if (res) res->reset();

    size_t size = inbuf->GetSize();

    for (size_t i = 0; i < size; ++i) {
        outputAccumulator[i] = 0.f;
        windowAccumulator[i] = 0.f;
    }

    // Avoid dividing opening sample (which will be discarded anyway) by zero
    windowAccumulator[0] = 1.f;
    
    prev_increment = 0;
    slicecnt = 0;
    incnt = 0;
    outcnt = 0;

}

}
