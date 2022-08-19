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

#include "phasevocoderimpl.h"

#include <set>

namespace audiomod
{

class resampler;

class phasevocodercore::Impl::channelinfo
{
public:        
    /**
     * Construct a channelinfo structure.
     *
     * The sizes passed in here are for the time-domain analysis
     * window and FFT calculation, and most of the buffer sizes also
     * depend on them.  In practice they are always powers of two, the
     * window and FFT sizes are either equal or generally in a 2:1
     * relationship either way, and except for very extreme stretches
     * the FFT size is either 1024, 2048 or 4096.
     *
     * The outbuf size depends on other factors as well, including
     * the pitch scale factor and any maximum processing block
     * size specified by the user of the code.
     */
    channelinfo(size_t windowSize,
                size_t fftSize,
                size_t outbufSize);


    ~channelinfo();

    /**
     * Reset buffers
     */
    void reset();

    /**
     * get the realsize and buffer size of this object
     */
    int GetRealSize() const {return realsize_;}
    int GetBufSize() const {return bufsize_;}
    

    data_type *mag;
    data_type *phase;
    data_type *real;
    data_type *imag;

    data_type *prev_phase;
    data_type *prev_error;
    data_type *prev_outphase;
    data_type *locked_phase; // for phase-locked vocoder

    size_t prev_increment; // only used in RT mode

    float *outputAccumulator;
    float *windowAccumulator;

    float *interfacebuffer;
    data_type *internalbuffer; // owned by FFT object, only used for time domain FFT i/o
    data_type *envelope; // for cepstral formant shift

    size_t slicecnt;
    size_t incnt;
    size_t outcnt;

    FFT *fft;

    resampler *res;

    circularqueue<float> *inbuf;
    circularqueue<float> *outbuf;

private:
    
    int realsize_;
    int bufsize_;
};        

}

