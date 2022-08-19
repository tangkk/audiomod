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

#pragma once

#include "../system/sys.h"

#define NO_EXCEPTIONS // temp fixed def

namespace audiomod {

class resamplerImpl;

class resampler
{
public:
    enum RSMode { Best, FastestTolerable, Fastest };
    // enum Exception { ImplementationError };

    /**
     * Construct a resampler with the given quality level and channel
     * count.  maxBufferSize gives a bound on the maximum incount size
     * that may be passed to the resample function before the
     * resampler needs to mem_reallocate its internal buffers.
     */
    resampler(RSMode quality, int channels, int maxBufferSize = 0);
    ~resampler();

    /**
     * Resample the given multi-channel buffers, where incount is the
     * number of frames in the input buffers.  Returns the number of
     * frames written to the output buffers.
     */
    int doresample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final = false);

    int numchannels() const;

    void reset();

protected:
    resamplerImpl *rs;
    int method_;
};

}

