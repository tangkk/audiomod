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

#include "dynamicfilter.h"

dynamicfilter::dynamicfilter(int sampleRate, int numChannels, float cutoffFreq, float q)
:biquadfilter(sampleRate, numChannels, biquadfilter::peaking, cutoffFreq, q), compressor(numChannels, sampleRate) {
}

dynamicfilter::~dynamicfilter() {
    
}

float dynamicfilter::process(float x, int channel) {
    compressor::process(x, channel);
    biquadfilter::setGain(dbMakeUpGain_ - dby_l[channel]);
    return biquadfilter::process(x, channel);
}

void dynamicfilter::processBlock(float *const * bufferData, int num_samples) {
    for (int i=0; i<biquadfilter::num_channels_; i++) {
        for (int j=0; j<num_samples; j++) {
            process(bufferData[i][j], i);
        }
    }
}
