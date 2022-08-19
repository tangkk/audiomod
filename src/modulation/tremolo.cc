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

#include "tremolo.h"
#include <cmath>

#include "../common/dsp/lfo.h"

tremolo::tremolo(int sampleRate, int numChannels, float frequency, float depth, int waveform)
:frequency_(frequency), depth_(depth), waveform_(waveform) {  
    lfoPhase_ = 0.0;

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;

    inverseSampleRate_ = 1.0 / sample_rate_;
}

tremolo::~tremolo() {
    
}

void tremolo::processBlock (float *const * bufferData, int num_samples) {
    float *left;
    float *right;
    
    if (num_channels_ > 1) {
        left = bufferData[0];
        right = bufferData[1];
        process(num_samples, left);
    } else {
        left = bufferData[0];
        right = nullptr;
        process(num_samples, left, right);
    }
}

void tremolo::process(int numSamples, float* left, float *right) {
    for(int i = 0; i < numSamples; i++) {
        float factor = 1.0f - depth_*lfo_zero2one(lfoPhase_, waveform_);
        left[i] *= factor;
        if(num_channels_ > 1 && right != nullptr)
            right[i] *= factor;
        lfoPhase_ += frequency_ * inverseSampleRate_;
        if(lfoPhase_ >= 1.0)
            lfoPhase_ -= 1.0;
    }
}

