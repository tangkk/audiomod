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

#include "ringmod.h"

#include <cmath>

#include "../common/dsp/lfo.h"

//==============================================================================
ringmod::ringmod(int sampleRate, int numChannels, float sweepwidth, float carrierfreq, float lfofreq) {
    // Set default values:
    carrierFrequency_ = carrierfreq;
    sweepWidth_ = sweepwidth;
    lfoFrequency_ = lfofreq;
    waveform_ = kWaveformSine;
    
    lfoPhase_ = carrierPhase_ = 0.0;

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;

    inverseSampleRate_ = 1.0 / sample_rate_;
 
}

ringmod::~ringmod() {
}


void ringmod::processBlock (float *const * bufferData, int num_samples) {
    
    int channel;
    float cph, lph;
    
    for (channel = 0; channel < num_channels_; ++channel) {
        // channelData is an array of length numSamples which contains the audio for one channel
        float* channelData = bufferData[channel];
        
        cph = carrierPhase_;
        lph = lfoPhase_;
        
        for (int i = 0; i < num_samples; ++i) {
            const float in = channelData[i];

            // Ring modulation is easy! Just multiply the waveform by a periodic carrier
            channelData[i] = in * sinf(2.0 * M_PI * cph);

            // Update the carrier and LFO phases, keeping them in the range 0-1
            lph += lfoFrequency_*inverseSampleRate_;
            if(lph >= 1.0) lph -= 1.0;
            cph += (carrierFrequency_ + sweepWidth_*lfo_neg2one(lph, waveform_))*inverseSampleRate_;
            if(cph >= 1.0) cph -= 1.0;
        }
    }
    
    carrierPhase_ = cph;
    lfoPhase_ = lph;
    
}
