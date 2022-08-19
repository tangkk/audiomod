
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
#include "autowah.h"

#include <cmath>

// The filter will produce a resonant peak of amplitude Q; bring everything
// down somewhat to compensate, though try to maintain some perceptual balance
// of being similar loudness. (This factor has been chosen somewhat arbitrarily.)
const float kWahwahFilterGain = 0.5;

//==============================================================================
autowah::autowah(int sampleRate, int numChannels) {
    // Set default values:
    baseFrequency_ = 600.0;
    q_ = 5.0;
    lfoFrequency_ = 2.0;
    lfoWidth_ = 1000.0;
    envelopeWidth_ = 0.0;
    envelopeAttack_ = 0.005;
    envelopeDecay_ = 0.1;

    lfoPhase_ = 0;
    
    // Initialise the filters later when we know how many channels
    wahFilters_ = nullptr;
    envelopes_ = nullptr;

    numEnvelopes_ = 0;
    attackMultiplier_ = 1.0;
    decayMultiplier_ = 0.0;

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;
    
    inverseSampleRate_ = 1.0 / sample_rate_; // start with a sensible default

    allocateFilters();

    if(envelopeDecay_ == 0.0)
        decayMultiplier_ = 0.0;
    else
        decayMultiplier_ = pow(1.0 / M_E, inverseSampleRate_ / envelopeDecay_);
    if(envelopeAttack_ == 0.0)
        attackMultiplier_ = 0.0;
    else
        attackMultiplier_ = pow(1.0 / M_E, inverseSampleRate_ / envelopeAttack_);
}

autowah::~autowah() {
    deallocateFilters();
}


void autowah::processBlock (float *const * bufferData, int num_samples) {
    // Helpful information about this block of samples:
    float ph;

    const int numInputChannels = num_channels_;
    const int numOutputChannels = num_channels_;
    
    // Go through each channel and put it through the resonant lowpass filter, updating
    // the coefficients as we go along. Each channel is processed identically in this effect.
    
    for(int channel = 0; channel < num_channels_; ++channel) {
        // channelData is an array of length num_samples which contains the audio for one channel
        float* channelData = bufferData[channel];
        ph = lfoPhase_;
        
        for (int sample = 0; sample < num_samples; ++sample) {
            const float in = channelData[sample];
            float centreFrequency = baseFrequency_;
            
            // Calculate the envelope of the signal. Do this even if we're not currently
            // changing the frequeny based on it, since it involves maintaining a history
            // of the signal's behaviour.
            
            if(channel < numEnvelopes_) {   // Safety check
                if(fabs(in) > envelopes_[channel]) {
                    envelopes_[channel] += (1.0 - attackMultiplier_) * (fabs(in) - (float)envelopes_[channel]);
                }
                else
                    envelopes_[channel] *= decayMultiplier_;
            }
            
            // Calculate the centre frequency of the filter based on the LFO and the
            // signal envelope
            if(lfoWidth_ > 0.0) {
                centreFrequency += lfoWidth_ * (0.5f + 0.5f*sinf(2.0 * M_PI * ph));
            }
            if(envelopeWidth_ > 0.0 && channel < numEnvelopes_) {
                centreFrequency += envelopeWidth_ * envelopes_[channel];
            }
            

            wahFilters_->setCutoff(centreFrequency);
            
            // Process one sample and store it back in place. See juce_IIRFilter.cpp for the
            // application of the IIR filter.
            channelData[sample] = wahFilters_->process(in,channel);
            
            // Update the LFO phase, keeping it in the range 0-1
            ph += lfoFrequency_*inverseSampleRate_;
            if(ph >= 1.0)
                ph -= 1.0;
        }
    }
    
    lfoPhase_ = ph;
    

}

void autowah::allocateFilters()
{
    // Prevent leaks from reallocation
    if(wahFilters_ != nullptr || envelopes_ != nullptr)
        deallocateFilters();
    
    // Create as many filters as we have input channels


    if (wahFilters_ == nullptr) {
        wahFilters_ = new biquadfilter(sample_rate_, num_channels_, biquadfilter::lowPass, baseFrequency_, q_);
    }
    
    numEnvelopes_ = num_channels_;
    envelopes_ = (float *)malloc(numEnvelopes_ * sizeof(float));
    if(envelopes_ == nullptr)
        numEnvelopes_ = 0;
    else {
        for(int i = 0; i < numEnvelopes_; i++)
            envelopes_[i] = 0.0;
    }
}

void autowah::deallocateFilters()
{

    if (wahFilters_ != nullptr) {
        delete wahFilters_;
        wahFilters_ = nullptr;
    }

    if(envelopes_ != nullptr) {
        free(envelopes_);
        envelopes_ = nullptr;
        numEnvelopes_ = 0;
    }
}

