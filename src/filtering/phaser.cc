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
#include "phaser.h"
#include <math.h>

#include "../common/dsp/lfo.h"

//==============================================================================
phaser::phaser(int sampleRate, int numChannels) {
    // Set default values:
    baseFrequency_ = 2000.0;
    sweepWidth_ = 600.0;
    depth_ = 0.6;
    feedback_ = 0.6;
    lfoFrequency_ = 6.0;
    waveform_ = kWaveformSine;
    stereo_ = 0;
    
    // Start with no filters (at least until we have some channels)
    allpassFilters_ = 0;
    filtersPerChannel_ = 4;
    totalNumFilters_ = 0;
    lastFilterOutputs_ = 0;
    numLastFilterOutputs_ = 0;

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;
    
    lfoPhase_ = 0.0;
    inverseSampleRate_ = 1.0/sample_rate_;
    sampleCount_ = 0;
    filterUpdateInterval_ = 8;

    allocateFilters();
  
}

phaser::~phaser() {
    deallocateFilters();
}


void phaser::processBlock (float *const * bufferData, int num_samples) {
    float ph, channel0EndPhase = lfoPhase_;
    unsigned int sc;
    
    // Go through each channel of audio that's passed in, applying one or more allpass filters
    // to each. Each channel will be treated identically in a (non-stereo) phaser, but we have
    // to have separate filter objects for each channel since the filters store the last few samples
    // passed through them.
    
    // Filters are stored with all channel 0 filters first, then all channel 1 filters, etc.
    
    for(int channel = 0; channel < num_channels_; ++channel) {
        // channelData is an array of length numSamples which contains the audio for one channel
        float* channelData = bufferData[channel];
        
        ph = lfoPhase_;
        sc = sampleCount_;
        
        // For stereo phasing, keep the channels 90 degrees out of phase with each other
        if(stereo_ != 0 && channel != 0)
            ph = fmodf(ph + 0.25, 1.0);
        
        for (int sample = 0; sample < num_samples; ++sample) {
            float out = channelData[sample];
            
            // If feedback is enabled, include the feedback from the last sample in the
            // input of the allpass filter chain. This is actually not accurate to how
            // analog phasers work because there is a sample of delay between output and
            // input, which adds a further phase shift of up to 180 degrees at half the
            // sampling frequency. To truly model an analog phaser with feedback involves
            // modelling a delay-free loop, which is beyond the scope of this example.
            
            if(feedback_ != 0.0 && channel < numLastFilterOutputs_)
                out += feedback_ * lastFilterOutputs_[channel];
            
            // See OnePoleAllpassFilter.cpp for calculation of coefficients and application
            // of filter to audio data. The filter processes the audio buffer in place,
            // putting the output sample in place of the input.
            
            for(int j = 0; j < filtersPerChannel_; ++j) {
                // Safety check
                // if(channel * filtersPerChannel_ + j >= totalNumFilters_)
                //     continue;
                
                // First, update the current allpass filter coefficients depending on the parameter
                // settings and the LFO phase
                
                // Recalculating the filter coefficients is much more expensive than calculating
                // a sample. Only update the coefficients at a fraction of the sample rate; since
                // the LFO moves slowly, the difference won't generally be audible.
                if(sc % filterUpdateInterval_ == 0) {
                    allpassFilters_[j]->setCutoff(
                       baseFrequency_ + sweepWidth_*lfo_zero2one(ph, waveform_));
                }
                out = allpassFilters_[j]->process(out,channel);
            }
            
            if(channel < numLastFilterOutputs_)
                lastFilterOutputs_[channel] = out;
            
            // Add the allpass signal to the output, though maintaining constant level
            // depth = 0 --> input only ; depth = 1 --> evenly balanced input and output
            channelData[sample] = (1.0f-0.5f*depth_)*channelData[sample] + 0.5f*depth_*out;
        
            // Update the LFO phase, keeping it in the range 0-1
            ph += lfoFrequency_*inverseSampleRate_;
            if(ph >= 1.0)
                ph -= 1.0;
            sc++;
        }

        // Use channel 0 only to keep the phase in sync between calls to processBlock()
        // Otherwise quadrature phase on multiple channels will create problems.
        if(channel == 0)
            channel0EndPhase = ph;
    }
    
    lfoPhase_ = channel0EndPhase;
    sampleCount_ = sc;
    

}

void phaser::allocateFilters() {
    // Create any filters we need; depends on number of channels and number of
    // filters per channel
    // totalNumFilters_ = num_channels_ * filtersPerChannel_;
    totalNumFilters_ = filtersPerChannel_;
    if(totalNumFilters_ > 0) {
        allpassFilters_ = (biquadfilter**)malloc(totalNumFilters_ * sizeof(biquadfilter*));
        for(int i = 0; i < totalNumFilters_; i++)
            allpassFilters_[i] = new biquadfilter(sample_rate_, num_channels_, biquadfilter::allpass, baseFrequency_);
    }
    
    numLastFilterOutputs_ = num_channels_;
    lastFilterOutputs_ = (float *)malloc(numLastFilterOutputs_ * sizeof(float));
    for(int i = 0; i < numLastFilterOutputs_; i++)
        lastFilterOutputs_[i] = 0.0f;
    
    // Coefficients of allpass filters will get updated in processBlock()
}

void phaser::deallocateFilters() {
    // Release the filters that were created in prepareToPlay()
    for(int i = 0; i < totalNumFilters_; i++)
        delete allpassFilters_[i];
    if(totalNumFilters_ != 0)
        free(allpassFilters_);
    totalNumFilters_ = 0;
    allpassFilters_ = 0;
    
    if(numLastFilterOutputs_ != 0)
        free(lastFilterOutputs_);
    numLastFilterOutputs_ = 0;
    lastFilterOutputs_ = 0;
}

// Release and recreate the filters in one atomic operation:
// the ScopedLock will not let the audio thread run between
// release and allocation
void phaser::reallocateFilters() {
    deallocateFilters();
    allocateFilters();
}

