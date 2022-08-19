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
#include "chorus.h"

#include <cmath>
#include <cstring>

#include "../common/dsp/lfo.h"


const float chorus::kMaximumDelay = 0.05;
const float chorus::kMaximumSweepWidth = 0.05;

//==============================================================================
chorus::chorus(int sampleRate, int numChannels) {

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;

    // Set default values:
    delay_ = .03;
    sweepWidth_ = .02;
    depth_ = 1.0;
    frequency_ = 0.2;
    waveform_ = kWaveformSine;
    interpolation_ = kInterpolationLinear;
    numVoices_ = 2;
    stereo_ = num_channels_ == 2;
    
    delayBufferLength_ = 1;
    lfoPhase_ = 0.0;
    inverseSampleRate_ = 1.0 / sample_rate_;
    
    // Start the circular buffer pointer at the beginning
    delayWritePosition_ = 0;

    delayBufferLength_ = (int)((kMaximumDelay + kMaximumSweepWidth)*sample_rate_) + 3;

    delayBuffer_ = new float*[num_channels_];
    for(int i = 0; i < num_channels_; i++) {
        delayBuffer_[i] = new float[delayBufferLength_];
        memset(delayBuffer_[i], 0, sizeof(float) * delayBufferLength_);
    }

    lfoPhase_ = 0.0;
    
}

chorus::~chorus() {
    for(int i = 0; i < num_channels_; i++) {
        if (delayBuffer_[i]) {
            delete[] delayBuffer_[i];
            delayBuffer_ = nullptr;
        }
    }
    if (delayBuffer_) {
        delete[] delayBuffer_;
        delayBuffer_ = nullptr;
    }
}


void chorus::processBlock (float *const * bufferData, int num_samples) {
    // Helpful information about this block of samples:
    const int numInputChannels = num_channels_;     // How many input channels for our effect?
    const int numOutputChannels = numInputChannels;   // How many output channels for our effect?
    const int numSamples = num_samples;          // How many samples in the buffer for this block?
    
    int channel, dpw; // dpr = delay read pointer; dpw = delay write pointer
    float dpr, currentDelay, ph;
    
    // Go through each channel of audio that's passed in. In this example we apply identical
    // effects to each channel, regardless of how many input channels there are. For some effects, like
    // a stereo chorus or panner, you might do something different for each channel.
    
    for (channel = 0; channel < numInputChannels; ++channel) {
        // channelData is an array of length numSamples which contains the audio for one channel
        float* channelData = bufferData[channel];
        
        // delayData is the circular buffer for implementing delay on this channel
        float* delayData = delayBuffer_[channel];
        
        // Make a temporary copy of any state variables declared in PluginProcessor.h which need to be
        // maintained between calls to processBlock(). Each channel needs to be processed identically
        // which means that the activity of processing one channel can't affect the state variable for
        // the next channel.
        
        dpw = delayWritePosition_;
        ph = lfoPhase_;
        
        for (int i = 0; i < numSamples; ++i) {
            const float in = channelData[i];
            float interpolatedSample = 0.0;
            float phaseOffset = 0.0;
            float weight;
            
            // chorus can have more than 2 voices (where the original, undelayed signal counts as a voice).
            // In this implementation, all voices use the same LFO, but with different phase offsets. It
            // is also possible to use different waveforms and different frequencies for each voice.
            
            for(int j = 0; j < numVoices_ - 1; ++j) {
                if(stereo_ != 0 && numVoices_ > 2) {
                    // A stereo chorus pans each voice to a different location in the stereo field.
                    // How this is done depends on the number of voices:
                    // -- 2 voices: N/A (need at least 2 delayed voices for stereo chorus)
                    // -- 3 voices: 1 voice left, 1 voice right (0, 1)
                    // -- 4 voices: 1 voice left, 1 voice centre, 1 voice right (0, 0.5, 1)
                    // -- 5 voices: 1 voice left, 1 voice left-centre,
                    //              1 voice right-centre, 1 voice right (0, 0.33, 0.66, 1)
                    
                weight = (float)j/(float)(numVoices_ - 2);
                    
                    // Left and right channels are mirrors of each other in weight
                    if(channel != 0)
                        weight = 1.0 - weight;
                }
                else
                    weight = 1.0;

                // Add the voice to the mix if it has nonzero weight
                if(weight != 0.0) {
                    // Recalculate the read pointer position with respect to the write pointer. A more efficient
                    // implementation might increment the read pointer based on the derivative of the LFO without
                    // running the whole equation again, but this format makes the operation clearer.
                    
                    currentDelay = delay_ + sweepWidth_*lfo_zero2one(fmodf(ph + phaseOffset, 1.0f), waveform_);
                    dpr = fmodf((float)dpw - (float)(currentDelay * sample_rate_) + (float)delayBufferLength_,
                                (float)delayBufferLength_);
                    
                    // In this example, the output is the input plus the contents of the delay buffer (weighted by delayMix)
                    // The last term implements a tremolo (variable amplitude) on the whole thing.
          
                    if(interpolation_ == kInterpolationLinear) {
                        // Find the fraction by which the read pointer sits between two
                        // samples and use this to adjust weights of the samples
                        float fraction = dpr - floorf(dpr);
                        int previousSample = (int)floorf(dpr);
                        int nextSample = (previousSample + 1) % delayBufferLength_;
                        interpolatedSample = fraction*delayData[nextSample]
                            + (1.0f-fraction)*delayData[previousSample];
                    } else if(interpolation_ == kInterpolationCubic) {
                        // Cubic interpolation will produce cleaner results at the expense
                        // of more computation. This code uses the Catmull-Rom variant of
                        // cubic interpolation. To reduce the load, calculate a few quantities
                        // in advance that will be used several times in the equation:
                        
                        int sample1 = (int)floorf(dpr);
                        int sample2 = (sample1 + 1) % delayBufferLength_;
                        int sample3 = (sample2 + 1) % delayBufferLength_;
                        int sample0 = (sample1 - 1 + delayBufferLength_) % delayBufferLength_;
                        
                        float fraction = dpr - floorf(dpr);
                        float frsq = fraction*fraction;
                        
                        float a0 = -0.5f*delayData[sample0] + 1.5f*delayData[sample1]
                                    - 1.5f*delayData[sample2] + 0.5f*delayData[sample3];
                        float a1 = delayData[sample0] - 2.5f*delayData[sample1]
                                    + 2.0f*delayData[sample2] - 0.5f*delayData[sample3];
                        float a2 = -0.5f*delayData[sample0] + 0.5f*delayData[sample2];
                        float a3 = delayData[sample1];
                        
                        interpolatedSample = a0*fraction*frsq + a1*frsq + a2*fraction + a3;
                    } else { // Nearest neighbour interpolation
                        // Find the nearest input sample by rounding the fractional index to the
                        // nearest integer. It's possible this will round it to the end of the buffer,
                        // in which case we need to roll it back to the beginning.
                        int closestSample = (int)floorf(dpr + 0.5f);
                        if(closestSample == delayBufferLength_)
                            closestSample = 0;
                        interpolatedSample = delayData[closestSample];
                    }

                    // Store the output sample in the buffer, which starts by containing the input sample
                    channelData[i] += depth_ * weight * interpolatedSample;
                }
                
                // 3-voice chorus uses two voices in quadrature phase (90 degrees apart). Otherwise,
                // spread the voice phases evenly around the unit circle. (For 2-voice chorus, this
                // code doesn't matter since the loop only runs once.)
                if(numVoices_ < 3)
                    phaseOffset += 0.25f;
                else
                    phaseOffset += 1.0f / (float)(numVoices_ - 1);
            }
            
            // Store the current input in the delay buffer (no feedback in a chorus, unlike a flanger).
            delayData[dpw] = in;
            
            // Increment the write pointer at a constant rate. The read pointer will move at different
            // rates depending on the settings of the LFO, the delay and the sweep width.
            
            if (++dpw >= delayBufferLength_)
                dpw = 0;

            // Update the LFO phase, keeping it in the range 0-1
            ph += frequency_*inverseSampleRate_;
            if(ph >= 1.0)
                ph -= 1.0;
        }
    }
    
    // Having made a local copy of the state variables for each channel, now transfer the result
    // back to the main state variable so they will be preserved for the next call of processBlock()
    
    delayWritePosition_ = dpw;
    lfoPhase_ = ph;
    

}



