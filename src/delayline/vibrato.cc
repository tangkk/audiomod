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
#include "vibrato.h"

#include <cmath>
#include <cstring>

#include "../common/dsp/lfo.h"

vibrato::vibrato(int sampleRate, int numChannels, float sweepWidth, float frequency)
:delayBufferChannels_(numChannels), sweepWidth_(sweepWidth), frequency_(frequency) {
    waveform_ = kWaveformSine;
    interpolation_ = kInterpolationLinear;
    
    lfoPhase_ = 0.0;
    sample_rate_ = sampleRate;
    num_channels_ = numChannels;

    delayWritePosition_ = 0;
    delayBuffer_ = new float*[numChannels];
    delayBufferLength_ = (int)(0.05 * sampleRate) + 3;
    for(int i = 0; i < numChannels; i++) {
        delayBuffer_[i] = new float[delayBufferLength_];
        memset(delayBuffer_[i], 0, sizeof(float) * delayBufferLength_);
    }

}

vibrato::~vibrato() {
    for(int i = 0; i < delayBufferChannels_; i++) {
        if (delayBuffer_[i]) {
            delete[] delayBuffer_[i];
            delayBuffer_[i] = nullptr;
        }
    }
    if (delayBuffer_) {
        delete[] delayBuffer_;
        delayBuffer_ = nullptr;
    }
}

void vibrato::processBlock(float *const * bufferData, int num_samples) {
    float *left;
    float *right;
    
    if (num_channels_ > 1) {
        left = bufferData[0];
        right = bufferData[1];
        process(num_samples, left, right);
    } else {
        left = bufferData[0];
        right = nullptr;
        process(num_samples, left);
    }
}

void vibrato::process(int numSamples, float* left, float *right) {
    int dpw;
    float dpr, currentDelay, ph;
    dpw = delayWritePosition_;
    ph = lfoPhase_;
        
    for (int i = 0; i < numSamples; ++i) {
        float interpolatedSampleL = 0.0;
        float interpolatedSampleR = 0.0;
            
        // Recalculate the read pointer position with respect to the write pointer. A more efficient
        // implementation might increment the read pointer based on the derivative of the LFO without
        // running the whole equation again, but this format makes the operation clearer.
            
        currentDelay = sweepWidth_*lfo_zero2one(ph, waveform_);
            
        // Subtract 3 samples to the delay pointer to make sure we have enough previously written
        // samples to interpolate with
        dpr = fmodf((float)dpw - (float)(currentDelay * sample_rate_) + (float)delayBufferLength_ - 3.0,
                        (float)delayBufferLength_);
            
            // In this example, the output is the input plus the contents of the delay buffer (weighted by delayMix)
            // The last term implements a tremolo (variable amplitude) on the whole thing.
            
        if(interpolation_ == kInterpolationLinear) {
            // Find the fraction by which the read pointer sits between two
            // samples and use this to adjust weights of the samples
            float fraction = dpr - floorf(dpr);
            int previousSample = (int)floorf(dpr);
            int nextSample = (previousSample + 1) % delayBufferLength_;
            interpolatedSampleL = fraction * delayBuffer_[0][nextSample]
                + (1.0f-fraction) * delayBuffer_[0][previousSample];
            if(num_channels_ > 1 && right != nullptr)
                interpolatedSampleR = fraction * delayBuffer_[1][nextSample]
                + (1.0f-fraction) * delayBuffer_[1][previousSample];
        }
        else if(interpolation_ == kInterpolationCubic) {
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
                
            float a0 = -0.5f * delayBuffer_[0][sample0] + 1.5f * delayBuffer_[0][sample1]
                - 1.5f * delayBuffer_[0][sample2] + 0.5f * delayBuffer_[0][sample3];
            float a1 = delayBuffer_[0][sample0] - 2.5f * delayBuffer_[0][sample1]
                + 2.0f * delayBuffer_[0][sample2] - 0.5f * delayBuffer_[0][sample3];
            float a2 = -0.5f * delayBuffer_[0][sample0] + 0.5f * delayBuffer_[0][sample2];
            float a3 = delayBuffer_[0][sample1];
                
            interpolatedSampleL = a0*fraction*frsq + a1*frsq + a2*fraction + a3;
            if(num_channels_ > 1 && right != nullptr) {
                a0 = -0.5f * delayBuffer_[1][sample0] + 1.5f * delayBuffer_[1][sample1]
                    - 1.5f * delayBuffer_[1][sample2] + 0.5f * delayBuffer_[1][sample3];
                a1 = delayBuffer_[1][sample0] - 2.5f * delayBuffer_[1][sample1]
                    + 2.0f * delayBuffer_[1][sample2] - 0.5f * delayBuffer_[1][sample3];
                a2 = -0.5f * delayBuffer_[1][sample0] + 0.5f * delayBuffer_[1][sample2];
                a3 = delayBuffer_[1][sample1];
                interpolatedSampleR = a0*fraction*frsq + a1*frsq + a2*fraction + a3;
            }
        }
        else // Nearest neighbour interpolation
        {
            // Find the nearest input sample by rounding the fractional index to the
            // nearest integer. It's possible this will round it to the end of the buffer,
            // in which case we need to roll it back to the beginning.
            int closestSample = (int)floorf(dpr + 0.5);
            if(closestSample == delayBufferLength_)
                closestSample = 0;
            interpolatedSampleL = delayBuffer_[0][closestSample];
            if(num_channels_ > 1 && right != nullptr)
                interpolatedSampleR = delayBuffer_[0][closestSample];
        }
            
        // Store the current information in the delay buffer. With feedback, what we read is
        // included in what gets stored in the buffer, otherwise it's just a simple delay line
        // of the input signal.
        delayBuffer_[0][dpw] = left[i];
        if(num_channels_ > 1 && right != nullptr)
            delayBuffer_[1][dpw] = right[i];
            
        // Increment the write pointer at a constant rate. The read pointer will move at different
        // rates depending on the settings of the LFO, the delay and the sweep width.
            
        if (++dpw >= delayBufferLength_)
            dpw = 0;
            
        // Store the output sample in the buffer, replacing the input. In the vibrato effect,
        // the delaye sample is the only component of the output (no mixing with the dry signal)
        left[i] = interpolatedSampleL;
        if(num_channels_ > 1 && right != nullptr)
            right[i] = interpolatedSampleR;
        // Update the LFO phase, keeping it in the range 0-1
        ph += frequency_ / sample_rate_;
        if(ph >= 1.0)
            ph -= 1.0;
    }
    // Having made a local copy of the state variables for each channel, now transfer the result
    // back to the main state variable so they will be preserved for the next call of processBlock()
    
    delayWritePosition_ = dpw;
    lfoPhase_ = ph;
}


