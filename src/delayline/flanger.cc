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
#include "flanger.h"

#include <cmath>
#include <cstring>

#include "../common/dsp/lfo.h"

// const float flanger::kTwoPi = 2.0f * M_PI;

flanger::flanger(int sampleRate, int numChannels, float delaytime, float currentmix, float currentfeedback):
        delay_time(delaytime), currentMix(currentmix), currentFeedback(currentfeedback) {
    float maxDelayTime = 1.0;
    delayBufferChannels = numChannels;
    delayBufferSamples = (int)(maxDelayTime * (float)sampleRate) + 1;
    delayBuffer_ = new float*[delayBufferChannels];
    for(int i = 0; i < delayBufferChannels; i++) { 
        delayBuffer_[i] = new float[delayBufferSamples];
        memset(delayBuffer_[i], 0, sizeof(float) * delayBufferSamples);
    }
    delayWritePosition = 0;

    width = 0.001;
    lfoPhase = 90;
    freq = 6;

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;
}

flanger::~flanger() {
    for(int i = 0; i < delayBufferChannels; i++) {
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


void flanger::updateParameters(int sampleRate, int block_size, float newDelayTime, float newWidth, float newFrequency,
                        float newCurrentMix, float newCurrentFeedback, float newNumChannels) {
    float maxDelayTime = 1.0;
    delayBufferSamples = (int)(maxDelayTime * (float)sampleRate) + 1;
    if (delayBufferSamples < 1)
        delayBufferSamples = 1;
    
    delayBuffer_ = new float*[num_channels_];
    for(int i = 0; i < num_channels_; i++) {
        delayBuffer_[i] = new float[delayBufferSamples];
        memset(delayBuffer_[i], 0, sizeof(float) * delayBufferSamples);
    }
    
    delayWritePosition = 0;
    delay_time = newDelayTime/1000;
    currentMix = newCurrentMix;
    currentFeedback = newCurrentFeedback;
    width = newWidth/1000;
    lfoPhase = 0.0;
    freq = newFrequency;
    num_channels_ = newNumChannels;
}


void flanger::processBlock(float *const * bufferData, int num_samples) {
    float delay_sample;
    delay_sample = delay_time * sample_rate_;

    int localWritePosition;
    float phase;
    float phaseMain;
    const int numSamples = num_samples;
    const int numInputChannels = num_channels_;
    const int numOutputChannels = num_channels_;

    for (int channel = 0; channel < numInputChannels; ++channel) {
        float* channelData = bufferData[channel];
        float* delayData = delayBuffer_[channel];
        localWritePosition = delayWritePosition;
        phase = lfoPhase;
        if (num_channels_==2 && channel != 0)
            phase = fmodf (phase + 0.25f, 1.0f);

        for (int sample = 0; sample < numSamples; ++sample) {

            const float in = channelData[sample];
            float out = 0.0f;

            float localDelayTime = (delay_time + width * lfo_zero2one(phase, 0)) * sample_rate_;

            float readPosition =
            fmodf ((float)localWritePosition - localDelayTime + (float)delayBufferSamples, delayBufferSamples);
            int localReadPosition = floorf (readPosition);

            float fraction = readPosition - (float)localReadPosition;
            float delayed1 = delayData[(localReadPosition + 0)];
            float delayed2 = delayData[(localReadPosition + 1) % delayBufferSamples];
            out = delayed1 + fraction * (delayed2 - delayed1);

            channelData[sample] = in + out * currentMix;
            delayData[localWritePosition] = in + out * currentFeedback;

            phase += freq * 1.0f / sample_rate_;
            if (phase >= 1.0f)
                phase -= 1.0f;

            if (++localWritePosition >= delayBufferSamples)
                localWritePosition -= delayBufferSamples;
        }
        if (channel == 0)
            phaseMain = phase;
    }

    delayWritePosition = localWritePosition;

    lfoPhase = phaseMain;
}

