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
#include "delay.h"
#include <cstring>
#include <math.h>

delay::delay(int sampleRate, int numChannels, float currentDelayTime, float currentMix, float currentFeedback)
        :delayBufferChannels(numChannels),
        delayTime(currentDelayTime),
        mix(currentMix),
        feedback(currentFeedback) {
    float maxDelayTime = 1.0;
    delayBufferSamples = (int)(maxDelayTime * (float)sampleRate) + 1;
    buffers = new float*[delayBufferChannels];
    for(int i = 0; i < delayBufferChannels; i++) { 
        buffers[i] = new float[delayBufferSamples];
        memset(buffers[i], 0, sizeof(float) * delayBufferSamples);
    }
    delayWritePosition = 0;
    delaySamples = delayTime * sampleRate;

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;
}

delay::~delay() {
    for(int i = 0; i < delayBufferChannels; i++) {
        if (buffers[i] != nullptr) {
            delete[] buffers[i];
        }
    }
    if (buffers != nullptr) {
        delete[] buffers;
    }
}

void delay::process(int numSamples, float* left, float *right, bool interleaved) {
    int localWritePosition = delayWritePosition;
    for (int i = 0; i < numSamples; ++i) {
        int index = i;
        if (interleaved) {
            index = i * num_channels_;
        }
        
        float readPosition = fmodf ((float)localWritePosition - delaySamples + (float)delayBufferSamples, delayBufferSamples);
        int localReadPosition = floorf (readPosition);
            
        if (localReadPosition != localWritePosition) {
            float fraction = readPosition - (float)localReadPosition;
            int delayed2Offset = interleaved ? numSamples : 1;
            
            float delayed1L = buffers[0][(localReadPosition + 0)];
            float delayed2L = buffers[0][(localReadPosition + delayed2Offset) % delayBufferSamples];
            float outL = delayed1L + fraction * (delayed2L - delayed1L);
            buffers[0][localWritePosition] = left[index] + outL * feedback;
            left[index] = left[index] + mix * outL;

            if(num_channels_ > 1 && right != nullptr) {
                float delayed1R = buffers[1][(localReadPosition + 0)];
                float delayed2R = buffers[1][(localReadPosition + delayed2Offset) % delayBufferSamples];
                float outR = delayed1R + fraction * (delayed2R - delayed1R);
                buffers[1][localWritePosition] = right[index] + outR * feedback;
                right[index] = right[index] + mix * outR;
            }
        }
        if (++localWritePosition >= delayBufferSamples)
            localWritePosition -= delayBufferSamples;
    }
    delayWritePosition = localWritePosition;
}

void delay::processBlock(float *const * bufferData, int num_samples) {
    float *left;
    float *right;
    
    if (num_channels_ > 1) {
        left = bufferData[0];
        right = bufferData[1];
        process(num_samples, left, right);
    } else {
        left = bufferData[0];
        right = nullptr;
        process(num_samples, left, right);
    }
}

void delay::clearBuffer() {
    for(int i = 0; i < delayBufferChannels; i++) {
        memset(buffers[i], 0, sizeof(float) * delayBufferSamples);
    }
}
