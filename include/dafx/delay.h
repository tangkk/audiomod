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

#include <memory>

#include "modbase.h"

class delay : public modbase
{
public:
    delay(int sampleRate, int numChannels, float currentDelayTime, float currentMix, float currentFeedback);
    ~delay();

    /**
     * set mod params
     * @param params the key-val params to be set
     */ 
    void setParams(std::map<std::string, float> params) {

    }

    /**
     * get mod params
     * @param params the key-val params to be returned
     */ 
    void getParams(std::map<std::string, float> &params) {

    }

    void processBlock (float *const * bufferData, int num_samples);

    /****************************/

    void setDelayTime(float newDelayTime)
    {
        delaySamples *= newDelayTime / delayTime;
        delayTime = newDelayTime;
    }
    float getDelayTime() const { return delayTime; }
    void setMix(float newMix) { mix = newMix; }
    float getMix()const { return mix; }
    void setFeedback(float newFeedback) { feedback = newFeedback; }
    float getFeedback()const { return feedback; }
    void clearBuffer(void);

private:
    void process(int numSamples, float* left, float *right = nullptr, bool interleaved = false);

    int delayBufferSamples;
    int delayBufferChannels;
    int delayWritePosition;
    float delayTime;
    float delaySamples;
    float mix;
    float feedback;
    float** buffers;
};
