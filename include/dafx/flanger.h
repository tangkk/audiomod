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

#include "modbase.h"

class flanger : public modbase {
public:
    flanger(int sampleRate, int numChannels, float delaytime, float currentmix, float currentfeedback);
    ~flanger();

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
    
    void updateParameters(int sampleRate, int block_size, float newDelayTime, float newWidth, float newFrequency, float newCurrentMix, float newCurrentFeedback, float newNumChannels);
    
    void processBlock(float *const * bufferData, int num_samples);
    
private:

    static const float kTwoPi;
    
    float **delayBuffer_;
    int delayBufferSamples;
    int delayBufferChannels;
    int delayWritePosition;
    float delay_time;
    float currentMix;
    float currentFeedback;
    
    float width;
    float lfoPhase;
    float freq;
    
};
