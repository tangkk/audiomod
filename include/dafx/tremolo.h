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

class tremolo : public modbase
{
public:
    
    tremolo(int sampleRate, int numChannels, float frequency, float depth, int waveform = 0);
    ~tremolo();

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
    
    void setFrequency(float frequency) { frequency_ = frequency; }
    float getFrequency()const { return frequency_; }
    void setDepth(float depth) { depth_ = depth; }
    float getDepth()const { return depth_; }
private:

    void process (int numSamples, float* left, float *right = nullptr);


    // Adjustable parameters:
    float frequency_;  // LFO frequency (Hz)
    float depth_;      // Depth of effect (0-1)
    int waveform_;      // What shape should be used for the LFO
    float lfoPhase_;     // Phase of the low-frequency oscillator
    float inverseSampleRate_; // It's more efficient to multiply than divide, so
    // cache the inverse of the sample rate
};
