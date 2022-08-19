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

class vibrato : public modbase{

public:
    vibrato(int sampleRate, int numChannels, float sweepWidth = 0.01, float frequency = 3.0);
    ~vibrato();

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

    void processBlock(float *const * bufferData, int num_samples);
    
    void setWidth(float width) { sweepWidth_ = width; }
    float getWidth()const { return sweepWidth_; }
    void setFrequency(float frequency) { frequency_ = frequency; }
    float getFrequency()const { return frequency_; }
private:

    void process (int numSamples, float* left, float *right = nullptr);

    enum Interpolation
    {
        kInterpolationNearestNeighbour = 0,
        kInterpolationLinear,
        kInterpolationCubic,
        kNumInterpolations
    };

    // Adjustable parameters:
    float sweepWidth_; // Amount of change from min to max delay
    float frequency_;  // LFO frequency (Hz)
    int   waveform_;   // What shape should be used for the LFO
    int   interpolation_; // What type of interpolation to use
    // Circular buffer variables for implementing delay
    float** delayBuffer_;
    int delayBufferChannels_;
    int delayBufferLength_;
    int delayWritePosition_;
    
    float lfoPhase_;   // Phase of the low-frequency oscillator
    // cache the inverse of the sample rate
};
