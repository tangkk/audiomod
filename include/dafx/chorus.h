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

# pragma once

#include "modbase.h"


//==============================================================================
/**
*/
class chorus : public modbase
{
public:
    //==============================================================================
    chorus(int sampleRate, int numChannels);
    ~chorus();

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

private:
    
    enum Parameters
    {
        kDelayParam = 0,
        kSweepWidthParam,
        kDepthParam,
        kFrequencyParam,
        kWaveformParam,
        kInterpolationParam,
        kNumVoicesParam,
        kStereoParam,
        kNumParameters
    };
    
    enum Interpolation
    {
        kInterpolationNearestNeighbour = 1,
        kInterpolationLinear,
        kInterpolationCubic,
        kNumInterpolations
    };
    
    static const float kMaximumDelay;
    static const float kMaximumSweepWidth;
    
    // Adjustable parameters:
    float delay_;      // Minimum length of delay line in seconds
    float sweepWidth_; // Amount of change from min to max delay
    float depth_;      // Mix level of delayed signal (0-1)
    float frequency_;  // LFO frequency (Hz)
    int   waveform_;   // What shape should be used for the LFO
    int   interpolation_; // What type of interpolation to use
    int   numVoices_;  // How many voices to use in the chorus (2-5)
    int   stereo_;     // Whether to use stereo (quadrature-phase) chorus
        
    // Circular buffer variables for implementing delay
    // float** delayBuffer_; // one or two channels
    float** delayBuffer_;

    int delayBufferLength_;
    int delayWritePosition_;
    
    float lfoPhase_;   // Phase of the low-frequency oscillator
    float inverseSampleRate_; // It's more efficient to multiply than divide, so
                               // cache the inverse of the sample rate
    
};

