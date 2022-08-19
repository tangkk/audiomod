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

#include "biquadfilter.h"
#include "modbase.h"


//==============================================================================
/**
*/
class phaser : public modbase
{
public:
    //==============================================================================
    phaser(int sampleRate, int numChannels);
    ~phaser();

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
    
    
    enum Parameters
    {
        kBaseFrequencyParam = 0,
        kSweepWidthParam,
        kDepthParam,
        kFeedbackParam,
        kLFOFrequencyParam,
        kFiltersParam,
        kWaveformParam,
        kStereoParam,
        kNumParameters
    };

    
private:
    // Adjustable parameters:
    float baseFrequency_; // Lowest frequency of allpass filters
    float sweepWidth_;    // Amount of change from min to max delay
    float depth_;         // Mix level for phase-shifted signal (0-1)
    float feedback_;      // Feedback level for feedback phaser (0-<1)
    float lfoFrequency_;  // LFO frequency (Hz)
    int   filtersPerChannel_; // How many allpass filters to use
    int   waveform_;      // What shape should be used for the LFO
    int   stereo_;        // Whether to use stereo phasing

    void allocateFilters();   // Create the filter objects...
    void deallocateFilters(); // Delete them...
    void reallocateFilters(); // Delete and rebuild in one combined operation
    
    float lfoPhase_;   // Phase of the low-frequency oscillator
    float inverseSampleRate_; // It's more efficient to multiply than divide, so
                               // cache the inverse of the sample rate
    unsigned int sampleCount_; // Counter to keep track of samples elapsed, to
                               // keep track of filter updates
    unsigned int filterUpdateInterval_; // How often to update filter coefficients
    
    // Bank of allpass filters that do the phasing; N filters x M channels
    biquadfilter **allpassFilters_;
    
    // Storage of the last output sample from each bank of filters, for use in
    // feedback loop
    float *lastFilterOutputs_;
    int numLastFilterOutputs_;

    int totalNumFilters_;

};


