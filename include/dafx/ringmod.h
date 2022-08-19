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

class ringmod : public modbase
{
public:

    //==============================================================================
    ringmod(int sampleRate, int numChannels, float sweepwidth = 0.1, float carrierfreq = 200.0, float lfofreq = 2);
    ~ringmod();

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
        kCarrierFrequencyParam = 0,
        kSweepWidthParam,
        kLFOFrequencyParam,
        kWaveformParam,
        kNumParameters
    };
    
    
    // Adjustable parameters:
    float carrierFrequency_; // Frequency of main modulator (Hz)
    float sweepWidth_;    // Amount of change from min to max delay
    float lfoFrequency_;  // LFO frequency (Hz)
    int   waveform_;      // What shape should be used for the LFO
    
    float lfoPhase_;     // Phase of the low-frequency oscillator
    float carrierPhase_; // Phase of the main (carrier) oscillator
    
    float inverseSampleRate_; // It's more efficient to multiply than divide, so
                               // cache the inverse of the sample rate

};

