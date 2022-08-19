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
class autowah : public modbase
{
public:
    enum Parameters
    {
        kBaseFrequencyParam = 0, /* Centre frequency in Hz */
        kQParam, /* Q of the resonant filter */
        kLFOFrequencyParam,
        kLFOWidthParam,
        kEnvelopeWidthParam,
        kEnvelopeAttackParam,
        kEnvelopeDecayParam,
        kNumParameters
    };
    //==============================================================================
    autowah(int sampleRate, int numChannels);
    ~autowah();

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
    
    float baseFrequency_, q_;
    float lfoFrequency_, lfoWidth_;
    float envelopeWidth_, envelopeAttack_, envelopeDecay_;
    // Methods for creating and releasing memory associated with filters
    void allocateFilters();
    void deallocateFilters();
    
    float lfoPhase_;   // Phase of the low-frequency oscillator

    float *envelopes_; // Values of signal envelopes for each channel
    int numEnvelopes_;
    
    // Convert the attack and decay time constants to a multiplier for
    // a first-order lowpass filter
    float attackMultiplier_, decayMultiplier_;
    
    biquadfilter *wahFilters_;
    float inverseSampleRate_; // Save the inverse of the sample rate for faster calculation
    
};

