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

#include "biquadfilter.h"

class equalizer : public modbase
{
public:
    enum EqualizerFilterType {
        HighPassFilter = 0,
        LowShelfFilter,
        PeakingFilter_0,
        PeakingFilter_1,
        PeakingFilter_2,
        PeakingFilter_3,
        HighShelfFilter,
        LowPassFilter
    };

    equalizer(int sampleRate, int numChannels, float *paramlist = nullptr);
    ~equalizer();

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
    bool outputReady() {return outready_;}

    bool getFilterStatus(EqualizerFilterType type) const;
    float getGainValue(EqualizerFilterType type) const;
    float getQValue(EqualizerFilterType type) const;
    float getCutoffFreq(EqualizerFilterType type) const;

    int setFilterEnable(bool enable, EqualizerFilterType type);
    int setGainValue(float gain, EqualizerFilterType type);
    int setQValue(float q, EqualizerFilterType type);
    int setCutoffFreq(float cutoffFreq, EqualizerFilterType type);

private:
    void allocateEqualizerFilters();
    void deallocateEqualizerFilters();
    void filterProcessBlock(float *const * bufferData, int num_samples, biquadfilter* filter);


    // std::string config_path_;

    biquadfilter *HighPassFilter_;

    biquadfilter *LowShelfFilter_;
    
    biquadfilter *PeakingFilter_0_;

    biquadfilter *PeakingFilter_1_;

    biquadfilter *PeakingFilter_2_;

    biquadfilter *PeakingFilter_3_;

    biquadfilter *HighShelfFilter_;

    biquadfilter *LowPassFilter_;

    bool highPass_useflag_ = false;
    float highPass_cutoffFreq_;
    float highPass_q_;
    float highPass_gain_;

    bool lowShelf_useflag_ = false;
    float lowShelf_cutoffFreq_;
    float lowShelf_q_;
    float lowShelf_gain_;

    bool peaking_0_useflag_ = false;
    float peaking_0_cutoffFreq_;
    float peaking_0_q_;
    float peaking_0_gain_;

    bool peaking_1_useflag_ = false;
    float peaking_1_cutoffFreq_;
    float peaking_1_q_;
    float peaking_1_gain_;

    bool peaking_2_useflag_ = false;
    float peaking_2_cutoffFreq_;
    float peaking_2_q_;
    float peaking_2_gain_;

    bool peaking_3_useflag_ = false;
    float peaking_3_cutoffFreq_;
    float peaking_3_q_;
    float peaking_3_gain_;

    bool highShelf_useflag_ = false;
    float highShelf_cutoffFreq_;
    float highShelf_q_;
    float highShelf_gain_;

    bool lowPass_useflag_ = false;
    float lowPass_cutoffFreq_;
    float lowPass_q_;
    float lowPass_gain_;

    bool outready_;
};