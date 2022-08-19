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
#include "equalizer.h"
#include <stdio.h>

equalizer::equalizer(int sampleRate, int numChannels, float *paramlist){
    // Init params.
    sample_rate_ = sampleRate;
    num_channels_ = numChannels;

    HighPassFilter_ = nullptr;

    LowShelfFilter_ = nullptr;
    
    PeakingFilter_0_ = nullptr;

    PeakingFilter_1_ = nullptr;

    PeakingFilter_2_ = nullptr;

    PeakingFilter_3_ = nullptr;

    HighShelfFilter_ = nullptr;

    LowPassFilter_ = nullptr;

    // HighPassFilter
    highPass_useflag_ = true;
    highPass_cutoffFreq_ = 200;
    highPass_q_ = 0.3;
    highPass_gain_ = 1.0;

    if (paramlist != nullptr) {
        highPass_useflag_ = paramlist[0] > 0;
        highPass_cutoffFreq_ = paramlist[1];
        highPass_q_ = paramlist[2];
        highPass_gain_ = paramlist[3];
    }

    // LowShelfFilter
    lowShelf_useflag_ = false;
    lowShelf_cutoffFreq_ = 400;
    lowShelf_q_ = 0.3;
    lowShelf_gain_ = -1.5;

    if (paramlist != nullptr) {
        lowShelf_useflag_ = paramlist[4] > 0;
        lowShelf_cutoffFreq_ = paramlist[5];
        lowShelf_q_ = paramlist[6];
        lowShelf_gain_ = paramlist[7];
    }

    // PeakingFilter_0
    peaking_0_useflag_ = false;
    peaking_0_cutoffFreq_ = 1000;
    peaking_0_q_ = 0.3;
    peaking_0_gain_ = 1.5;

    if (paramlist != nullptr) {
        peaking_0_useflag_ = paramlist[8] > 0;
        peaking_0_cutoffFreq_ = paramlist[9];
        peaking_0_q_ = paramlist[10];
        peaking_0_gain_ = paramlist[11];
    }

    // PeakingFilter_1
    peaking_1_useflag_ = false;
    peaking_1_cutoffFreq_ = 2000;
    peaking_1_q_ = 0.3;
    peaking_1_gain_ = 1.5;

    if (paramlist != nullptr) {
        peaking_1_useflag_ = paramlist[12] > 0;
        peaking_1_cutoffFreq_ = paramlist[13];
        peaking_1_q_ = paramlist[14];
        peaking_1_gain_ = paramlist[15];
    }

    // PeakingFilter_2
    peaking_2_useflag_ = false;
    peaking_2_cutoffFreq_ = 3000;
    peaking_2_q_ = 0.3;
    peaking_2_gain_ = 1.5;

    if (paramlist != nullptr) {
        peaking_2_useflag_ = paramlist[16] > 0;
        peaking_2_cutoffFreq_ = paramlist[17];
        peaking_2_q_ = paramlist[18];
        peaking_2_gain_ = paramlist[19];
    }

    // PeakingFilter_3
    peaking_3_useflag_ = false;
    peaking_3_cutoffFreq_ = 4000;
    peaking_3_q_ = 0.3;
    peaking_3_gain_ = 1.5;

    if (paramlist != nullptr) {
        peaking_3_useflag_ = paramlist[20] > 0;
        peaking_3_cutoffFreq_ = paramlist[21];
        peaking_3_q_ = paramlist[22];
        peaking_3_gain_ = paramlist[23];
    }

    // HighShelfFilter
    highShelf_useflag_ = false;
    highShelf_cutoffFreq_ = 5000;
    highShelf_q_ = 0.3;
    highShelf_gain_ = -1.5;

    if (paramlist != nullptr) {
        highShelf_useflag_ = paramlist[24] > 0;
        highShelf_cutoffFreq_ = paramlist[25];
        highShelf_q_ = paramlist[26];
        highShelf_gain_ = paramlist[27];
    }

    // LowPassFilter
    lowPass_useflag_ = false;
    lowPass_cutoffFreq_ = 6000;
    lowPass_q_ = 0.3;
    lowPass_gain_ = 1.0;

    if (paramlist != nullptr) {
        lowPass_useflag_ = paramlist[28] > 0;
        lowPass_cutoffFreq_ = paramlist[29];
        lowPass_q_ = paramlist[30];
        lowPass_gain_ = paramlist[31];
    }

    allocateEqualizerFilters();

    outready_ = false;
}

equalizer::~equalizer() {
    deallocateEqualizerFilters();
}

void equalizer::allocateEqualizerFilters(){
    if(HighPassFilter_ != nullptr || LowShelfFilter_ != nullptr || PeakingFilter_0_ != nullptr || PeakingFilter_1_ != nullptr || PeakingFilter_2_ != nullptr || PeakingFilter_3_ != nullptr || HighShelfFilter_ != nullptr || LowPassFilter_ != nullptr){
        deallocateEqualizerFilters();
    }

    if (HighPassFilter_ == nullptr) {
        HighPassFilter_ = new biquadfilter(sample_rate_, num_channels_, biquadfilter::highPass, highPass_cutoffFreq_, highPass_q_, highPass_gain_);
    }

    if (LowShelfFilter_ == nullptr) {
        LowShelfFilter_ = new biquadfilter(sample_rate_, num_channels_, biquadfilter::lowShelf, lowShelf_cutoffFreq_, lowShelf_q_, lowShelf_gain_);
    }

    if (PeakingFilter_0_ == nullptr) {
        PeakingFilter_0_ = new biquadfilter(sample_rate_, num_channels_, biquadfilter::peaking, peaking_0_cutoffFreq_, peaking_0_q_, peaking_0_gain_);
    }

    if (PeakingFilter_1_ == nullptr) {
        PeakingFilter_1_ = new biquadfilter(sample_rate_, num_channels_, biquadfilter::peaking, peaking_1_cutoffFreq_, peaking_1_q_, peaking_1_gain_);
    }

    if (PeakingFilter_2_ == nullptr) {
        PeakingFilter_2_ = new biquadfilter(sample_rate_, num_channels_, biquadfilter::peaking, peaking_2_cutoffFreq_, peaking_2_q_, peaking_2_gain_);
    }

    if (PeakingFilter_3_ == nullptr) {
        PeakingFilter_3_ = new biquadfilter(sample_rate_, num_channels_, biquadfilter::peaking, peaking_3_cutoffFreq_, peaking_3_q_, peaking_3_gain_);
    }

    if (HighShelfFilter_ == nullptr) {
        HighShelfFilter_ = new biquadfilter(sample_rate_, num_channels_, biquadfilter::highShelf, highShelf_cutoffFreq_, highShelf_q_, highShelf_gain_);
    }

    if (LowPassFilter_ == nullptr) {
        LowPassFilter_ = new biquadfilter(sample_rate_, num_channels_, biquadfilter::lowPass, lowPass_cutoffFreq_, lowPass_q_, lowPass_gain_);
    }
        
}

void equalizer::deallocateEqualizerFilters(){
    if (HighPassFilter_ != nullptr) {
        delete HighPassFilter_;
        HighPassFilter_ = nullptr;
    }

    if (LowShelfFilter_ != nullptr) {
        delete LowShelfFilter_;
        LowShelfFilter_ = nullptr;
    }

    if (PeakingFilter_0_ != nullptr) {
        delete PeakingFilter_0_;
        PeakingFilter_0_ = nullptr;
    }

    if (PeakingFilter_1_ != nullptr) {
        delete PeakingFilter_1_;
        PeakingFilter_1_ = nullptr;
    }

    if (PeakingFilter_2_ != nullptr) {
        delete PeakingFilter_2_;
        PeakingFilter_2_ = nullptr;
    }

    if (PeakingFilter_3_ != nullptr) {
        delete PeakingFilter_3_;
        PeakingFilter_3_ = nullptr;
    }

    if (HighShelfFilter_ != nullptr) {
        delete HighShelfFilter_;
        HighShelfFilter_ = nullptr;
    }

    if (LowPassFilter_ != nullptr) {
        delete LowPassFilter_;
        LowPassFilter_ = nullptr;
    }
}

int equalizer::setCutoffFreq(float cutoffFreq, EqualizerFilterType type) {
    switch (type) {
        case HighPassFilter: 
            if(HighPassFilter_ != nullptr){
                HighPassFilter_->setCutoff(cutoffFreq);
            }
            break;
        case LowShelfFilter: 
            if(LowShelfFilter_ != nullptr){
                LowShelfFilter_->setCutoff(cutoffFreq);
            }           
            break;
        case PeakingFilter_0: 
            if(PeakingFilter_0_ != nullptr){
                PeakingFilter_0_->setCutoff(cutoffFreq);
            }   
            break;
        case PeakingFilter_1: 
            if(PeakingFilter_1_ != nullptr){
                PeakingFilter_1_->setCutoff(cutoffFreq);
            }   
            break;
        case PeakingFilter_2: 
            if(PeakingFilter_2_ != nullptr){
                PeakingFilter_2_->setCutoff(cutoffFreq);
            }   
            break;
        case PeakingFilter_3: 
            if(PeakingFilter_3_ != nullptr){
                PeakingFilter_3_->setCutoff(cutoffFreq);
            }   
            break;
        case HighShelfFilter: 
            if(HighShelfFilter_ != nullptr){
                HighShelfFilter_->setCutoff(cutoffFreq);
            }   
            break;
        case LowPassFilter: 
            if(LowPassFilter_ != nullptr){
                LowPassFilter_->setCutoff(cutoffFreq);
            }   
            break;
        default:            
            break;
    }
    return 1;
}

int equalizer::setQValue(float q, EqualizerFilterType type) {
    switch (type) {
        case HighPassFilter: 
            if(HighPassFilter_ != nullptr){
                HighPassFilter_->setQ(q);
            }
            break;
        case LowShelfFilter: 
            if(LowShelfFilter_ != nullptr){
                LowShelfFilter_->setQ(q);
            }           
            break;
        case PeakingFilter_0: 
            if(PeakingFilter_0_ != nullptr){
                PeakingFilter_0_->setQ(q);
            }   
            break;
        case PeakingFilter_1: 
            if(PeakingFilter_1_ != nullptr){
                PeakingFilter_1_->setQ(q);
            }   
            break;
        case PeakingFilter_2: 
            if(PeakingFilter_2_ != nullptr){
                PeakingFilter_2_->setQ(q);
            }   
            break;
        case PeakingFilter_3: 
            if(PeakingFilter_3_ != nullptr){
                PeakingFilter_3_->setQ(q);
            }   
            break;
        case HighShelfFilter: 
            if(HighShelfFilter_ != nullptr){
                HighShelfFilter_->setQ(q);
            }   
            break;
        case LowPassFilter: 
            if(LowPassFilter_ != nullptr){
                LowPassFilter_->setQ(q);
            }   
            break;
        default:            
            break;
    }
    return 1;
}

int equalizer::setGainValue(float gain, EqualizerFilterType type) {
    switch (type) {
        case HighPassFilter: 
            if(HighPassFilter_ != nullptr){
                HighPassFilter_->setGain(gain);
            }
            break;
        case LowShelfFilter: 
            if(LowShelfFilter_ != nullptr){
                LowShelfFilter_->setGain(gain);
            }           
            break;
        case PeakingFilter_0: 
            if(PeakingFilter_0_ != nullptr){
                PeakingFilter_0_->setGain(gain);
            }   
            break;
        case PeakingFilter_1: 
            if(PeakingFilter_1_ != nullptr){
                PeakingFilter_1_->setGain(gain);
            }   
            break;
        case PeakingFilter_2: 
            if(PeakingFilter_2_ != nullptr){
                PeakingFilter_2_->setGain(gain);
            }   
            break;
        case PeakingFilter_3: 
            if(PeakingFilter_3_ != nullptr){
                PeakingFilter_3_->setGain(gain);
            }   
            break;
        case HighShelfFilter: 
            if(HighShelfFilter_ != nullptr){
                HighShelfFilter_->setGain(gain);
            }   
            break;
        case LowPassFilter: 
            if(LowPassFilter_ != nullptr){
                LowPassFilter_->setGain(gain);
            }   
            break;
        default:            
            break;
    }
    return 1;
}

int equalizer::setFilterEnable(bool enable, EqualizerFilterType type) {
    switch (type) {
        case HighPassFilter: 
            highPass_useflag_ = enable;
            break;
        case LowShelfFilter: 
            lowShelf_useflag_ = enable;
            break;
        case PeakingFilter_0: 
            peaking_0_useflag_ = enable;
            break;
        case PeakingFilter_1: 
            peaking_1_useflag_ = enable;
            break;
        case PeakingFilter_2: 
            peaking_2_useflag_ = enable;
            break;
        case PeakingFilter_3: 
            peaking_3_useflag_ = enable;
            break;
        case HighShelfFilter: 
            highShelf_useflag_ = enable;
            break;
        case LowPassFilter: 
            lowPass_useflag_ = enable;
            break;
        default:            
            break;
    }
    return 1;
}

float equalizer::getCutoffFreq(EqualizerFilterType type) const {
    float ret = 0.0f;
    switch (type) {
        case HighPassFilter: 
            if(HighPassFilter_ != nullptr){
                ret = HighPassFilter_->getCutoff();
            }
            break;
        case LowShelfFilter: 
            if(LowShelfFilter_ != nullptr){
                ret = LowShelfFilter_->getCutoff();
            }           
            break;
        case PeakingFilter_0: 
            if(PeakingFilter_0_ != nullptr){
                ret = PeakingFilter_0_->getCutoff();
            }   
            break;
        case PeakingFilter_1: 
            if(PeakingFilter_1_ != nullptr){
                ret = PeakingFilter_1_->getCutoff();
            }   
            break;
        case PeakingFilter_2: 
            if(PeakingFilter_2_ != nullptr){
                ret = PeakingFilter_2_->getCutoff();
            }   
            break;
        case PeakingFilter_3: 
            if(PeakingFilter_3_ != nullptr){
                ret = PeakingFilter_3_->getCutoff();
            }   
            break;
        case HighShelfFilter: 
            if(HighShelfFilter_ != nullptr){
                ret = HighShelfFilter_->getCutoff();
            }   
            break;
        case LowPassFilter: 
            if(LowPassFilter_ != nullptr){
                ret = LowPassFilter_->getCutoff();
            }   
            break;
        default:            
            break;
    }
    return ret;
}

float equalizer::getQValue(EqualizerFilterType type) const {
    float ret = 0.0f;
    switch (type) {
        case HighPassFilter: 
            if(HighPassFilter_ != nullptr){
                ret = HighPassFilter_->getQ();
            }
            break;
        case LowShelfFilter: 
            if(LowShelfFilter_ != nullptr){
                ret = LowShelfFilter_->getQ();
            }           
            break;
        case PeakingFilter_0: 
            if(PeakingFilter_0_ != nullptr){
                ret = PeakingFilter_0_->getQ();
            }   
            break;
        case PeakingFilter_1: 
            if(PeakingFilter_1_ != nullptr){
                ret = PeakingFilter_1_->getQ();
            }   
            break;
        case PeakingFilter_2: 
            if(PeakingFilter_2_ != nullptr){
                ret = PeakingFilter_2_->getQ();
            }   
            break;
        case PeakingFilter_3: 
            if(PeakingFilter_3_ != nullptr){
                ret = PeakingFilter_3_->getQ();
            }   
            break;
        case HighShelfFilter: 
            if(HighShelfFilter_ != nullptr){
                ret = HighShelfFilter_->getQ();
            }   
            break;
        case LowPassFilter: 
            if(LowPassFilter_ != nullptr){
                ret = LowPassFilter_->getQ();
            }   
            break;
        default:            
            break;
    }
    return ret;
}

float equalizer::getGainValue(EqualizerFilterType type) const {
    float ret = 0.0f;
    switch (type) {
        case HighPassFilter: 
            if(HighPassFilter_ != nullptr){
                ret = HighPassFilter_->getGain();
            }
            break;
        case LowShelfFilter: 
            if(LowShelfFilter_ != nullptr){
                ret = LowShelfFilter_->getGain();
            }           
            break;
        case PeakingFilter_0: 
            if(PeakingFilter_0_ != nullptr){
                ret = PeakingFilter_0_->getGain();
            }   
            break;
        case PeakingFilter_1: 
            if(PeakingFilter_1_ != nullptr){
                ret = PeakingFilter_1_->getGain();
            }   
            break;
        case PeakingFilter_2: 
            if(PeakingFilter_2_ != nullptr){
                ret = PeakingFilter_2_->getGain();
            }   
            break;
        case PeakingFilter_3: 
            if(PeakingFilter_3_ != nullptr){
                ret = PeakingFilter_3_->getGain();
            }   
            break;
        case HighShelfFilter: 
            if(HighShelfFilter_ != nullptr){
                ret = HighShelfFilter_->getGain();
            }   
            break;
        case LowPassFilter: 
            if(LowPassFilter_ != nullptr){
                ret = LowPassFilter_->getGain();
            }   
            break;
        default:            
            break;
    }
    return ret;
}

bool equalizer::getFilterStatus(EqualizerFilterType type) const {
    bool ret = false;
    switch (type) {
        case HighPassFilter: 
            if(HighPassFilter_ != nullptr){
                ret = highPass_useflag_;
            }
            break;
        case LowShelfFilter: 
            if(LowShelfFilter_ != nullptr){
                ret = lowShelf_useflag_;
            }           
            break;
        case PeakingFilter_0: 
            if(PeakingFilter_0_ != nullptr){
                ret = peaking_0_useflag_;
            }   
            break;
        case PeakingFilter_1: 
            if(PeakingFilter_1_ != nullptr){
                ret = peaking_1_useflag_;
            }   
            break;
        case PeakingFilter_2: 
            if(PeakingFilter_2_ != nullptr){
                ret = peaking_2_useflag_;
            }   
            break;
        case PeakingFilter_3: 
            if(PeakingFilter_3_ != nullptr){
                ret = peaking_3_useflag_;
            }   
            break;
        case HighShelfFilter: 
            if(HighShelfFilter_ != nullptr){
                ret = highShelf_useflag_;
            }   
            break;
        case LowPassFilter: 
            if(LowPassFilter_ != nullptr){
                ret = lowPass_useflag_;
            }   
            break;
        default:            
            break;
    }
    return ret;
}

void equalizer::filterProcessBlock(float *const * bufferData, int num_samples, biquadfilter* filter){
    for (int i=0; i<num_channels_; i++) {
        for (int j=0; j<num_samples; j++) {
            bufferData[i][j] = filter->process(bufferData[i][j],i);
        }
    }
}

void equalizer::processBlock(float *const * bufferData, int num_samples) {
    if(highPass_useflag_){
        filterProcessBlock(bufferData, num_samples, HighPassFilter_);
    }

    if(lowShelf_useflag_){
        filterProcessBlock(bufferData, num_samples, LowShelfFilter_);
    }

    if(peaking_0_useflag_){
        filterProcessBlock(bufferData, num_samples, PeakingFilter_0_);
    }

    if(peaking_1_useflag_){
        filterProcessBlock(bufferData, num_samples, PeakingFilter_1_);
    }

    if(peaking_2_useflag_){
        filterProcessBlock(bufferData, num_samples, PeakingFilter_2_);
    }

    if(peaking_3_useflag_){
        filterProcessBlock(bufferData, num_samples, PeakingFilter_3_);
    }

    if(highShelf_useflag_){
        filterProcessBlock(bufferData, num_samples, HighShelfFilter_);
    }
    
    if(lowPass_useflag_){
        filterProcessBlock(bufferData, num_samples, LowPassFilter_);
    }

    outready_ = true;
}