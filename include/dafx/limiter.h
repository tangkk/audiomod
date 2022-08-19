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

#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <memory>

#include "modbase.h"

//add -0.01 dB to avoid clipping
#define LIMIT_OFFSET 0.01

class limiter : public modbase {
public:
    limiter(int sampleRate, int numChannels, float dBThreshold = -10.0, float dBMakeUpGain = 6.0, float attackTimeMs = 0.0, float releaseTimeMs = 100.0);
    ~limiter();

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
    
    void resetBuffer();

    float getdBMakeUpGain() const { return log10(makeUpGain_) * 20.0; }
    void setMakeUpGain(float dBMakeUpGain) { makeUpGain_ = pow(10.0, dBMakeUpGain / 20.0); }

    float getdBThreshold() const { return log10(threshold_) * 20.0 + LIMIT_OFFSET; }
    void setThreshold(float dBThreshold){ threshold_ = pow(10.0, (dBThreshold - LIMIT_OFFSET) / 20.0); }

    float getAhead() const { return ahead_; }
    void setAhead(float ahead) { ahead_ = ahead; }

    float getAttack() const {return tauAttack_; }
    void setAttack (float attack) {
        tauAttack_ = attack;
        alphaAttack_ = exp(-1.0 / (sample_rate_ * 0.001 * tauAttack_));
    }

    float getRelease() const { return tauRelease_; }
    void setRelease(float release) {
        tauRelease_ = release;
        alphaRelease_ = exp(-1.0 / (sample_rate_ * 0.001 * tauRelease_));
    }

    void processBlock(float *const * bufferData, int num_samples);
    
private:

    float process(float x, int channel);
    void process(float* buffer, int channel, int numSamples);

    float makeUpGain_, threshold_, ahead_;
    float tauAttack_, alphaAttack_, tauRelease_, alphaRelease_;
    float *xPeaks_;
    float *gains_;
    std::vector<std::list<float> > buffers_;
};
