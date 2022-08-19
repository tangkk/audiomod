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

#include <iostream>
#include <cmath>
#include <memory>
#include <mutex>

#include "modbase.h"

class compressor : public modbase
{
public:
    compressor(int sampleRate, int numChannels, float dBThreshold = -10, float ratio = 6.0, float dBMakeUpGain = 6.0, float attackTimeMs = 10.0, float releaseTimeMs = 100.0);
    virtual ~compressor();

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

    void setThreshold(float dBThreshold) {
        std::unique_lock<std::mutex> lock(coffupdate_mutex);
        dbThreshold_ = dBThreshold;
    }
    float getThreshold() const {return dbThreshold_;}

    void setRatio(float ratio) {
        std::unique_lock<std::mutex> lock(coffupdate_mutex);
        ratio_ = ratio;
    }
    float getRatio() const {return ratio_;}

    void setGain(float dbMakeUpGain) {
        std::unique_lock<std::mutex> lock(coffupdate_mutex);
        dbMakeUpGain_ = dbMakeUpGain;
    }
    float getGain() const { return dbMakeUpGain_;}

    void setAttack(float attack) {
        std::unique_lock<std::mutex> lock(coffupdate_mutex);
        tauAttack_ = attack;
        alphaAttack_ = exp(-1/(0.001 * sample_rate_ * tauAttack_));
    }
    float getAttack()const { return tauAttack_;}

    void setRelease(float release) {
        std::unique_lock<std::mutex> lock(coffupdate_mutex);
        tauRelease_ = release;
        alphaRelease_ = exp(-1/(0.001 * sample_rate_ * tauRelease_));
    }
    float getRelease()const { return tauRelease_; }

    void processBlock(float *const * bufferData, int num_samples);

protected:

    float process(float x, int channel);
    void process(float* buffer, int channel, int numSamples);

    float dbThreshold_, dbMakeUpGain_, ratio_;
    float tauAttack_, tauRelease_, alphaAttack_, alphaRelease_;
    
    float *dbx_g;
    float *dbx_l;
    float *dby_g;
    float *dby_l;
    float *dbyL_prev;
    float *c;

    std::mutex coffupdate_mutex;
};
