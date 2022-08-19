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
#include <memory>
#include <mutex>

#include "modbase.h"

class biquadfilter : public modbase
{
public:

    enum Type {
        highPass = 0,
        lowShelf,
        peaking,
        notch,
        highShelf,
        lowPass,
        
        bandpass_constant_skirt,
        bandpass_constant_zero,
        allpass
    };

    biquadfilter(int sampleRate, int numChannels, Type type, float cutoffFreq, float q = 5.0, float dBGain = 0.0);
    ~biquadfilter();

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

    void updateParams(Type type, float cutoffFreq, float q, float dBGain = 0.0) {
        if (type != type_) {
            setType(type);
        }
        if (cutoffFreq != cutoffFreq_) {
            setCutoff(cutoffFreq);
        }
        if (q != q_) {
            setQ(q);
        }
        if (dBGain != dBGain_) {
            setGain(dBGain);
        }
    }

    
    float process(float x, int channel);

    Type getType() const {return type_;}
    void setType(Type type);
    void setAllParams(Type type, float cutoffFreq, float q, float dBGain); 
    float getCutoff() const {return cutoffFreq_;}
    void setCutoff(float cutoffFreq);
    float getQ() const { return q_;}
    void setQ(float q);
    float getGain() const { return dBGain_; }
    void setGain(float dBGain);
    float getCutoff();
    float getQ();
    float getGain();
protected:
    Type type_;
    void computeCoeffs();
    
    float cutoffFreq_;
    float sampleRate_;
    float q_;
    float dBGain_;
    float b0, b1, b2, a0, a1, a2;
    float *x1s;
    float *x2s;
    float *y1s;
    float *y2s;

    bool outready_;

    std::mutex coffupdate_mutex;
};