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

class distortion : public modbase
{
public:
    //======================================
    
    distortion(int sampleRate, int num_channels);
    ~distortion(){};

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

    enum distortionTypeIndex {
        distortionTypeHardClipping = 0,
        distortionTypeSoftClipping = 1,
        distortionTypeExponential = 2,
        distortionTypeFullWaveRectifier = 3,
        distortionTypeHalfWaveRectifier = 4,
    };
    
    
    //======================================
    void updateParameters(int sampleRate, int num_channels, int newBlock_size, float newOutputGain, int newThreshold, float newDrive, int newTypeDistortion);
    
    
    void processBlock (float *const * bufferData, int num_samples);
    
private:
    
    float outputGain;
    int threshold;
    float drive;
    int typeDistortion;
};
