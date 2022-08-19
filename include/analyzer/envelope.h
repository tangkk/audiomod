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
#include <vector>
#include "modbase.h"

namespace audiomod {


class envelope : public modbase_analyzer {
public:
    envelope(int sampleRate, int numChannels);

    ~envelope();

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

    float getScalarMeasurement() const;

    void processInData (float *const * inData, int num_in_samples);

    void getOutData(float *const * outData, int num_out_symbols);

private:
    
    float thisAmp_;
    std::vector<float> ampVec;
    
};

}