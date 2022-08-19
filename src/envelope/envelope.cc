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


#include "envelope.h"
#include <iostream>
#include <cmath>
#include <cstring>


namespace audiomod {

envelope::envelope(int sampleRate, int numChannels) {

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;

    thisAmp_ = 0;
}

envelope::~envelope() {
    ampVec.clear();
    thisAmp_ = 0;
}

void envelope::processInData (float *const * inData, int num_in_samples) {
    float tmpAmp = 0;

    // compute RMS value
    for (int i=0; i < num_channels_; i++) {
        for (int j=0; j<num_in_samples; j++) {
            tmpAmp += inData[i][j] * inData[i][j];
        }
    }

    thisAmp_ = sqrt(tmpAmp / num_in_samples / num_channels_);
    ampVec.push_back(thisAmp_);
}

void envelope::getOutData(float *const * outData, int num_out_symbols) {
    // do nothing here for the moment
    outData[0][0] = thisAmp_;
}

float envelope::getScalarMeasurement() const {
    float sum = 0;
    for (int i=0; i<ampVec.size(); i++) {
        sum += ampVec[i];
    }
    return sum / ampVec.size();
}

} // end of namespace