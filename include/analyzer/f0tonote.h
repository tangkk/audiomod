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
#include "../src/pyin/MonoNote.h"

namespace audiomod {


class f0tonote : public modbase_analyzer {
public:
    f0tonote(int sampleRate, int numChannels, int blockSize, int stepSize, float onsetSens = 0.7);

    ~f0tonote();

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

    void getOutData(float *const * outData, int num_out_symbols, std::vector<std::string> *labels);

private:
    int numNotes;
    MonoNote *mn;
    std::vector<std::vector<float> > pYINNoteSequence; // a sequence of <onset, duration, notes>, in terms of floats
    int m_inputSampleRate;
    int m_stepSize;
    float m_pruneThresh;
    float m_onsetSensitivity;
    float m_frameDur;
    
};

}