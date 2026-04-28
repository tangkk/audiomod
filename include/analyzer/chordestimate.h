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
#include <string>
#include "modbase.h"
#include "../src/chromathing/Chordino.h" 
#include "../src/common/dsp/FFT.h"

namespace audiomod {


class chordestimate : public modbase_analyzer {
public:
    chordestimate(int sampleRate, int numChannels, int blockSize, int stepSize);

    ~chordestimate();

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
    Chordino* Ch;
    std::vector<Chordino::Feature> ChordSequence; // a sequence of <onset, chords>, in terms of floats
    std::vector<Chordino::Feature> ChordNoteSequence;
    FFT *fft;
    float *realOut;
    float *imagOut;
    float realtime_f;
    float samplerate_f;
};

}