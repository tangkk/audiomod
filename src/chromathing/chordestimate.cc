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

#include "chordestimate.h"
#include "RealTime.h"

namespace audiomod {

chordestimate::chordestimate(int sampleRate, int numChannels, int blockSize, int stepSize) {
    Ch = new Chordino(sampleRate);
    Ch->initialise(numChannels, stepSize, blockSize);
    Ch->getParameterDescriptors();
    Ch->getOutputDescriptors();
    Ch->reset();

    fft = new FFT(blockSize);
    realOut = new float[blockSize];
    imagOut = new float[blockSize];

    realtime_f = 0;
    samplerate_f = sampleRate;

}

chordestimate::~chordestimate() {
    Ch->reset();
    if (Ch != nullptr) {
        delete Ch;
        Ch = nullptr;
    }

    if (fft != nullptr) {
        delete fft;
        delete[] realOut;
        delete[] imagOut;
    }
    fft = nullptr;
    realOut = nullptr;
    imagOut = nullptr;
}

void chordestimate::processInData(float *const * inData, int num_in_samples) {
    // Ch->process(inData, Vamp::RealTime(0,0));
    // FFT the inData's left channel
    realtime_f += float(num_in_samples) / samplerate_f;
    fft->forward(inData[0], realOut, imagOut);
    Ch->process(realOut, imagOut, Vamp::RealTime::fromSeconds(realtime_f));
}

void chordestimate::getOutData(float *const * outData, int num_out_symbols, std::vector<std::string> *labels) {
    Chordino::FeatureSet fs = Ch->getRemainingFeatures();
    // m_outputChords is on index 0 (fs[0])
    ChordSequence = fs[0]; // vector of features
    ChordNoteSequence = fs[1]; // vector of features (start, dur, notenum)
    int noteseq_len = ChordNoteSequence.size();

    // printf("1\n");
    for (const auto &x : ChordSequence) {
        std::cerr << x.timestamp << ',' << x.label << std::endl;
        labels->push_back(x.timestamp.toString());
        labels->push_back(x.label);
    }

    outData[0][0] = noteseq_len;
    int idx = 0;
    for (const auto &x : ChordNoteSequence) {
        float st = std::stof(x.timestamp.toString());
        float dur = std::stof(x.duration.toString());
        int notenum = x.values[0];
        outData[0][idx+1] = st;
        outData[0][idx+2] = dur;
        outData[0][idx+3] = notenum;
        idx += 3;
    }
    // printf("2\n");
}

float chordestimate::getScalarMeasurement() const {
    return ChordSequence.size();
}

}