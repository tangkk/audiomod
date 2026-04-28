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

#include "pyin.h"

namespace audiomod {

pyin::pyin(int sampleRate, int numChannels, int blockSize, int stepSize, float onsetSens) {
    // note that pyin only accept mono input, if numChannels > 1, take left only
    pYA = new pYINAnalyzer(sampleRate, blockSize, stepSize, onsetSens);
}

pyin::~pyin() {
    if (pYA != nullptr) {
        delete pYA;
        pYA = nullptr;
    }
}

void pyin::processInData (float *const * inData, int num_in_samples) {
    std::cerr << "pYIN processInData..." << num_in_samples << std::endl;
    std::vector<pYINInfo> pYINRes = pYA->pYINProcess(inData[0], num_in_samples, true);
    std::cerr << "pYINRes.size():" << pYINRes.size() << std::endl;
    for(std::vector<pYINInfo>::iterator I=pYINRes.begin(),E=pYINRes.end(); I!=E; ++I) {
        // std::cerr << (*I).iStart << " , " << (*I).iDuration << " , " << (*I).iNote << std::endl;
        std::vector<float> tmp = {(*I).iStart, (*I).iDuration, (*I).iNote};
        pYINNoteSequence.push_back(tmp);
    }
    numNotes = pYINNoteSequence.size();
}

void pyin::getOutData(float *const * outData, int num_out_symbols, std::vector<std::string> *labels) {

    int i = 0;
    for (const auto &x : pYINNoteSequence) {
        std::cerr << x[0] << "," << x[1] << "," << x[2] << std::endl;
        outData[0][i++] = x[0];
        outData[0][i++] = x[1];
        outData[0][i++] = x[2];
    }
    pYINNoteSequence.clear();
}

float pyin::getScalarMeasurement() const {
    return numNotes;
}

}