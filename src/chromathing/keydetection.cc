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

#include "keydetection.h"
#include "RealTime.h"

namespace audiomod {

keydetection::keydetection(int sampleRate, int numChannels, int blockSize, int stepSize) {
    K = new KeyDetector(sampleRate);
    K->initialise(numChannels, stepSize, blockSize);
    K->getParameterDescriptors();
    K->getOutputDescriptors();
    K->reset();
    realtime_f = 0;
    samplerate_f = sampleRate;
}

keydetection::~keydetection() {
    K->reset();
    if (K != nullptr) {
        delete K;
        K = nullptr;
    }
}

void keydetection::processInData(float *const * inData, int num_in_samples) {
    realtime_f += float(num_in_samples) / samplerate_f;
    K->process(inData, Vamp::RealTime::fromSeconds(realtime_f));
}

void keydetection::getOutData(float *const * outData, int num_out_symbols, std::vector<std::string> *labels) {
    KeyDetector::FeatureSet fs = K->getRemainingFeatures(); // global tuning should be at fs[0]
    // "Estimated key (from C major = 1 to B major = 12 and C minor = 13 to B minor = 24)"
    // 1    2   3   4   5   6   7   8   9   10  11  12
    // C    Db  D   Eb  E   F   Gb  G   Ab  A   Bb  B
    // 13   14  15  16  17  18  19  20  21  22  23  24
    for (int i=0; i<fs[2].size(); i++) {
        // outData[0][2*i] = fs[2][i].timestamp.sec; // FIXME: just turn it to float please
        // outData[0][2*i+1] = fs[2][i].values[0];
        labels->push_back(fs[2][i].timestamp.toString());
        labels->push_back(std::to_string(int(fs[2][i].values[0])));
        labels->push_back(std::to_string(fs[2][i].values[1]));
        labels->push_back(fs[2][i].label);
    }
    labels->push_back(Vamp::RealTime::fromSeconds(realtime_f).toString());
    labels->push_back("0");
    labels->push_back("0");
    labels->push_back("End");

}

float keydetection::getScalarMeasurement() const {
    return K->getNumKeys();
}

}