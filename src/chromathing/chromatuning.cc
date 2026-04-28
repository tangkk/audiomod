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

#include "chromatuning.h"
#include "RealTime.h"

namespace audiomod {

chromatuning::chromatuning(int sampleRate, int numChannels, int blockSize, int stepSize) {
    T = new Tuning(sampleRate);
    T->initialise(numChannels, stepSize, blockSize);
    T->getParameterDescriptors();
    T->getOutputDescriptors();
    T->reset();

    fft = new FFT(blockSize);
    realOut = new float[blockSize];
    imagOut = new float[blockSize];
}

chromatuning::~chromatuning() {
    T->reset();
    if (T != nullptr) {
        delete T;
        T = nullptr;
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

void chromatuning::processBlock(float *const * inData, int num_in_samples) {
    // FFT the inData's left channel
    fft->forward(inData[0], realOut, imagOut);
    T->process(realOut, imagOut, Vamp::RealTime(0,0));
}

float chromatuning::getScalarMeasurement() const {
    Tuning::FeatureSet fs = T->getRemainingFeatures(); // global tuning should be at fs[0]
    std::cerr << "global Tuning:" << fs[0][0].values[0] << std::endl; //cumulativetuning
    return fs[0][0].values[0];
}

}