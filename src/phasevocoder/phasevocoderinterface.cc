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

#include "phasevocoderimpl.h"

using namespace std;

namespace audiomod {


phasevocodercore::phasevocodercore(size_t sampleRate,
                                         size_t channels,
                                         Options options,
                                         double initialTimeRatio,
                                         double initialPitchScale) :
    m_pv(new Impl(sampleRate, channels, options,
                 initialTimeRatio, initialPitchScale)) {
}

phasevocodercore::~phasevocodercore() {
    delete m_pv;
}


double
phasevocodercore::getTimeRatio() const {
    return m_pv->getTimeRatio();
}

double
phasevocodercore::getPitchScale() const {
    return m_pv->getPitchScale();
}

void
phasevocodercore::processNormal(const float *const *input, size_t num_samples) {
    m_pv->processNormal(input, num_samples);
}

void
phasevocodercore::processConstant(const float *const *input, size_t num_samples) {
    m_pv->processConstant(input, num_samples);
}

void 
phasevocodercore::processVocoder(const float *const *input, size_t num_samples, int carrierType) {
    m_pv->processVocoder(input, num_samples, static_cast<phasevocodercore::Impl::carrierType>(carrierType));
}

int
phasevocodercore::numsamples_available() const {
    return m_pv->numsamples_available();
}

int
phasevocodercore::numsamples_availableCarrier() const {
    return m_pv->numsamples_availableCarrier();
}

size_t
phasevocodercore::retrieve(float *const *output, size_t num_samples) const {
    return m_pv->retrieve(output, num_samples);
}

size_t
phasevocodercore::retrieveCarrier(float *const *output, size_t num_samples) const {
    return m_pv->retrieveCarrier(output, num_samples);
}

size_t
phasevocodercore::getInputHopSize() const {
    return m_pv->getInputHopSize();
}

size_t
phasevocodercore::getChannelCount() const {
    return m_pv->getChannelCount();
}


void
phasevocodercore::setDebugLevel(int level) {
    m_pv->setDebugLevel(level);
}

void
phasevocodercore::setDefaultDebugLevel(int level) {
    Impl::setDefaultDebugLevel(level);
}

void
phasevocodercore::setDefaultFftSize(int fftsize) {
    Impl::setDefaultFftSize(fftsize);
}

void
phasevocodercore::setDefaultHopSize(int hopsize) {
    Impl::setDefaultHopSize(hopsize);
}

void
phasevocodercore::setDefaultCoreMode(int coremode) {
    Impl::setDefaultCoreMode(coremode);
}

} // end of namespace audiomod

