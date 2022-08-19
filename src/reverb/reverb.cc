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

#include "reverb.h"
#include "revmodel.h"

reverb::reverb(int sampleRate, int numChannels, float roomsize, float damp, float width, float dry, float wet) {
    sample_rate_ = sampleRate;
    num_channels_ = numChannels;

    rev = new revmodel(sample_rate_);

    rev->setroomsize(roomsize);
    rev->setdamp(damp);
    rev->setwidth(width);
    rev->setdry(dry);
    rev->setwet(wet);
    

}

reverb::~reverb() {
    if (rev != nullptr) {
        delete rev;
        rev = nullptr;
    }
}

void reverb::processBlock(float *const * bufferData, int num_samples) {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    if (num_channels_ == 1) {
        rev->processreplace(bufferData[0], nullptr, bufferData[0], nullptr, num_samples);
    } else if (num_channels_ == 2) {
        rev->processreplace(bufferData[0], bufferData[1], bufferData[0], bufferData[1], num_samples);
    }
}

void	reverb::setroomsize(float value) {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    rev->setroomsize(value);
}
float	reverb::getroomsize() const {
    return rev->getroomsize();
}
void	reverb::setdamp(float value) {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    rev->setdamp(value);
}
float	reverb::getdamp() const {
    return rev->getdamp();
}
void	reverb::setwet(float value) {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    rev->setwet(value);
}
float	reverb::getwet() const {
    return rev->getwet();
}
void	reverb::setdry(float value) {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    rev->setdry(value);
}
float	reverb::getdry() const {
    return rev->getdry();
}
void	reverb::setwidth(float value) {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    rev->setwidth(value);
}
float	reverb::getwidth() const {
    return rev->getwidth();
}
void	reverb::setmode(float value) {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    rev->setmode(value);
}
float	reverb::getmode() const {
    return rev->getmode();
}