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
#include "limiter.h"

limiter::limiter(int sampleRate, int numChannels, float dBThreshold, float dBMakeUpGain, float attackTimeMs, float releaseTimeMs) :
makeUpGain_(pow(10.0, dBMakeUpGain / 20.0)),
threshold_(pow(10.0, (dBThreshold - LIMIT_OFFSET) / 20.0)), ahead_(6.0),
alphaAttack_(exp(-1.0 / (sampleRate * 0.001 * attackTimeMs))),
alphaRelease_(exp(-1.0 / (sampleRate * 0.001 * releaseTimeMs))) {

    xPeaks_ = new float[numChannels];
    gains_ = new float[numChannels];
    buffers_.resize(numChannels);

    num_channels_ = numChannels;
    sample_rate_ = sampleRate;
    tauAttack_ = attackTimeMs;
    tauRelease_ = releaseTimeMs;
    for(int i = 0; i < num_channels_; i++) {
        xPeaks_[i] = pow(10.0, -120.0 / 20.0);
        gains_[i]  = 1.f;
        buffers_[i].resize((int)(sample_rate_ * 0.001 * ahead_) + 1, 0.f);
    }
}

limiter::~limiter() {
    delete[] xPeaks_;
    delete[] gains_;
}

float limiter::process(float x, int channel) {
    x *= makeUpGain_;
    float x_abs = fabs(x);
    if (x_abs < 0.000001)
        x_abs = 0.000001;
    float alpha = x_abs > xPeaks_[channel] ? alphaAttack_ : alphaRelease_;
    xPeaks_[channel] = alpha * xPeaks_[channel] + (1.0 - alpha) * x_abs;
    float gain = std::fmin(1.f, threshold_ / xPeaks_[channel]);
    alpha = gain < gains_[channel] ? alphaAttack_ : alphaRelease_;
    gains_[channel] = alpha * gains_[channel] + (1.0 - alpha) * gain;
    float y = buffers_[channel].front() * gains_[channel];
    buffers_[channel].pop_front();
    buffers_[channel].push_back(x);

    if (y > 1) y = 1;
    if (y < -1) y = -1;
    
    return y;
}

void limiter::process(float* buffer, int channel, int numSamples) {
    for (int i=0; i<numSamples; i++) {
        buffer[i] = process(buffer[i], channel);
    }
}

void limiter::processBlock(float *const *bufferData, int num_samples) {
    for (int i=0; i<num_channels_; i++) {
        float *channelData = bufferData[i];
        process(channelData, i, num_samples);
    }
}

void limiter::resetBuffer() {
    for(int i = 0; i < num_channels_; i++) {
        buffers_[i].clear();
        buffers_[i].resize((int)(sample_rate_ * 0.001 * ahead_) + 1, 0.f);
    }
}
    

