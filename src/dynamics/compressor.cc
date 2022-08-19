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
#include "compressor.h"

compressor::compressor(int sampleRate, int numChannels, float dBThreshold, float ratio, float dBMakeUpGain, float attackTimeMs, float releaseTimeMs):
    dbThreshold_(dBThreshold),
    dbMakeUpGain_(dBMakeUpGain),
    ratio_(ratio),
    alphaAttack_(exp(-1.0 / (sampleRate * 0.001 * attackTimeMs))),
    alphaRelease_(exp(-1.0 / (sampleRate * 0.001 * releaseTimeMs))) {
    
    dbx_g = new float[numChannels];
    dbx_l = new float[numChannels];
    dby_g = new float[numChannels];
    dby_l = new float[numChannels];
    dbyL_prev = new float[numChannels];
    c = new float[numChannels];

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;
    tauAttack_ = attackTimeMs;
    tauRelease_ = releaseTimeMs;

    alphaAttack_ = exp(-1/(0.001 * sampleRate * tauAttack_));
    alphaRelease_ = exp(-1/(0.001 * sampleRate * tauRelease_));
    for (int i = 0 ; i < numChannels ; ++i) {
        dbx_g[i] = 0;    dby_g[i] = 0;
        dbx_l[i] = 0;    dby_l[i] = 0;
        c[i] = 0; dbyL_prev[i]=0;
    }
}

compressor::~compressor() {
    delete[] dbx_g;
    delete[] dbx_l;
    delete[] dby_g;
    delete[] dby_l;
    delete[] dbyL_prev;
    delete[] c;
}

float compressor::process(float x, int channel) {
    if (fabs(x) < 0.000001)
        dbx_g[channel] = -120;
    else
        dbx_g[channel] = 20*log10(fabs(x));
    //Gain computer- static apply input/output curve
    if (dbx_g[channel] >= dbThreshold_)
        dby_g[channel] = dbThreshold_ + (dbx_g[channel] - dbThreshold_) / ratio_;
    else
        dby_g[channel] = dbx_g[channel];
    dbx_l[channel] = dbx_g[channel] - dby_g[channel];

    //Ballistics- smoothing of the gain
    float alpha = dbx_l[channel] > dbyL_prev[channel] ? alphaAttack_ : alphaRelease_;
    dby_l[channel] = alpha * dbyL_prev[channel] + (1 - alpha) * dbx_l[channel] ;
    
    //find control
    c[channel] = pow(10, (dbMakeUpGain_ - dby_l[channel]) / 20);
    dbyL_prev[channel] = dby_l[channel];
    return x * c[channel];
}

void compressor::process(float* buffer, int channel, int numSamples) {
    for(int i = 0; i < numSamples; i++) {
        buffer[i] = process(buffer[i], channel);
    }
}

void compressor::processBlock(float *const *bufferData, int num_samples) {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    for (int i=0; i<num_channels_; i++) {
        float *channelData = bufferData[i];
        process(channelData, i, num_samples);
    }
}
