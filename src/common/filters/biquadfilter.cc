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

/*
  ==============================================================================
    reference: https://www.w3.org/2011/audio/audio-eq-cookbook.html
    All filter transfer functions were derived from analog prototypes 
    (that are shown below for each equalizer (EQ) filter type) 
    and had been digitized using the Bilinear Transform (BLT).
  ==============================================================================
*/

#include "biquadfilter.h"
#include <cmath>
#include <stdio.h>

biquadfilter::biquadfilter(int sampleRate, int numChannels, Type type, float cutoffFreq, float q, float dBGain)
:type_(type), cutoffFreq_(cutoffFreq), sampleRate_(sampleRate), q_(q), dBGain_(dBGain) {
    x1s = new float[numChannels];
    x2s = new float[numChannels];
    y1s = new float[numChannels];
    y2s = new float[numChannels];
    for(int i = 0; i < numChannels; i++)
        x1s[i] = x2s[i] = y1s[i] = y2s[i] = 0.0;
    computeCoeffs();

    sample_rate_ = sampleRate;
    num_channels_ = numChannels;

    outready_ = true;
}

biquadfilter::~biquadfilter() {
    delete[] x1s;
    delete[] x2s;
    delete[] y1s;
    delete[] y2s;
}

float biquadfilter::process(float x, int channel) {
    float y = (b0 * x + b1 * x1s[channel] + b2 * x2s[channel] - a1 * y1s[channel] - a2 * y2s[channel]) / a0;
    // printf("%f\n",x);
    // printf("%f, %f, %f, %f, %f, %f\n", b0, b1, b2, a0, a1, a2);
    x2s[channel] = x1s[channel];
    y2s[channel] = y1s[channel];
    x1s[channel] = x;
    y1s[channel] = y;
    // printf("%f\n",y);
    return y;
}

void biquadfilter::processBlock (float *const * bufferData, int num_samples) {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    for (int i=0; i<num_channels_; i++) {
        for (int j=0; j<num_samples; j++) {
            bufferData[i][j] = process(bufferData[i][j],i);
        }
    }
}

void biquadfilter::setCutoff(float cutoffFreq) {
    cutoffFreq_ = cutoffFreq;
    computeCoeffs();
}

void biquadfilter::setQ(float q) {
    q_ = q;
    computeCoeffs();
}

void biquadfilter::setGain(float dBGain) {
    dBGain_ = dBGain;
    computeCoeffs();
}

float biquadfilter::getCutoff() {
    return cutoffFreq_;
}

float biquadfilter::getQ() {
    return q_;
}

float biquadfilter::getGain() {
    return dBGain_;
}

void biquadfilter::setType(Type type) {
    type_ = type;
    computeCoeffs();
}

void biquadfilter::setAllParams(Type type, float cutoffFreq, float q, float dBGain) {
    type_ = type;
    cutoffFreq_ = cutoffFreq;
    q_ = q;
    dBGain_ = dBGain;
    computeCoeffs();
}

void biquadfilter::computeCoeffs() {
    std::unique_lock<std::mutex> lock(coffupdate_mutex);
    float a = pow(10.0, dBGain_ / 40.0);
    float omega = 2 * M_PI * cutoffFreq_ / sampleRate_;
    float alpha = sin(omega) / 2.0 / q_;
    
    switch (type_) {
        case lowPass: // H(s) = 1 / (s^2 + s/Q + 1)
            b0 = b2 = (1.0 - cos(omega)) / 2.0;
            b1 = 1.0 - cos(omega);
            a0 = 1.0 + alpha;
            a1 = -2.0 * cos(omega);
            a2 = 1 - alpha;
            break;
        case highPass: // H(s) = s^2 / (s^2 + s/Q + 1)
            b0 = b2 = (1.0 + cos(omega)) / 2.0;
            b1 = -(1.0 + cos(omega));
            a0 = 1.0 + alpha;
            a1 = -2.0 * cos(omega);
            a2 = 1 - alpha;
            break;
        case bandpass_constant_skirt: // H(s) = s / (s^2 + s/Q + 1)
            b0 = sin(omega) / 2;
            b1 = 0;
            b2 = - sin(omega) / 2;
            a0 = 1 + alpha;
            a1 = -2 * cos(omega);
            a2 = 1 - alpha;
            break;
        case bandpass_constant_zero: // H(s) = (s/Q) / (s^2 + s/Q + 1)
            b0 = alpha;
            b1 = 0;
            b2 = -alpha;
            a0 = 1 + alpha;
            a1 = -2 * cos(omega);
            a2 = 1 - alpha;
            break;
        case notch: // H(s) = (s^2 + 1) / (s^2 + s/Q + 1)
            b0 = 1;
            b1 = -2 * cos(omega);
            b2 = 1;
            a0 = 1 + alpha;
            a1 = -2 * cos(omega);
            a2 = 1 - alpha;
            break;
        case allpass: // H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)
            b0 = 1 - alpha;
            b1 = -2 * cos(omega);
            b2 = 1 + alpha;
            a0 = 1 + alpha;
            a1 = -2 * cos(omega);
            a2 = 1 - alpha;
            break;
        case peaking: // H(s) = (s^2 + s*A/Q + 1) / (s^2 + s/(AQ) + 1)
            b0 = 1 + alpha * a;
            b1 = -2 * cos(omega);
            b2 = 1 - alpha * a;
            a0 = 1 + alpha / a;
            a1 = -2 * cos(omega);
            a2 = 1 - alpha / a;
            break;
        case lowShelf: // H(s) = A * ((s^2 + sqrt(A)/Q * s + A) / (A*s^2 + sqrt(A)/Q + 1))
            b0 = a * (a + 1 - (a - 1) * cos(omega) + 2 * sqrt(a) * alpha);
            b1 = 2 * a * (a - 1 - (a + 1) * cos(omega));
            b2 = a * (a + 1 - (a - 1) * cos(omega) - 2 * sqrt(a) * alpha);
            a0 = a + 1  + (a - 1) * cos(omega) + 2 * sqrt(a) * alpha;
            a1 = -2 * (a - 1 + (a + 1) * cos(omega));
            a2 = a + 1  + (a - 1) * cos(omega) - 2 * sqrt(a) * alpha;
            break;
        case highShelf: // H(s) = A * ((A*s^2 + sqrt(A)/Q * s + 1) / (s^2 + sqrt(A)/Q + A))
            b0 = a * (a + 1 + (a - 1) * cos(omega) + 2 * sqrt(a) * alpha);
            b1 = -2 * a * (a - 1 + (a + 1) * cos(omega));
            b2 = a * (a + 1 + (a - 1) * cos(omega) - 2 * sqrt(a) * alpha);
            a0 = a + 1  - (a - 1) * cos(omega) + 2 * sqrt(a) * alpha;
            a1 = 2 * (a - 1 - (a + 1) * cos(omega));
            a2 = a + 1  - (a - 1) * cos(omega) - 2 * sqrt(a) * alpha;
            break;

        default:
            b0 = b1 = b2 = a0 = a1 = a2 = 0.0;
            break;
    }
}
