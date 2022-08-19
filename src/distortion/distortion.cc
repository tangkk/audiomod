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
#include "distortion.h"

#include <cmath>


distortion::distortion(int sampleRate, int num_channels):
                    outputGain(0.1),
                    drive(0.5),
                    typeDistortion(0) {
    num_channels_ = num_channels;
    sample_rate_ = sampleRate;
    threshold = sample_rate_;
}


void distortion::updateParameters(int sampleRate, int num_channels, int newBlock_size,
           float newOutputGain, int newThreshold, float newDrive, int newTypeDistortion) {
    num_channels_ = num_channels;
    sample_rate_ = sampleRate;
    outputGain = newOutputGain;
    threshold = newThreshold;
    drive = newDrive;
    typeDistortion = newTypeDistortion;

}


void distortion::processBlock (float *const * bufferData, int num_samples) {
    for (int channel = 0; channel < num_channels_; ++channel) {
        float* channelData = bufferData[channel];
        
        float out;
        
        for (int sample = 0; sample < num_samples; ++sample) {
            const float in = channelData[sample] * drive;
            switch (typeDistortion) {
                case distortionTypeHardClipping: {
                    
                    float thr = pow(10, -1*threshold/20);
                    if (in > thr)
                        out = thr;
                    else if (in < -thr)
                        out = -thr;
                    else
                        out = in;
                    break;
                }
                case distortionTypeSoftClipping: {
                    float thr1 = pow(10, -1*threshold/20)/2;
                    float thr2 = pow(10, -1*threshold/20);
                    if (in > thr2)
                        out = 1.0f;
                    else if (in > thr1)
                        out = 1.0f - powf (2.0f - 1 / thr1 * in, 2.0f) * (1 - 2 * thr1);
                    else if (in < -thr2)
                        out = -1.0f;
                    else if (in < -thr1)
                        out = -1.0f + powf (2.0f + 1 / thr1 * in, 2.0f) * (1 - 2 * thr1);
                    else
                        out = 2.0f * in;
                    out *= 0.5f;
                    break;
                }
                case distortionTypeExponential: {
                    if (in > 0.0f)
                        out = 1.0f - expf (-in);
                    else
                        out = -1.0f + expf (in);
                    break;
                }
                case distortionTypeFullWaveRectifier: {
                    out = fabsf (in);
                    break;
                }
                case distortionTypeHalfWaveRectifier: {
                    if (in > 0.0f)
                        out = in;
                    else
                        out = 0.0f;
                    break;
                }
            }
            channelData[sample] = out * outputGain;
        }
    }
    
}
