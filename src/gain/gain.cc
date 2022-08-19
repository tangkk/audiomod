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

#include "gain.h"

#include <cmath>


gain::gain(int sampleRate, int num_channels, float initgain) {
    num_channels_ = num_channels;
    sample_rate_ = sampleRate;
    gain_ = initgain;
}


void gain::processBlock (float *const * bufferData, int num_samples) {
    for (int channel = 0; channel < num_channels_; ++channel) {
        for (int j=0; j<num_samples; j++) {
            bufferData[channel][j] *= gain_;
            if (bufferData[channel][j] > 1) bufferData[channel][j] = 1;
            if (bufferData[channel][j] < -1) bufferData[channel][j] = -1;
        }
    }
    
}
