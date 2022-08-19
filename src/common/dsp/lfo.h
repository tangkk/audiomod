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

#pragma once

#include <cmath>

enum BasicWaveforms {
    kWaveformSine = 1,
    kWaveformTriangle,
    kWaveformSquare,
    kWaveformSawtooth,
    kWaveformInverseSawtooth,
    kWaveformSquareSlopedEdges,
    kNumWaveforms
};

float lfo_neg2one(float phase, int waveform);
float lfo_zero2one(float phase, int waveform);