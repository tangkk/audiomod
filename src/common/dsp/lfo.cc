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

#include "lfo.h"

#include <cmath>


float lfo_neg2one(float phase, int waveform) {
    switch(waveform) {
        case kWaveformTriangle:
            if(phase < 0.25f)
                return 4.0f*phase;
            else if(phase < 0.75f)
                return 1.0f - 4.0f*(phase - 0.25f);
            else
                return -1.0f + 4.0f*(phase - 0.75f);
        case kWaveformSquare:
            if(phase < 0.5f)
                return 1.0f;
            else
                return -1.0f;
        case kWaveformSawtooth:
            if(phase < 0.5f)
                return 2.0f*phase;
            else
                return 2.0f*phase - 2.0f;
        case kWaveformInverseSawtooth:
            if(phase < 0.5f)
                return -2.0f*phase;
            else
                return 2.0f - 2.0f*phase;
        case kWaveformSquareSlopedEdges:
            if(phase < 0.48f)
                return 1.0f;
            else if(phase < 0.5f)
                return 1.0f - 50.0f*(phase - 0.48f);
            else if(phase < 0.98f)
                return -1.0f;
            else
                return 50.0f*(phase - 0.98f) - 1.0f;
        case kWaveformSine:
        default:
            return sinf(2.0 * M_PI * phase);
    }
}

float lfo_zero2one(float phase, int waveform) {
    switch(waveform) {
        case kWaveformTriangle:
            if(phase < 0.25f)
                return 0.5f + 2.0f*phase;
            else if(phase < 0.75f)
                return 1.0f - 2.0f*(phase - 0.25f);
            else
                return 2.0f*(phase-0.75f);
        case kWaveformSquare:
            if(phase < 0.5f)
                return 1.0f;
            else
                return 0.0f;
        case kWaveformSawtooth:
            if(phase < 0.5f)
                return 0.5f + phase;
            else
                return phase - 0.5f;
        case kWaveformInverseSawtooth:
            if (phase < 0.5f)
                return 0.5f - phase;
            else
                return 1.5f - phase;
        case kWaveformSquareSlopedEdges:
            if(phase < 0.48f)
                return 1.0f;
            else if(phase < 0.5f)
                return 1.0f - 50.0f*(phase - 0.48f);
            else if(phase < 0.98f)
                return 0.0f;
            else
                return 50.0f*(phase - 0.98f);
        case kWaveformSine:
        default:
            return 0.5f + 0.5f*sinf(2.0 * M_PI * phase);
    }
}