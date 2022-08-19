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

#include "rosenberg.h"
#include <cmath>

rosenberg::rosenberg(float sample_rate, float freq, float alpha, float beta) {
    m_samplerate = sample_rate;
    m_freq = freq;
    period_in_samples = round(1.f / freq * sample_rate);

    // TODO: value protection (0-1)
    m_alpha = alpha;
    m_beta = beta;
    m_phase = 0;
    phase_in_samples = 0;

    n1 = round(alpha * period_in_samples);
    inv_n1 = 1.f / static_cast<float>(n1);
    n2 = round(beta * period_in_samples);
    inv_2n2 = 0.5 / static_cast<float>(n2);

}

rosenberg::~rosenberg() {}

float rosenberg::getNextSample() {
    float res = 0;
    if (phase_in_samples <= n1) {
        res = 0.5 * (1 - cosf(M_PI * phase_in_samples * inv_n1));
    } else if (phase_in_samples - n1 <= n2) {
        res = cosf(M_PI * (phase_in_samples - n1) * inv_2n2);
    } else {
        res = 0;
    }

    if (++phase_in_samples > period_in_samples) phase_in_samples = 0;
    m_phase = static_cast<float>(phase_in_samples) / static_cast<float>(period_in_samples);

    return res;
}

