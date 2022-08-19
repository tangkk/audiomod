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

#include "rosenbergchord.h"
#include "rosenberg.h"

rosenbergchord::rosenbergchord(float sample_rate, float alpha, float beta, const std::vector<float> &freqs) {

    for (auto &freq : freqs) {
        rsbchord.emplace_back(new rosenberg(sample_rate, freq, alpha, beta));
    }
    chordsize = rsbchord.size();

}

rosenbergchord::~rosenbergchord() {
    for (auto &rsb : rsbchord) {
        if (rsb != nullptr) {
            delete rsb;
            rsb = nullptr;
        }
    }
    rsbchord.clear();
}

float rosenbergchord::getNextSample() {
    float res = 0;
    for (auto &rsb : rsbchord) {
        res += rsb->getNextSample() / chordsize;
    }
    return res;
}