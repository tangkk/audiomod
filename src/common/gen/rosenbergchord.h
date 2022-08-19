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

// this defines chords generated as three or more rosenberg glottal pulse trains

#include <vector>

class rosenberg;

class rosenbergchord {
public:
    rosenbergchord(float sample_rate, float alpha, float beta, const std::vector<float> &freqs);
    ~rosenbergchord();

    void reset();

    // TODO: implement this
    void updateChord(const std::vector<float> &freqs);

    float getNextSample(); 

private:
    std::vector<rosenberg*> rsbchord;
    int chordsize;

};