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

// this defines the rosenberg glottal pulse generation routine
// refer to https://www.mathworks.com/matlabcentral/fileexchange/45316-rosenburg-glottal-pulse

class rosenberg {

public:
    rosenberg(float sample_rate, float freq, float alpha, float beta);
    ~rosenberg();

    // TODO: implement the following functions
    void updateFreq(float freq);
    void updateSampleRate(float sample_rate);
    void updatePhase(float phase);
    void updateAlpha(float alpha);
    void updateBeta(float beta);

    void reset();

    float getNextSample(); // generate the next sample according to the current state

private:
    float m_samplerate; // the digital sample rate of the generator
    float m_freq; // the frequency of the glottal pulse train, 1/freq should be the period
    int period_in_samples; // 1/freq*samplerate

    // glotal opening and closing raio, the rest (1 - alpha - beta) should be silence
    float m_alpha; // glottal opening ratio, should be 0-0.5
    float m_beta; // glottal closing ratio, should be 0-0.5
    int n1; // number of samples in opening phase
    float inv_n1;
    int n2; // number of samples in closing phase
    float inv_2n2;

    // phase
    float m_phase; // the current phase of generation (0-1)
    int phase_in_samples;


};