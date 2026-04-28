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

#include "f0tonote.h"

namespace audiomod {

f0tonote::f0tonote(int sampleRate, int numChannels, int blockSize, int stepSize, float onsetSens) {
    // note that pyin only accept mono input, if numChannels > 1, take left only
    mn = new MonoNote();
    m_inputSampleRate = sampleRate;
    m_stepSize = stepSize;
    m_pruneThresh = 0.1;
    m_onsetSensitivity = onsetSens;
    m_frameDur = (float) m_stepSize / m_inputSampleRate;
}

f0tonote::~f0tonote() {
    if (mn != nullptr) {
        delete mn;
        mn = nullptr;
    }
}

void f0tonote::processInData (float *const * inData, int num_in_samples) {
    // std::cerr << "f0tonote processInData..." << num_in_samples << std::endl;
    std::vector<float> mpOut(inData[0], inData[0] + num_in_samples); // assuming this is only 1-d (because the input is actually an f0 vector)
    std::vector<float> m_level(inData[1], inData[1] + num_in_samples);

    // std::cerr << "f0tonote 1..." << num_in_samples << std::endl;
    std::vector<std::vector<std::pair<float, float> > > smoothedPitch;
    for (size_t iFrame = 0; iFrame < mpOut.size(); ++iFrame) {
        std::vector<std::pair<float, float> > temp;
        if (mpOut[iFrame] > 0)
        {
            float tempPitch = 12 * std::log(mpOut[iFrame]/440)/std::log(2.) + 69;
            temp.push_back(std::pair<float,float>(tempPitch, .9));
        }
        smoothedPitch.push_back(temp); // if 
        // std::cerr << "temp:" << temp[0].first << temp[0].second << std::endl;
    }
    // std::cerr << "f0tonote 2..." << num_in_samples << std::endl;
    vector<MonoNote::FrameOutput> mnOut = mn->process(smoothedPitch);
    // std::cerr << "f0tonote mnOut.size()..." << mnOut.size() << std::endl;
    // std::cerr << "f0tonote smoothedPitch.size()..." << smoothedPitch.size() << std::endl;

    int onsetFrame = 0;
    bool isVoiced = 0;
    bool oldIsVoiced = 0;
    // size_t nFrame = m_pitchProb.size();
    size_t nFrame = num_in_samples; // this is not samples, but actually frames

    float minNoteFrames = (m_inputSampleRate*m_pruneThresh) / m_stepSize;
    
    std::vector<float> notePitchTrack; // collects pitches for one note at a time
    std::vector<float> values;
    // std::cerr << "nFrame..." << nFrame << std::endl;
    // std::cerr << "minNoteFrames..." << minNoteFrames << std::endl;

    for (size_t iFrame = 0; iFrame < nFrame; ++iFrame) {
        // NOTE: we omit the m_level condition here, which is the yo.rms condition
        // TODO: let's add back the m_level v.s. onset_sensitivity detection condition here (iFrame >= nFrame-2) && ((m_level[iFrame]/m_level[iFrame+2]) > m_onsetSensitivity)
        // isVoiced = mnOut[iFrame].noteState < 3 && smoothedPitch[iFrame].size() > 0;

        isVoiced = mnOut[iFrame].noteState < 3
                   && smoothedPitch[iFrame].size() > 0
                   && (iFrame >= nFrame-2 || ((m_level[iFrame]/m_level[iFrame+2]) > m_onsetSensitivity)); // orig code, maybe 0.8 for sensitivity?
                //    && true;

        // std::cerr << "isVoiced..." << iFrame << "," << isVoiced << "," << mnOut[iFrame].noteState << std::endl;

		if (isVoiced && iFrame != nFrame-1)
        {
            if (oldIsVoiced == 0) // beginning of a note
            {
                onsetFrame = iFrame;
            }
            float pitch = smoothedPitch[iFrame][0].first;
            notePitchTrack.push_back(pitch); // add to the note's pitch track
        } else { // not currently voiced
            if (oldIsVoiced == 1) // end of note
            {
                // std::cerr << notePitchTrack.size() << " " << minNoteFrames << std::endl;
                if (notePitchTrack.size() >= minNoteFrames)
                {
                    std::sort(notePitchTrack.begin(), notePitchTrack.end());
                    float medianPitch = notePitchTrack[notePitchTrack.size()/2];
                    float medianFreq = std::pow(2,(medianPitch - 69) / 12) * 440;
                    values.clear();
                    values.push_back(onsetFrame * m_frameDur);
                    values.push_back((iFrame - onsetFrame) * m_frameDur);
                    // values.push_back(medianFreq);
                    values.push_back(medianPitch);
                    // timestamp = m_timestamp[onsetFrame];
                    // duration = m_timestamp[iFrame] - m_timestamp[onsetFrame];
					//printf("%d, %d, %d\n", m_oNotes, f.timestamp.nsec, f.duration.nsec);
                    // fs[m_oNotes].push_back(f);
                    pYINNoteSequence.push_back(values);
                }
                notePitchTrack.clear();
            }
        }
        oldIsVoiced = isVoiced;
	}
    // std::cerr << "f0tonote 4..." << num_in_samples << std::endl;
    numNotes = pYINNoteSequence.size();
    // std::cerr << "f0tonote numNotes..." << numNotes << std::endl;
}

void f0tonote::getOutData(float *const * outData, int num_out_symbols, std::vector<std::string> *labels) {

    int i = 0;
    for (const auto &x : pYINNoteSequence) {
        // std::cerr << x[0] << "," << x[1] << "," << x[2] << std::endl;
        outData[0][i++] = x[0];
        outData[0][i++] = x[1];
        outData[0][i++] = x[2];
    }
    pYINNoteSequence.clear();
}

float f0tonote::getScalarMeasurement() const {
    return numNotes;
}

}