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


#include "modbase.h"

// mode define
#define CONSTANT -1
#define NORMAL_SHIFT 0
#define GENDER_CHANGE 1
#define FORMANT_PRESERVE 2
#define VOCODER_ROSENBERG 3
#define VOCODER_CHORD 4
#define NORMAL_STRETCH 5
#define ROBOTIC 6
#define WHISPER 7

// #define OUTBUF_DELAY 1024

#define NORMAL_PV 0
#define PHASE_LOCKED 1 
#define INT_RATIO 2

namespace audiomod {

class phasevocodercore;

class phasevocoder : public modbase, public modbase_offline {
public:

  /**
   * construct a phasevocoder structure
   * 
   * @param timeratio time-stretch ratio. i.e. outputlen / inputlen = timeratio
   * @param pitchshift pitch-shift amount (in semitone)
   * @param mode 0: pitch-shift, 1: gender change, 2: formant preserve, 3: rosenberg vocoder, 4: chord vocoder, 5, time-stretch
   * @param coremode 0: normal phase vocoder, 1: phase-locked vocoder, 2: int-ratio phase vocoder
   */
  // phasevocoder(int sampleRate, int numChannels, float timeratio, float pitchshift, int mode = NORMAL_SHIFT, int coremode = PHASE_LOCKED, int fftsize = 2048, int hopsize = 256); // force hopsize
  phasevocoder(int sampleRate, int numChannels, float timeratio, float pitchshift, int mode = NORMAL_SHIFT, int coremode = PHASE_LOCKED, int fftsize = 2048, int hopsize = 0); // auto hopsize

  ~phasevocoder();

  /**
   * set mod params
   * @param params the key-val params to be set
   */ 
  void setParams(std::map<std::string, float> params) {

  }

  /**
   * get mod params
   * @param params the key-val params to be returned
   */ 
  void getParams(std::map<std::string, float> &params) {

  }

  void processBlock (float *const * bufferData, int num_samples); // for real-time process
  
  void processInData (float *const * inData, int num_in_samples); // for non-real-time process

  void getOutData(float *const * outData, int num_out_samples);

  bool outputReady() {return outready_;}

private:

  // copy and assignment protection
  phasevocoder(const phasevocoder&) = delete;
  void operator=(const phasevocoder&) = delete;

  // on init
  void init();
  
  // private process methods

  int processBlockNormal(float *const * bufferData, int num_samples);

  int processBlockConstant(float *const * bufferData, int num_samples);

  int processBlockVocoder(float *const *bufferData, int num_samples, int carrierType);
  
  phasevocodercore *ts;
  int options_; // phasevocodercore::Options

  int m_mode;
  int m_log;

  float timeratio_;
  float pitchscale_;

  int defaultfftsize_;
  int defaulthopsize_;
  int defaultcoremode_;

  int sample_rate_;
  int num_channels_;

  bool outready_;

};


} // end of namespace audiomod