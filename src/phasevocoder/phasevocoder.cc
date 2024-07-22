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

#include "phasevocoder.h"
#include "phasevocoderinterface.h"

#include <cmath>
//#include <sys/time.h>

namespace audiomod {

phasevocoder::phasevocoder(int sampleRate, int numChannels, float timeratio, float pitchshift, int mode, int coremode, int fftsize, int hopsize) {
    timeratio_ = timeratio;
    pitchscale_ = pitchshift != 0 ? std::pow(2.0, pitchshift / 12) : 1.0;

    options_ = 0;

    // mode selection
    if (mode == GENDER_CHANGE) {
        options_ |= phasevocodercore::OptionGenderChange;
    }
    if (mode == FORMANT_PRESERVE) {
        options_ |= phasevocodercore::OptionFormantPreserved;
    }
    if (mode == ROBOTIC) {
        options_ |= phasevocodercore::OptionRobotic;
    }
    if (mode == WHISPER) {
        options_ |= phasevocodercore::OptionWhisper;
    }


    sample_rate_ = sampleRate;
    num_channels_ = numChannels;
    // defaultfftsize_ = 2048; // hard coded
    defaultfftsize_ = fftsize;
    // defaulthopsize_ = 256; // hard coded
    defaulthopsize_ = hopsize;
    defaultcoremode_ = coremode;

    m_mode = mode;
    m_log = 1; // open debug message

    outready_ = false;

    init();
    
}

phasevocoder::~phasevocoder() {
    if (ts != nullptr) {
        delete ts;
        ts = nullptr;
    }
}

void phasevocoder::init() {
    phasevocodercore::setDefaultDebugLevel(m_log); // better set before instantiating objects
    phasevocodercore::setDefaultFftSize(defaultfftsize_);
    phasevocodercore::setDefaultHopSize(defaulthopsize_);
    phasevocodercore::setDefaultCoreMode(defaultcoremode_);

    if (defaultcoremode_ == 0) {
        printf("use simple phase-vocoder\n");
    } else if (defaultcoremode_ == 1) {
        printf("use phase-locked phase-vocoder\n");
    } else if (defaultcoremode_ == 2) {
        printf("use int-ratio phase-vocoder\n");
    } else {
        printf("use simple phase-vocoder\n");
    }
    ts = new phasevocodercore(sample_rate_, num_channels_, options_, timeratio_, pitchscale_);
}

void phasevocoder::processInData(float *const * inData, int num_in_samples) {

    int numres = 0;
    if (m_mode == NORMAL_STRETCH || m_mode == NORMAL_SHIFT || m_mode == GENDER_CHANGE || m_mode == FORMANT_PRESERVE || m_mode == ROBOTIC || m_mode == WHISPER) {
        ts->processNormal(inData, num_in_samples);
        numres = ts->numsamples_available();
    } else if (m_mode == VOCODER_ROSENBERG) {
        ts->processVocoder(inData, num_in_samples, 0);
        numres = ts->numsamples_availableCarrier();
    } else if (m_mode == VOCODER_CHORD) {
        ts->processVocoder(inData, num_in_samples, 1);
        numres = ts->numsamples_availableCarrier();
    } else if (m_mode == CONSTANT) {
        ts->processConstant(inData, num_in_samples);
        numres = ts->numsamples_available();
    } else {
        // do nothing here
    }

    
    num_res_ = numres;
}

void phasevocoder::getOutData(float *const * outData, int num_out_samples) {
    if (num_out_samples > num_res_) {
        num_out_samples = num_res_;
    }
    if (m_mode == CONSTANT || m_mode == NORMAL_STRETCH || m_mode == NORMAL_SHIFT || m_mode == GENDER_CHANGE || m_mode == FORMANT_PRESERVE || m_mode == ROBOTIC || m_mode == WHISPER) {
        ts->retrieve(outData, num_out_samples);
    } else if (m_mode == VOCODER_ROSENBERG || m_mode == VOCODER_CHORD) {
        ts->retrieveCarrier(outData, num_out_samples);
    } else {
        // do nothing here
    }

    outready_ = true; // offline mode, always ready

}

void phasevocoder::processBlock(float *const * bufferData, int num_samples) {

    static int debug_cnt = 0;
    if (debug_cnt == 0) {
        printf("num_samples:%d\n", num_samples);
    }

    int ret = 0;
    if (m_mode == NORMAL_SHIFT || m_mode == GENDER_CHANGE || m_mode == FORMANT_PRESERVE || m_mode == ROBOTIC || m_mode == WHISPER) {
        ret = processBlockNormal(bufferData, num_samples);
    } else if (m_mode == VOCODER_ROSENBERG) {
        ret = processBlockVocoder(bufferData, num_samples, 0);
    } else if (m_mode == VOCODER_CHORD) {
        ret = processBlockVocoder(bufferData, num_samples, 1);
    } else if (m_mode == CONSTANT) {
        ret = processBlockConstant(bufferData, num_samples);
    } else {
        // do nothing here
    }

    if (ret < 0) {
        printf("out of output data, skip this block! debug_cnt:%d, numres:%d\n", debug_cnt, num_res_);
        outready_ = false;
    } else {
        outready_ = true;
    }
    
    debug_cnt++;
}

int phasevocoder::processBlockNormal(float *const * bufferData, int num_samples) {

    static size_t sum_out = 0;
    static size_t sum_in = 0;

    int count = num_samples;
    sum_in += count;


    // printf("******************\n");
    // printf("bufferData_in: %d\n ", count);
    ts->processNormal(bufferData, count);
    int numres = ts->numsamples_available();
    // printf("availFinal: %d\n",numres);

    num_res_ = numres;

    // FIXME: this may lead to more and more delay in realtime application
    // if (numres >= count && outbuf_ready) {
    if (numres >= count) {
        ts->retrieve(bufferData, count);
        // printf("retrieved: %d\n",count);
        sum_out += count;
        return 0;
    } else {
        return -1;
    }
}

int phasevocoder::processBlockConstant(float *const * bufferData, int num_samples) {
    static size_t sum_out = 0; // TEMP var
    static size_t sum_in = 0;

    int count = num_samples;
    sum_in += count;

    // printf("******************\n");
    // printf("bufferData_in: %d, ", count);
    ts->processConstant(bufferData, count);
    int numres = ts->numsamples_available();
    // printf("availFinal: %d\n",numres);

    num_res_ = numres;

    // FIXME: this may lead to more and more delay in realtime application
    // if (numres >= count && outbuf_ready) {
    if (numres >= count) {
        ts->retrieve(bufferData, count);
        // printf("retrieved: %d\n",count);
        sum_out += count;
        return 0;
    } else {
        return -1;
    }
}

int phasevocoder::processBlockVocoder(float *const * bufferData, int num_samples, int carrierType) {
    static size_t sum_out = 0; // TEMP var
    static size_t sum_in = 0;

    int count = num_samples;
    sum_in += count;

    // printf("******************\n");
    // printf("bufferData_in: %d, ", count);
    ts->processVocoder(bufferData, count, carrierType);
    int numres = ts->numsamples_availableCarrier();
    // printf("availFinal: %d\n",numres);

    num_res_ = numres;

    // FIXME: this may lead to more and more delay in realtime application
    // if (numres >= count && outbuf_ready) {
    if (numres >= count) {
        ts->retrieveCarrier(bufferData, count);
        // printf("retrieved: %d\n",count);
        sum_out += count;
        return 0;
    } else {
        return -1;
    }
}
} // end of namespace audiomod
