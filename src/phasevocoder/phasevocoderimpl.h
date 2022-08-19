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

#include "phasevocoderinterface.h"

#include <algorithm>
#include <set>

#include "../common/base/circularqueue.h"
#include "../common/dsp/windowfunc.h"
#include "../common/dsp/FFT.h"
#include "../common/system/sys.h"


class rosenberg;
class rosenbergchord;

namespace audiomod
{


typedef float data_type;

// class AudioCurveCalculator;
// class StretchCalculator;

class phasevocodercore::Impl
{
public:
    Impl(size_t sampleRate, size_t channels, Options options,
         float initialTimeRatio, float initialPitchScale);
    ~Impl();
    
    // void reset();
    // void setTimeRatio(float ratio);
    // void setPitchScale(float scale);

    float getTimeRatio() const;
    float getPitchScale() const;

    // size_t getLatency() const;

    // void study(const float *const *input, size_t num_samples, bool final);
    void processNormal(const float *const *input, size_t num_samples);

    /*
    process without stretch ratio, normal phase vocoder setting (no "final" argument)
    the code is mostly simplified
    */
    void processConstant(const float *const *input, size_t num_samples);

    enum carrierType {
        rosenbergType = 0,
        rosenbergChordType
    };
    void processVocoder(const float *const *input, size_t num_samples, carrierType ct);

    int numsamples_available() const;
    int numsamples_availableCarrier() const;

    size_t retrieve(float *const *output, size_t num_samples) const;
    size_t retrieveCarrier(float *const *output, size_t num_samples) const;

    // float getFrequencyCutoff(int n) const;
    // void setFrequencyCutoff(int n, float f);

    size_t getInputHopSize() const {
        return m_hopsize;
    }


    // std::vector<float> getPhaseResetCurve() const;
    // std::vector<int> getExactTimePoints() const;

    size_t getChannelCount() const {
        return m_channels;
    }
    
    // void calculateStretch();

    void setDebugLevel(int level);
    static void setDefaultDebugLevel(int level) { m_defaultDebugLevel = level; }
    static void setDefaultFftSize(int fftsize) { m_defaultFftSize = fftsize; }
    static void setDefaultHopSize(int hopsize) { m_defaultHopsize = hopsize; }
    static void setDefaultCoreMode(int coremode) { m_defaultCoreMode = coremode; }

protected:
    size_t m_sampleRate;
    size_t m_channels;

    // void prepareChannelMS(size_t channel, const float *const *inputs,
    //                       size_t offset, size_t num_samples, float *prepared);

    size_t enbufferChannel(size_t channel, const float *const *inputs,
                          size_t offset, size_t num_samples);

    size_t enbufferChannelVocoder(size_t channel, const float *const *inputs,
                          size_t offset, size_t num_samples, carrierType ct);

    // void processSlices(size_t channel, bool &any, bool &last);
    int processOneSlice(); // across all channels, for real time use, return outframes
    int processSliceForChannel(size_t channel, size_t phaseIncrement,
                                size_t shiftIncrement); // return outframes

    int processOneSliceConstant(); // across all channels, for real time use, return outframes
    int processOneSliceVocoder(); // across all channels, for real time use, return outframes
    // int processOneSliceVocoderRealImag(); // across all channels, for real time use, return outframes

    bool inbufReady(size_t channel);
    void calculateIncrements(size_t &phaseIncrement,
                             size_t &shiftIncrement);

    void analyzeSlice(size_t channel);
    // void analyzeSliceRealImag(size_t channel);
    void analyzeCarrierSlice(size_t channel);
    // void analyzeCarrierSliceRealImag(size_t channel);

    void roboticSlice(size_t channel);
    void whisperSlice(size_t channel);

    void modifySliceSimple(size_t channel, size_t phaseIncrement);
    void modifySlicePhaseLocked(size_t channel, size_t phaseIncrement); // phase-locked vocoder
    void modifySliceIntRatio(size_t channel, size_t phaseIncrement); // integer ratio modification

    void modifySliceVocoder(size_t channel, size_t phaseIncrement); // note that we also keep phasereset for future use
    // void modifySliceVocoderRealImag(size_t channel, size_t phaseIncrement);

    void freqCompSlice(size_t channel, float freq_comp);
    void formantShiftSlice(size_t channel, float env_comp);
    void formantPreserveSlice(size_t channel);
    void maleToFemale(size_t channel);
    void femaleToMale(size_t channel);

    void synthesiseSlice(size_t channel, size_t shiftIncrement);
    void synthesiseSliceCarrier(size_t channel, size_t shiftIncrement);
    // void synthesiseSliceCarrierRealImag(size_t channel, size_t shiftIncrement);

    int writeSlice(size_t channel, size_t shiftIncrement, bool last); // return outframes
    int writeSliceCarrier(size_t channel, size_t shiftIncrement);

    void writeOutSimple(circularqueue<float> &to, float *from,
                     size_t qty, size_t &outcnt);

    void calculateSizes();
    void configure();

    float getHopSizeRatio() const;
    bool isIntRatio() const;
    
    size_t nextpow2(size_t value); // to next power of two

    /* tangkk: circular shift of the windowed segment before FFT (and after IFFT) provides zero-phase analysis */
    template <typename T, typename S>
    void fftshift(T *target, int targetSize,
                         S *src, // destructive to src
                         windowfunc<float> *window) {
        const int windowSize = window->GetSize();
        const int hs = targetSize / 2;
        if (windowSize == targetSize) {
            window->apply(src);
            vector_assign(target, src + hs, hs);
            vector_assign(target + hs, src, hs);
        } else {
            printf("fftshift error!\n");
            throw 6;
        }
    }

    template <typename T, typename S>
    void ifftshift(T *target, int targetSize,
                         S *src, // destructive to src
                         windowfunc<float> *window) {

        const int windowSize = window->GetSize();
        const int hs = targetSize / 2;
        if (windowSize == targetSize) {
            vector_assign(target, src + hs, hs);
            vector_assign(target + hs, src, hs);
            window->apply(target);
        } else {
            printf("ifftshift error!\n");
            throw 6;
        }
    }

    // bool resampleBeforeStretching() const {return false;}
    
    float m_timeRatio;
    float m_pitchScale;

    // n.b. either m_fftSize is an integer multiple of m_windowSize,
    // or vice versa
    size_t m_fftSize;
    size_t m_analyzeWindowSize; //!!! or use m_analyzeWindow->GetSize() throughout?  (tangkk: analysis window size)
    size_t m_synthesisWindowSize; //!!! or use m_synthesisWindow->GetSize() throughout? (tangkk: synthesis window size)
    size_t m_hopsize;
    size_t m_outbufSize;

    // size_t m_maxProcessSize;
    // size_t m_expectedInputDuration;

    // bool m_realtime;
    Options m_options;
    int m_logLevel;

    windowfunc<float> *m_analyzeWindow;
    // SincWindow<float> *m_afilter;
    windowfunc<float> *m_synthesisWindow;

    std::vector<rosenberg *> rsb_gen; // rosenberg glottal pulse generator, each per channel
    std::vector<rosenbergchord *> rsb_chordgen;

    class channelinfo; 
    std::vector<channelinfo *> m_audioData;
    std::vector<channelinfo *> m_carrierData;

    static int m_defaultDebugLevel;
    static size_t m_defaultHopsize; // could be reset during construction time
    static size_t m_defaultFftSize; // could be reset during construction time
    static int m_defaultCoreMode;

    // for phase-locked vocoder
    std::vector<int> peak_loc;
    std::vector<int> prev_peak_loc;

    // resample buffer
    float *resamplebuf;
    size_t resamplebufSize;
    void SetResampleBufSize(size_t resamplebufSize);

};

}

