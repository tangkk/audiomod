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

#include "phasevocoderimpl.h"

#ifndef _WIN32
#include <alloca.h>
#endif

#include <cassert>
#include <cmath>
#include <set>
#include <map>
#include <algorithm>
#include <vector>

#include "channelinfo.h"

#include "../common/dsp/resampler.h"
#include "../common/gen/rosenberg.h"
#include "../common/gen/rosenbergchord.h"


using std::cerr;
using std::endl;
using std::vector;
using std::map;
using std::set;
using std::max;
using std::min;

namespace audiomod {

size_t
// phasevocodercore::Impl::m_defaultHopsize = 512;
phasevocodercore::Impl::m_defaultHopsize = 256;
// phasevocodercore::Impl::m_defaultHopsize = 128;
// phasevocodercore::Impl::m_defaultHopsize = 64;
// phasevocodercore::Impl::m_defaultHopsize = 0; // calculate good increments automatically

size_t
phasevocodercore::Impl::m_defaultFftSize = 2048;
// phasevocodercore::Impl::m_defaultFftSize = 1024; // make it smaller
// phasevocodercore::Impl::m_defaultFftSize = 512;

int
phasevocodercore::Impl::m_defaultCoreMode = 1; // 0: normal phase vocoder; 1: phase-locked vocoder

int
phasevocodercore::Impl::m_defaultDebugLevel = 0;

static bool system_specific_initialised = false;

phasevocodercore::Impl::Impl(size_t sampleRate,
                                size_t channels,
                                Options options,
                                float initialTimeRatio,
                                float initialPitchScale) :
    m_sampleRate(sampleRate),
    m_channels(channels),
    m_timeRatio(initialTimeRatio),
    m_pitchScale(initialPitchScale),
    m_fftSize(m_defaultFftSize),
    m_analyzeWindowSize(m_defaultFftSize),
    m_synthesisWindowSize(m_defaultFftSize),
    m_hopsize(m_defaultHopsize),
    m_outbufSize(m_defaultFftSize * 2),

    m_options(options),
    m_logLevel(m_defaultDebugLevel),
    m_analyzeWindow(0),
    m_synthesisWindow(0)

    {

    if (!system_specific_initialised) {
        system_specific_initialise();
        system_specific_initialised = true;
    }

    configure(); // hyperparameters are computed and configured here
}

phasevocodercore::Impl::~Impl() {

    // FIXME: no delete protection?

    for (size_t chan = 0; chan < m_channels; ++chan) {
        delete m_audioData[chan];
        m_audioData[chan] = nullptr;
        delete m_carrierData[chan];
        m_audioData[chan] = nullptr;
    }
    m_audioData.clear();
    m_carrierData.clear();

    if (m_analyzeWindow != nullptr) {
        delete m_analyzeWindow;
        m_analyzeWindow = nullptr;
    }
    if (m_synthesisWindow != nullptr) {
        delete m_synthesisWindow;
        m_synthesisWindow = nullptr;
    }

    for (size_t chan = 0; chan < m_channels; ++chan) {
        delete rsb_gen[chan];
        rsb_gen[chan] = nullptr;
        delete rsb_chordgen[chan];
        rsb_chordgen[chan] = nullptr;
    }
    rsb_gen.clear();
    rsb_chordgen.clear();

    if (resamplebuf != nullptr) {
        mem_deallocate(resamplebuf);
        resamplebuf = nullptr;
    }
}


float
phasevocodercore::Impl::getTimeRatio() const {
    return m_timeRatio;
}

float
phasevocodercore::Impl::getPitchScale() const {
    return m_pitchScale;
}

float
phasevocodercore::Impl::getHopSizeRatio() const {
    return m_timeRatio * m_pitchScale;
}

bool
phasevocodercore::Impl::isIntRatio() const {
    float efr = getHopSizeRatio();
    if (abs(efr - floorf(efr)) <= 0.001) {
        return true;
    } else {
        return false;
    }
}

size_t
phasevocodercore::Impl::nextpow2(size_t value)
{
    if (!(value & (value - 1))) return value;
    int bits = 0;
    while (value) { ++bits; value >>= 1; }
    value = 1 << bits;
    return value;
}

void
phasevocodercore::Impl::calculateSizes()
{
    // size_t inputHopSize = m_defaultHopsize;
    size_t windowSize = m_fftSize;
    size_t inputHopSize = 0;
    size_t outputHopSize = 0;

    size_t pow2size = nextpow2(windowSize);
    if (windowSize != pow2size) {
        std::cerr << "WARNING: Window size supplied not power of 2! We force it to be one." << std::endl;
        windowSize = pow2size;
    }

    if (m_pitchScale <= 0.0) {
        // This special case is likelier than one might hope, because
        // of naive initialisations in programs that set it from a
        // variable
        std::cerr << "phasevocodercore: WARNING: Pitch scale must be greater than zero!\nResetting it from " << m_pitchScale << " to the default of 1.0: no pitch change will occur" << std::endl;
        m_pitchScale = 1.0;
    }
    if (m_timeRatio <= 0.0) {
        // Likewise
        std::cerr << "phasevocodercore: WARNING: Time ratio must be greater than zero!\nResetting it from " << m_timeRatio << " to the default of 1.0: no time stretch will occur" << std::endl;
        m_timeRatio = 1.0;
    }

    float hsratio = getHopSizeRatio();

    /* 
    NOTE THAT the m_defaultHopsize is actually "useless"
    because in all the cases below, inputHopSize and outputHopSize are set by 
    windowSize (fftsize) and windowIncrRatio, which is determined by "r" and "rsb"
    */

    if (m_defaultHopsize > 0) {
        cerr << "phasevocoder: force hopsize" << endl;
        inputHopSize = m_defaultHopsize;
        outputHopSize = int(floor(inputHopSize * hsratio));

    } else {
        cerr << "phasevocoder: auto hopsize" << endl;
        float windowIncrRatio = 4.5;
        if (hsratio < 1) { 
            if (hsratio == 1.0) windowIncrRatio = 4;
            else if (m_pitchScale < 1.0) windowIncrRatio = 4.5;
            else windowIncrRatio = 6;

            inputHopSize = int(windowSize / windowIncrRatio);
            outputHopSize = int(inputHopSize * hsratio);
        } else {
            if (hsratio == 1.0) windowIncrRatio = 4;
            else windowIncrRatio = 8;

            outputHopSize = int(windowSize / windowIncrRatio);
            inputHopSize = int(outputHopSize / hsratio);
        }
    }

    // windowSize normally won't be changed in the above process

    m_fftSize = windowSize;
    m_analyzeWindowSize = windowSize;
    m_synthesisWindowSize = windowSize;

    /* 
    NOTE here m_hopsize is assigned by inputHopSize, which is decided by window size and ratio
    NOT assigned by m_defaultHopsize
    m_defaultHopsize is actually useless
    */
    m_hopsize = inputHopSize;

    // When squashing, the greatest theoretically possible output
    // increment is the input increment.  When stretching adaptively
    // the sky's the limit in principle, but we expect
    // StretchCalculator to restrict itself to using no more than
    // twice the basic output increment (i.e. input increment times
    // ratio) for any slice.

    if (m_logLevel > 0) {
        cerr << "phasevocoder: sample rate = " << m_sampleRate << ", num channels = " << m_channels << endl;
        cerr << "phasevocoder: time ratio = " << m_timeRatio << ", pitch scale = " << m_pitchScale << ", effective ratio = " << getHopSizeRatio() << endl;
        cerr << "phasevocoder: analysis window size = " << m_analyzeWindowSize << ", synthesis window size = " << m_synthesisWindowSize << ", fft size = " << m_fftSize << ", input hopsize = " << inputHopSize << " (approx) output hopsize = " << outputHopSize << endl;
    }

    size_t m_maxProcessSize = std::max(m_analyzeWindowSize, m_synthesisWindowSize);

    // This headroom is so as to try to avoid reallocation when
    // the pitch scale changes
    m_outbufSize = hsratio > 1 ? m_maxProcessSize * 16 * hsratio : m_maxProcessSize * 16;

    if (m_logLevel > 0) {
        cerr << "phasevocoder: outbuf size = " << m_outbufSize << endl;
    }
}

void
phasevocodercore::Impl::configure() {
    calculateSizes();

    m_analyzeWindow = new windowfunc<float>(Hanning, m_analyzeWindowSize);
    m_synthesisWindow = new windowfunc<float>(Hanning, m_synthesisWindowSize);
        
    // ring buffers (carrier and modulator) initialized
    for (size_t chan = 0; chan < m_audioData.size(); ++chan) {
        delete m_audioData[chan];
    }
    m_audioData.clear();

    for (size_t chan = 0; chan < m_carrierData.size(); ++chan) {
        delete m_carrierData[chan];
    }
    m_carrierData.clear();

    size_t max_WindowSize = std::max(m_analyzeWindowSize, m_synthesisWindowSize);

    for (size_t chan = 0; chan < m_channels; ++chan) {
        m_audioData.push_back
            (new channelinfo(max_WindowSize, m_fftSize, m_outbufSize));
        m_carrierData.push_back
            (new channelinfo(max_WindowSize, m_fftSize, m_outbufSize));
    }

    // resample buffer initialized
    for (size_t chan = 0; chan < m_channels; ++chan) {
        // rbs is the amount of buffer space we think we'll need
        // for resampling; but mem_allocate a sensible amount in case
        // the pitch scale changes during use
        size_t rbs = lrintf(ceil((m_hopsize * m_timeRatio * 2) / m_pitchScale));
        if (rbs < m_hopsize * 16) rbs = m_hopsize * 16;
        SetResampleBufSize(rbs);
    }

    for (size_t chan = 0; chan < rsb_gen.size(); ++chan) {
        delete rsb_gen[chan];
    }
    rsb_gen.clear();
    for (size_t chan = 0; chan < rsb_gen.size(); ++chan) {
        delete rsb_chordgen[chan];
    }
    rsb_chordgen.clear();

    
    for (size_t chan = 0; chan < m_channels; ++chan) {
        // rosenberg generator
        rsb_gen.push_back(new rosenberg(m_sampleRate, 440, 0.01, 0.06));
        // rosenberg chord generator
        // std::vector<float> chordmaj = {440, 554.365, 659.255}; // A major chord on A4
        // rsb_chordgen = new rosenbergchord(m_sampleRate, 0.01, 0.06, chordmaj);
        std::vector<float> chordmin = {440, 523.251, 659.255}; // A minor chord on A4
        rsb_chordgen.push_back(new rosenbergchord(m_sampleRate, 0.01, 0.06, chordmin));
    }
    
}

void
phasevocodercore::Impl::SetResampleBufSize(size_t sz)
{
    resamplebuf = mem_allocate_and_zero<float>(sz);
    resamplebufSize = sz;
    if (m_logLevel > 0) {
        cerr << "phasevocoder: resamplebufSize = " << resamplebufSize << endl;
    }
}


void
phasevocodercore::Impl::setDebugLevel(int level) {
    m_logLevel = level;
}	

void
phasevocodercore::Impl::processNormal(const float *const *input, size_t num_samples) {
    bool allread = false;

	std::vector<size_t> num_readsamples;
	num_readsamples.resize(m_channels);
    //size_t num_readsamples[m_channels];
    for (size_t chan = 0; chan < m_channels; ++chan) {
        num_readsamples[chan] = 0;
    }

    while (!allread) {
        for (size_t chan = 0; chan < m_channels; ++chan) {
            num_readsamples[chan] += enbufferChannel(chan,
                                          input,
                                          num_readsamples[chan],
                                          num_samples - num_readsamples[chan]);
            if (num_readsamples[chan] < num_samples) {
                allread = false;
            } else {
                allread = true;
            }

        }

        processOneSlice();

    }

}


void phasevocodercore::Impl::processConstant(const float *const *input, size_t num_samples) {
    bool allread = false;

	std::vector<size_t> num_readsamples;
	num_readsamples.resize(m_channels);
    for (size_t chan = 0; chan < m_channels; ++chan) {
        num_readsamples[chan] = 0;
    }
    while (!allread) {
        for (size_t chan = 0; chan < m_channels; ++chan) {
            num_readsamples[chan] += enbufferChannel(chan,
                                          input,
                                          num_readsamples[chan],
                                          num_samples - num_readsamples[chan]);

            if (num_readsamples[chan] < num_samples) {
                allread = false;
            } else {
                allread = true;
            }
        }

        processOneSliceConstant();
    }
}

void phasevocodercore::Impl::processVocoder(const float *const *input, size_t num_samples, carrierType ct) {
    bool allread = false;

	std::vector<size_t> num_readsamples;
	num_readsamples.resize(m_channels);
    for (size_t chan = 0; chan < m_channels; ++chan) {
        num_readsamples[chan] = 0;
    }
    while (!allread) {
        for (size_t chan = 0; chan < m_channels; ++chan) {
            num_readsamples[chan] += enbufferChannelVocoder(chan,
                                          input,
                                          num_readsamples[chan],
                                          num_samples - num_readsamples[chan], ct);

            if (num_readsamples[chan] < num_samples) {
                allread = false;
            } else {
                allread = true;
            }
        }

        // processOneSliceVocoderRealImag();
        processOneSliceVocoder();
    }
}


} // end of PhaseVocoder namespace

