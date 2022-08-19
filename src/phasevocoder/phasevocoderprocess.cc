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
#include <deque>
#include <algorithm>

#include "channelinfo.h"


#include "../common/dsp/resampler.h"
#include "../common/gen/rosenberg.h"
#include "../common/gen/rosenbergchord.h"
#include "../common/system/vectorization.h"

using std::cerr;
using std::endl;

namespace audiomod {

    
size_t
phasevocodercore::Impl::enbufferChannel(size_t chan,
                                          const float *const *inputs,
                                          size_t offset,
                                          size_t num_samples) {
    channelinfo &ad = *m_audioData[chan];
    circularqueue<float> &inbuf = *ad.inbuf;

    size_t toWrite = num_samples;
    size_t writable = inbuf.GetWriteSpace();

    if (writable < toWrite) {
        toWrite = writable;
    }

    const float *input = 0;
    input = inputs[chan] + offset;
    inbuf.write(input, toWrite);
    ad.incnt += toWrite;
    return toWrite;

}

size_t
phasevocodercore::Impl::enbufferChannelVocoder(size_t chan,
                                          const float *const *inputs,
                                          size_t offset,
                                          size_t num_samples, carrierType ct) {
    
    channelinfo &ad = *m_audioData[chan]; // input (modulator) data
    circularqueue<float> &inbuf = *ad.inbuf;
    

    channelinfo &caad = *m_carrierData[chan]; // carrier data
    circularqueue<float> &cainbuf = *caad.inbuf;

    size_t toWrite = num_samples;
    size_t writable = inbuf.GetWriteSpace();
    size_t cawritable = cainbuf.GetWriteSpace();

    // take the minimum of three as "toWrite"
    if (writable < toWrite) {
        toWrite = writable;
    }
    if (cawritable < toWrite) {
        toWrite = cawritable;
    }

    const float *input = 0; // offset pointer
    input = inputs[chan] + offset;
    inbuf.write(input, toWrite);

    // write rosenberg samples to cainbuf
    float *cainput = new float[toWrite];

    if (ct == rosenbergType) {
        for (size_t i=0; i<toWrite; i++) {
            cainput[i] = rsb_gen[chan]->getNextSample() * 0.3; // half the amplitude
        }
    }
    else if (ct == rosenbergChordType) {
        for (size_t i=0; i<toWrite; i++) {
            cainput[i] = rsb_chordgen[chan]->getNextSample() * 0.3; // half the amplitude
        }
    }
    else {
        std::cerr << "carrierType not supported!" << std::endl;
    }
    
    cainbuf.write(cainput, toWrite);
    delete[] cainput;
    cainput = nullptr;

    ad.incnt += toWrite;
    caad.incnt += toWrite;
    
    return toWrite;
}

int
phasevocodercore::Impl::processOneSliceConstant() {
    for (size_t chan = 0; chan < m_channels; ++chan) {
        if (!inbufReady(chan)) { // test whether there's enough input data in the inbuf
            // printf("early return!\n");
            return -1;
        }
        channelinfo &ad = *m_audioData[chan];
        size_t ready = ad.inbuf->GetReadSpace();
        // assert(ready >= m_analyzeWindowSize || ad.inputSize >= 0);
        ad.inbuf->touchread(ad.interfacebuffer, std::min(ready, m_analyzeWindowSize));
        ad.inbuf->discard(m_hopsize);
        analyzeSlice(chan);
    }

    int outframes = 0;
    for (size_t chan = 0; chan < m_channels; ++chan) {
        synthesiseSlice(chan, m_hopsize);
        channelinfo &ad = *m_audioData[chan];
        size_t ws = ad.outbuf->GetWriteSpace();
        if (ws < m_hopsize) {
            // if the outbuf is initialized properly, this should not be entered
            // circularqueue<float> *oldbuf = ad.outbuf;
            // ad.outbuf = oldbuf->resized(oldbuf->GetSize() + (m_hopsize - ws));
            // m_emergencyScavenger.claim(oldbuf);
            cerr << "Buffer overrun on output for channel, return 0 " << chan << endl;
            return 0;
        }
        outframes = writeSlice(chan, m_hopsize, false);
        m_audioData[chan]->slicecnt++;
    }

    return outframes;

}

int
phasevocodercore::Impl::processOneSliceVocoder() {
    for (size_t chan = 0; chan < m_channels; ++chan) {
        if (!inbufReady(chan)) { // test whether there's enough input data in the inbuf
            // printf("early return!\n");
            return -1;
        }
        channelinfo &ad = *m_audioData[chan];
        size_t ready = ad.inbuf->GetReadSpace();
        // assert(ready >= m_analyzeWindowSize || ad.inputSize >= 0);
        ad.inbuf->touchread(ad.interfacebuffer, std::min(ready, m_analyzeWindowSize));
        ad.inbuf->discard(m_hopsize);
        analyzeSlice(chan);
    }

    // assume that rsb is always in sync with input data
    for (size_t chan = 0; chan < m_channels; ++chan) {
        // Shall we test cdinbuf if there's enough data?
        channelinfo &caad = *m_carrierData[chan];
        size_t ready = caad.inbuf->GetReadSpace();
        // assert(ready >= m_analyzeWindowSize || caad.inputSize >= 0);
        caad.inbuf->touchread(caad.interfacebuffer, std::min(ready, m_analyzeWindowSize));
        caad.inbuf->discard(m_hopsize);
        analyzeCarrierSlice(chan);
    }

    for (size_t chan = 0; chan < m_channels; ++chan) {
        modifySliceVocoder(chan, m_hopsize);

        synthesiseSliceCarrier(chan, m_hopsize);

        writeSliceCarrier(chan, m_hopsize);
        m_audioData[chan]->slicecnt++;
        m_carrierData[chan]->slicecnt++;
    }

    return 0;
}

// int
// phasevocodercore::Impl::processOneSliceVocoderRealImag() {
//     for (size_t chan = 0; chan < m_channels; ++chan) {
//         if (!inbufReady(chan)) { // test whether there's enough input data in the inbuf
//             // printf("early return!\n");
//             return -1;
//         }
//         channelinfo &ad = *m_audioData[chan];
//         size_t ready = ad.inbuf->GetReadSpace();
//         // assert(ready >= m_analyzeWindowSize || ad.inputSize >= 0);
//         ad.inbuf->touchread(ad.interfacebuffer, std::min(ready, m_analyzeWindowSize));
//         ad.inbuf->discard(m_hopsize);
//         analyzeSliceRealImag(chan);
//     }

//     // assume that rsb is always in sync with input data
//     for (size_t chan = 0; chan < m_channels; ++chan) {
//         // Shall we test cdinbuf if there's enough data?
//         channelinfo &caad = *m_carrierData[chan];
//         size_t ready = caad.inbuf->GetReadSpace();
//         // assert(ready >= m_analyzeWindowSize || caad.inputSize >= 0);
//         caad.inbuf->touchread(caad.interfacebuffer, std::min(ready, m_analyzeWindowSize));
//         caad.inbuf->discard(m_hopsize);

//         analyzeCarrierSliceRealImag(chan);
//     }

//     for (size_t chan = 0; chan < m_channels; ++chan) {
//         modifySliceVocoderRealImag(chan, m_hopsize);

//         synthesiseSliceCarrierRealImag(chan, m_hopsize);

//         writeSliceCarrier(chan, m_hopsize);
//         m_audioData[chan]->slicecnt++;
//         m_carrierData[chan]->slicecnt++;
//     }
//     return 0;
// }

int
phasevocodercore::Impl::processOneSlice() {
    //Profiler profiler("phasevocodercore::Impl::processOneSlice");

    // Process a single slice for all channels, provided there is
    // enough data on each channel for at least one slice.  This is
    // able to calculate increments as it goes along.

    // This is the normal process method in RT mode.

    for (size_t chan = 0; chan < m_channels; ++chan) {
        if (!inbufReady(chan)) { // test whether there's enough input data in the inbuf
            if (m_logLevel > 2) {
                cerr << "processOneSlice: out of input" << endl;
            }
            printf("out of input data, early return!\n");
            return -1;
        }
        channelinfo &ad = *m_audioData[chan];

        size_t ready = ad.inbuf->GetReadSpace();
        // printf("inbuf ready:%d, ", ready);
        // assert(ready >= m_analyzeWindowSize || ad.inputSize >= 0);
        ad.inbuf->touchread(ad.interfacebuffer, std::min(ready, m_analyzeWindowSize));
        ad.inbuf->discard(m_hopsize);
        analyzeSlice(chan);

    }
    
    size_t phaseIncrement, shiftIncrement;

    if (m_options & OptionRobotic || m_options & OptionWhisper) {
        phaseIncrement = m_hopsize;
        shiftIncrement = m_hopsize;
    } else {
        if (isIntRatio()) {
            phaseIncrement = m_hopsize * getHopSizeRatio();
            shiftIncrement = m_hopsize * getHopSizeRatio();
        } else {
            calculateIncrements(phaseIncrement, shiftIncrement);
        }
    }
    

    int outframes = 0;
    for (size_t chan = 0; chan < m_channels; ++chan) {
        outframes = processSliceForChannel(chan, phaseIncrement, shiftIncrement);
        m_audioData[chan]->slicecnt++;
    }

    return outframes;
}

bool
phasevocodercore::Impl::inbufReady(size_t chan) {

    channelinfo &ad = *m_audioData[chan];
    circularqueue<float> &inbuf = *ad.inbuf;

    size_t rs = inbuf.GetReadSpace();

    if (rs < m_analyzeWindowSize) {
        return false;
    } else {
        return true;
    }

}

int 
phasevocodercore::Impl::processSliceForChannel(size_t chan,
                                                  size_t phaseIncrement,
                                                  size_t shiftIncrement) {

    // Process a single slice on a single channel.  This assumes
    // enough input data is available; caller must have tested this
    // using e.g. inbufReady first.  Return true if this is
    // the last slice on the channel.

    channelinfo &ad = *m_audioData[chan];

    if (m_options & OptionRobotic) { // robotic mode
        roboticSlice(chan);
    } else if (m_options & OptionWhisper) { // whispering mode
        whisperSlice(chan);
    } else { // time-stretching and pitch-shifting mode
        if (m_defaultCoreMode == 0) {
            modifySliceSimple(chan, phaseIncrement);
        } else if (m_defaultCoreMode == 1) {
            modifySlicePhaseLocked(chan, phaseIncrement);
        } else if (m_defaultCoreMode == 2) {
            modifySliceIntRatio(chan, phaseIncrement);
        } else {
            modifySliceSimple(chan, phaseIncrement);
        }
    }

    synthesiseSlice(chan, shiftIncrement); // reads from ad.mag, ad.phase

    // tangkk: the required number of output samples (after resampling) at this pass
    // consider only the aftereffect of m_pitchScale because m_timeRatio has already been considered in shiftIncrement
    int required = int(shiftIncrement / m_pitchScale) + 1;

    int ws = ad.outbuf->GetWriteSpace();

    // printf("ws:%d, required:%d,", ws, required);

    /* ws should be >= required, otherwise enter the following routine */
    if (ws < required) {
        // if (m_logLevel > 0) {
        //     cerr << "Buffer overrun on output for channel " << chan << endl;
        // }
        // // if the outbuf is initialized properly, this should not be entered

        // // The only correct thing we can do here is resize the buffer.
        // // We can't wait for the client thread to read some data out
        // // from the buffer so as to make more space, because the
        // // client thread (if we are threaded at all) is probably stuck
        // // in a process() call waiting for us to stow away enough
        // // input increments to allow the process() call to complete.
        // // This is an unhappy situation.

        // circularqueue<float> *oldbuf = ad.outbuf;
        // ad.outbuf = oldbuf->resized(oldbuf->GetSize() + (required - ws));
        // m_emergencyScavenger.claim(oldbuf);

        cerr << "Buffer overrun on output for channel, return 0 " << chan << endl;
        return 0;
    }

    /* 
    write the processed result to ad.outbuf
    note that only shiftIncrement amount of samples are written out
    this amount might be less than or more than a window size
    thus the writing speed is non-consistent
    */
    // printf("shiftIncrement:%d, ", shiftIncrement);
    bool last = false; // FIXME: apply this information later
    int outframes = writeSlice(chan, shiftIncrement, last); 
    return outframes;
}

// FIXME: please deal with ratio change here
static int calculateThisIncrement(float ratio, size_t increment, size_t samplerate) {
    static float m_recovery = 0;
    static float m_prevRatio = ratio; // FIXME: should deal with previous ratio here
    // bool ratioChanged = (ratio != m_prevRatio);
    // m_prevRatio = ratio;
    static float m_divergence = 0;

    // if (ratioChanged) {
    // m_recovery (as well as m_divergence) is very important to compensate for non-integer ratio,
    // to preserve the integrety of samples
    m_recovery = m_divergence / (( samplerate / 10.0) / increment);
    // }

    int incr = lrint(increment * ratio - m_recovery);

    if (incr < lrint((increment * ratio) / 2)) {
        incr = lrint((increment * ratio) / 2);
    } else if (incr > lrint(increment * ratio * 2)) {
        incr = lrint(increment * ratio * 2);
    }

    float divdiff = (increment * ratio) - incr;
    float prevDivergence = m_divergence;
    m_divergence -= divdiff;
    if ((prevDivergence < 0 && m_divergence > 0) ||
        (prevDivergence > 0 && m_divergence < 0)) {
        m_recovery = m_divergence / ((samplerate / 10.0) / increment);
    }

    return incr;

}

void
phasevocodercore::Impl::calculateIncrements(size_t &phaseIncrementRtn, size_t &shiftIncrementRtn) {

    // Calculate the next upcoming phase and shift increment, on the
    // basis that both channels are in sync.  This is in contrast to
    // getIncrements, which requires that all the increments have been
    // calculated in advance but can then return increments
    // corresponding to different slices in different channels.

    // Requires frequency domain representations of channel data in
    // the mag and phase buffers in the channel.

    // This function is only used in real-time mode.

    phaseIncrementRtn = m_hopsize;
    shiftIncrementRtn = m_hopsize;

    if (m_channels == 0) return;

    channelinfo &ad = *m_audioData[0];

    size_t bc = ad.slicecnt;
    for (size_t chan = 1; chan < m_channels; ++chan) {
        if (m_audioData[chan]->slicecnt != bc) {
            cerr << "ERROR: phasevocodercore::Impl::calculateIncrements: Channels are not in sync" << endl;
            return;
        }
    }

    const int halfsize = m_fftSize/2 + 1;

    // Normally we would mix down the time-domain signal and apply a
    // single FFT, or else mix down the Cartesian form of the
    // frequency-domain signal.  Both of those would be inefficient
    // from this position.  Fortunately, the onset detectors should
    // work reasonably well (maybe even better?) if we just sum the
    // magnitudes of the frequency-domain channel signals and forget
    // about phase entirely.  Normally we don't expect the channel
    // phases to cancel each other, and broadband effects will still
    // be apparent.

    // float df = 0.f;
    // bool silent = false;

    // FIXME: is df really necessary?

    // int incr = m_stretchCalculator->calculateSingle
    //     (getHopSizeRatio(), df, m_hopsize);

    int incr = calculateThisIncrement(getHopSizeRatio(), m_hopsize, m_sampleRate);

    
    // The returned increment (incr) is the phase increment.  The shift
    // increment for one slice is the same as the phase increment for
    // the following slice (see comment below).  This means we don't
    // actually know the shift increment until we see the following
    // phase increment... which is a bit of a problem.

    // This implies we should use this increment for the shift
    // increment, and make the following phase increment the same as
    // it.  This means in RT mode we'll be one slice later with our
    // phase reset than we would be in non-RT mode.  The sensitivity
    // of the broadband onset detector may mean that this isn't a
    // problem -- test it and see.

    shiftIncrementRtn = incr;
    // shiftIncrementRtn = m_hopsize * getHopSizeRatio(); // naive approach: try that
    // printf("incr:%d, m_hopsize * getHopSizeRatio(): %f\n", incr, m_hopsize * getHopSizeRatio());

    if (ad.prev_increment == 0) {
        phaseIncrementRtn = shiftIncrementRtn;
    } else {
        phaseIncrementRtn = ad.prev_increment;
    }

    ad.prev_increment = shiftIncrementRtn;

}


void
phasevocodercore::Impl::analyzeSlice(size_t channel)
{
    channelinfo &ad = *m_audioData[channel];

    data_type *const R__ internalbuffer = ad.internalbuffer;
    float *const R__ interfacebuffer = ad.interfacebuffer;

    fftshift(internalbuffer, m_fftSize, interfacebuffer, m_analyzeWindow);

    ad.fft->forwardPolar(internalbuffer, ad.mag, ad.phase);
}

// void
// phasevocodercore::Impl::analyzeSliceRealImag(size_t channel)
// {

//     channelinfo &ad = *m_audioData[channel];

//     data_type *const R__ internalbuffer = ad.internalbuffer;
//     float *const R__ interfacebuffer = ad.interfacebuffer;

//     // ad.interfacebuffer is known to contain m_analyzeWindowSize samples

//     // we assume m_analyzeWindowSize always <= m_fftSize
//     // if (m_analyzeWindowSize > m_fftSize) {
//     //     m_afilter->apply(interfacebuffer);
//     // }

//     fftshift(internalbuffer, m_fftSize, interfacebuffer, m_analyzeWindow);

//     ad.fft->forward(internalbuffer, ad.real, ad.imag);
// }

void
phasevocodercore::Impl::analyzeCarrierSlice(size_t channel) {
    channelinfo &caad = *m_carrierData[channel];
    data_type *const R__ internalbuffer = caad.internalbuffer;
    float *const R__ interfacebuffer = caad.interfacebuffer;

    // we assume m_analyzeWindowSize always <= m_fftSize
    // if (m_analyzeWindowSize > m_fftSize) {
    //     m_afilter->apply(interfacebuffer);
    // }

    fftshift(internalbuffer, m_fftSize, interfacebuffer, m_analyzeWindow);

    caad.fft->forwardPolar(internalbuffer, caad.mag, caad.phase);
}

// void
// phasevocodercore::Impl::analyzeCarrierSliceRealImag(size_t channel) {
//     channelinfo &caad = *m_carrierData[channel];
//     data_type *const R__ internalbuffer = caad.internalbuffer;
//     float *const R__ interfacebuffer = caad.interfacebuffer;

//     // we assume m_analyzeWindowSize always <= m_fftSize
//     // if (m_analyzeWindowSize > m_fftSize) {
//     //     m_afilter->apply(interfacebuffer);
//     // }

//     fftshift(internalbuffer, m_fftSize, interfacebuffer, m_analyzeWindow);

//     caad.fft->forward(internalbuffer, caad.real, caad.imag);
// }

void
phasevocodercore::Impl::modifySliceIntRatio(size_t channel, size_t phaseIncrement)
{
    //Profiler profiler("phasevocodercore::Impl::modifySliceIntRatio");

    channelinfo &ad = *m_audioData[channel];

    const int halfsize = m_fftSize / 2; // half window size (hWLen)
    
    for (int i=0; i<halfsize; i++) {
        ad.phase[i] = ad.phase[i] * phaseIncrement / m_hopsize; // note this should be int ratio
        // ad.phase[i] = ad.phase[i] * getHopSizeRatio();
    }

} 

void
phasevocodercore::Impl::modifySlicePhaseLocked(size_t channel, size_t phaseIncrement) {
    //Profiler profiler("phasevocodercore::Impl::modifySlicePhaseLocked");

    channelinfo &ad = *m_audioData[channel];

    const int halfsize = m_fftSize / 2; // half window size (hWLen)

    // find spectral peaks - simple method
    peak_loc.clear();
    data_type *const R__ mag = ad.mag;
    

    int b=2;
    while (b + 2 < halfsize) { // note to prevent from array out of bounce
        if (mag[b] > mag[b-1] && mag[b] > mag[b-2] && mag[b] > mag[b+1] && mag[b] > mag[b+2]) {
            // local maxima found
            peak_loc.push_back(b);
            b+=3;
        } else {
            b+=1;
        }
    }

    // std::vector<float> magvec(mag, mag + halfsize);
    // findPeaks(magvec, peak_loc, 1, 0.5); // NOTE That this is float based


    static bool firstentry = true;

    // printf("process block......\n");

    if (firstentry) {
        // printf("process init......\n");
        // init case
        for (int i=0; i<halfsize; i++) {
            data_type tp = ad.phase[i];
            // data_type perr = 0.0;
            data_type outphase = tp;
            ad.prev_phase[i] = tp;
            ad.phase[i] = outphase;
            ad.prev_outphase[i] = outphase;
        }
    } else if (peak_loc.size() == 0 || prev_peak_loc.size() == 0) {
        // printf("process normal......\n");
        // back to normal case
        for (int i=0; i<halfsize; i++) {
            // data_type tp = ad.phase[i];
            // data_type perr = 0.0;
            // data_type outphase = tp;
            // // normal process here
            data_type omega = (2 * M_PI * m_hopsize * i) / (m_fftSize);
            // data_type pp = ad.prev_phase[i]; // pp = previous phase
            // data_type ep = pp + omega; // ep = evolved phase (prev_phase + phase shift)
            // perr = princarg(tp - ep); // diff between this phase and evolved phase. wrap it in the range of [-PI,PI]
            // data_type advance = phaseIncrement * ((omega + perr) / m_hopsize);
            data_type delta_phi = omega +  princarg(ad.phase[i] - ad.prev_phase[i] - omega);
            data_type advance = delta_phi * phaseIncrement / m_hopsize;
            data_type outphase = princarg(ad.prev_outphase[i] + advance);
            ad.prev_phase[i] = ad.phase[i];
            ad.phase[i] = outphase;
            ad.prev_outphase[i] = outphase;
        }
    } else {
        // printf("process phase-locked......\n");
        // phase-locked process here
        size_t prev_p = 0; // global to the following loop
        for (size_t p = 0; p < peak_loc.size(); p++) {
            int p2 = peak_loc[p]; // new peak
            // connect current peak to the previous closest peak
            while (prev_p < prev_peak_loc.size() - 1) {
                if (abs(p2-prev_peak_loc[prev_p+1]) < abs(p2-prev_peak_loc[prev_p])) {
                    prev_p += 1;
                } else {
                    break;
                }
            }
            if (prev_p == prev_peak_loc.size()) prev_p--;
            int p1 = prev_peak_loc[prev_p]; // now p1 and p2 are connected (they have closest peak loc)

            // printf("p2:%d, p1:%d\n", p2, p1);

            // propagate peak's phase assuming linear frequency
            float avg_p = (p1 + p2) * 0.5;
            data_type pomega = (2 * M_PI * m_hopsize * (avg_p - 1)) / (m_fftSize);
            // data_type pp = ad.prev_phase[p1];
            // data_type ep = pp + pomega;
            // data_type perr = princarg(ad.phase[p2] - ep);
            // data_type peak_delta_phi = pomega + perr;
            data_type peak_delta_phi = pomega + princarg(ad.phase[p2] - ad.prev_phase[p1] - pomega);
            data_type peak_target_phase = princarg(ad.prev_outphase[p1] + (peak_delta_phi * phaseIncrement) / m_hopsize);
            data_type peak_phase_rotation = princarg(peak_target_phase - ad.phase[p2]);

            // rotate phases of all bins around the current peak
            int bin1 = 0;
            int bin2 = 0;

            if (peak_loc.size() == 1) {
                bin1 = 0;
                bin2 = halfsize;
            } else if (p == 0) {
                bin1 = 0;
                bin2 = round((peak_loc[p+1] + p2) * 0.5);
            } else if (p == peak_loc.size()-1) {
                bin1 = round((peak_loc[p-1] + p2) * 0.5);
                bin2 = halfsize;
            } else { // boundaries already guarded by the above conditions
                bin1 = round((peak_loc[p-1] + p2) * 0.5); // NOTE here p-1
                bin2 = round((peak_loc[p+1] + p2) * 0.5); // NOTE here p+1
            }

            // printf("bin1:%d, bin2:%d\n", bin1, bin2);

            for (int i=bin1; i<bin2; i++) {
                ad.locked_phase[i] = princarg(ad.phase[i] + peak_phase_rotation);
            }
        }
        // assign locked_phase to Phase
        // vector_copy(ad.phase, ad.locked_phase, count);

        // update internal states
        for (int i=0; i<halfsize; i++) {
            ad.prev_phase[i] = ad.phase[i];
            ad.prev_outphase[i] = ad.locked_phase[i];
            ad.phase[i] = ad.locked_phase[i];
        }
    }

    // update global status
    prev_peak_loc = peak_loc;
    firstentry = false; // note that this is a static variable
    
}

void
phasevocodercore::Impl::modifySliceSimple(size_t channel, size_t phaseIncrement) {
    //Profiler profiler("phasevocodercore::Impl::modifySliceSimple");

    channelinfo &ad = *m_audioData[channel];

    const int halfsize = m_fftSize / 2; // half window size (hWLen)

    static bool firstentry = true;

    if (firstentry) {
        // printf("process init......\n");
        // init case
        for (int i=0; i<halfsize; i++) {
            data_type tp = ad.phase[i];
            // data_type perr = 0.0;
            data_type outphase = tp;
            ad.prev_phase[i] = tp;
            ad.phase[i] = outphase;
            ad.prev_outphase[i] = outphase;
        }
    } else if (peak_loc.size() == 0 || prev_peak_loc.size() == 0) {
        // printf("process normal......\n");
        // back to normal case
        for (int i=0; i<halfsize; i++) {
            // data_type tp = ad.phase[i];
            // data_type perr = 0.0;
            // data_type outphase = tp;
            // // normal process here
            data_type omega = (2 * M_PI * m_hopsize * i) / (m_fftSize);
            // data_type pp = ad.prev_phase[i]; // pp = previous phase
            // data_type ep = pp + omega; // ep = evolved phase (prev_phase + phase shift)
            // perr = princarg(tp - ep); // diff between this phase and evolved phase. wrap it in the range of [-PI,PI]
            // data_type advance = phaseIncrement * ((omega + perr) / m_hopsize);
            data_type delta_phi = omega +  princarg(ad.phase[i] - ad.prev_phase[i] - omega);
            data_type advance = delta_phi * phaseIncrement / m_hopsize;
            data_type outphase = princarg(ad.prev_outphase[i] + advance);
            ad.prev_phase[i] = ad.phase[i];
            ad.phase[i] = outphase;
            ad.prev_outphase[i] = outphase;
        }
    }

    firstentry = false; // note that this is a static variable

}

void
phasevocodercore::Impl::modifySliceVocoder(size_t channel, size_t phaseIncrement) {
    int num_bands = 512;
    int band_len = (int)floor(float(m_fftSize) / float(num_bands * 2));
    channelinfo &ad = *m_audioData[channel];
    channelinfo &caad = *m_carrierData[channel];

    for (int band_no = 0; band_no < num_bands; band_no++) {
        float mean_modul_mag = 0;
        for (int i = 0, j = band_no * band_len; i < band_len; i++, j++)
        {
            mean_modul_mag += ad.mag[j];
        }
        mean_modul_mag /= (band_len * 2);
        for (int i = 0, j = band_no * band_len; i < band_len; i++, j++) {
            caad.mag[j] *= mean_modul_mag;
        }
        // zero the DC offset introduced by the carrier (e.g. rosenberg pulse)
        caad.mag[0] = 0;
        caad.mag[caad.GetRealSize()-1] = 0;
    }
}

// void
// phasevocodercore::Impl::modifySliceVocoderRealImag(size_t channel, size_t phaseIncrement) {
//     int num_bands = 64;
//     int band_len = (int)floor(float(m_fftSize) / float(num_bands * 2));
//     channelinfo &ad = *m_audioData[channel];
//     channelinfo &caad = *m_carrierData[channel];

//     for (int band_no = 0; band_no < num_bands; band_no++) {
//         float mean_modul_mag = 0;
//         for (int i = 0, j = band_no * band_len; i < band_len; i++, j++)
//         {
//             mean_modul_mag += sqrt(pow(ad.real[j],2) + pow(ad.imag[j],2));
//         }
//         mean_modul_mag /= (band_len * 2);
//         for (int i = 0, j = band_no * band_len, k = m_fftSize - j; i < band_len; i++, j++, k--) {
//             caad.real[j] *= mean_modul_mag;
//             caad.imag[j] *= mean_modul_mag;
//             if (j==0) continue;
//             caad.real[k] = caad.real[j];
//             caad.imag[k] = -caad.imag[j]; // negative because this is equivalent to fftshift
//         }
//     }
//     ad.real[m_fftSize/2] *= caad.real[m_fftSize/2];
//     ad.imag[m_fftSize/2] *= caad.real[m_fftSize/2];
    
// }

void
phasevocodercore::Impl::roboticSlice(size_t channel) {
    channelinfo &ad = *m_audioData[channel];
    int realsize = ad.GetRealSize();
    for (int i=0; i<realsize; i++) {
        ad.phase[i] = 0;
    }
}

void
phasevocodercore::Impl::whisperSlice(size_t channel) {
    channelinfo &ad = *m_audioData[channel];
    int realsize = ad.GetRealSize();
    float two_pi = 2*M_PI;
    for (int i=0; i<realsize; i++) {
        ad.phase[i] = two_pi * (float)rand() / (float)RAND_MAX;
    }
}

void
phasevocodercore::Impl::maleToFemale(size_t channel) {
    // formantShiftSlice(channel, 0.85); // FIXME: magic number here
    freqCompSlice(channel, 0.85 * m_pitchScale); // FIXME: magic number here
}

void
phasevocodercore::Impl::femaleToMale(size_t channel) {
    // formantShiftSlice(channel, 1.17); // FIXME: magic number here
    freqCompSlice(channel, 1.17 * m_pitchScale); // FIXME: magic number here
}

void
phasevocodercore::Impl::formantPreserveSlice(size_t channel) {
    // formantShiftSlice(channel, 1);
    freqCompSlice(channel, m_pitchScale);
}

void
phasevocodercore::Impl::freqCompSlice(size_t channel, float freq_comp) {
    //Profiler profiler("phasevocodercore::Impl::freqCompSlice");

    channelinfo &ad = *m_audioData[channel];

    // data_type *const R__ tmpphi = ad.locked_phase;
    data_type *const R__ phi = ad.phase;
    data_type *const R__ mag = ad.mag;
    // data_type *const R__ envelope = ad.envelope;
    // data_type *const R__ internalbuffer = ad.internalbuffer;

    const int halfsize = m_fftSize / 2;
    // const data_type factor = 1.0 / m_fftSize;

    // vector_copy(envelope, mag, ad.GetRealSize()); // whole mag as envelope

    // vector_divide(mag, envelope, halfsize + 1); // whiten the magnitude spectrum by the spectral envelope

    // float envScale = m_pitchScale * env_pitch_scale;
    // printf("envScale:%f\n",envScale);
    // const float fixedgain = freq_comp; // compensate for the zeroed bins
    float absps = m_pitchScale > 1 ? m_pitchScale : 1 / m_pitchScale;
    // const float fixedgain = freq_comp * absps; // compensate for the zeroed bins
    const float fixedgain = absps;
    // printf("absps:%f, fixedgain:%f\n", absps, fixedgain);

    if (freq_comp > 1.0) {
        // spectral compression
        for (int target = 0; target <= halfsize; ++target) { // target counting up
            int source = lrint(target * freq_comp); // source > target
            data_type delta_omega = (2 * M_PI * m_hopsize * (target - source)) / (m_fftSize);
            // data_type delta_omega = (2 * M_PI * m_hopsize * getHopSizeRatio() * (target - source)) / (m_fftSize);
            
            if (source > halfsize) {
                // mag[target] = mag[source - halfsize];
                mag[target] = 0.0;
                // envelope[target] = 0.0;
                // tmpphi[target] = 0.0; // add phi consideration
                // phi[target] = phi[source - halfsize] + delta_omega;
                // phi[target] = rand() % (6);
                phi[target] = 0.0;
            } else {
                mag[target] = mag[source];
                // envelope[target] = envelope[source];
                // tmpphi[target] = phi[source];
                // tmpphi[target] = phi[source] + delta_omega;
                // phi[target] = princarg(phi[source] + delta_omega);
                phi[target] = phi[source] + delta_omega;
                // phi[target] = phi[source]; // add phi consideration
                // phi[target] = phi[source] * (envelope[target] / envelope[source]);
            }
        }
    } else {
        // spectral expansion
        for (int target = halfsize; target > 0; ) { // target counting down
            --target;
            int source = lrint(target * freq_comp); // source < target
            data_type delta_omega = (2 * M_PI * m_hopsize * (target - source)) / (m_fftSize);
            // data_type delta_omega = (2 * M_PI * m_hopsize * getHopSizeRatio() * (target - source)) / (m_fftSize);
            mag[target] = mag[source];
            // envelope[target] = envelope[source];
            // tmpphi[target] = phi[source];
            // tmpphi[target] = phi[source] + delta_omega;
            // phi[target] = princarg(phi[source] + delta_omega);
            phi[target] = phi[source] + delta_omega;
            // phi[target] = phi[source]; // add phi consideration
            // phi[target] = phi[source] * (envelope[target] / envelope[source]);
        }
    }

    // modify the whitened magnitude with the new envelope
    // vector_multiply(mag, envelope, halfsize + 1);
    // vector_multiply(phi, envelope, halfsize + 1);
    vector_gain(mag, fixedgain, halfsize + 1);

    // update phase information
    // for (int i=0; i<=halfsize; i++) {
    //     phi[i] = tmpphi[i];
    // }

}

void
phasevocodercore::Impl::formantShiftSlice(size_t channel, float env_comp)
{
    channelinfo &ad = *m_audioData[channel];

    // data_type *const R__ phi = ad.phase;
    data_type *const R__ mag = ad.mag;
    data_type *const R__ envelope = ad.envelope;
    data_type *const R__ internalbuffer = ad.internalbuffer;

    const int halfsize = m_fftSize / 2;
    const data_type factor = 1.0 / m_fftSize;

    // vector_copy(envelope, mag, ad.GetRealSize()); // this also works in gender change mode, why?

    // internalbuffer is the complex cepstrum, and it is in quefrency domain
    // basically it is equivalent to invfft(log(mag))
    ad.fft->inverseCepstral(mag, internalbuffer);

    // const int cutoff = m_sampleRate / 700; // try modify this one
    const int cutoff = 60;
    // const int cutoff = ad.GetRealSize(); // this preserve things much better (why this works?)


    internalbuffer[0] /= 2;
    internalbuffer[cutoff-1] /= 2;


    // only keep "cutoff" amount of cepstral coefficents
    for (size_t i = cutoff; i < m_fftSize; ++i) {
        internalbuffer[i] = 0.0;
    }

    vector_gain(internalbuffer, factor, cutoff);

    data_type *spare = (data_type *)alloca((halfsize + 1) * sizeof(data_type)); // to contain the phase information
    ad.fft->forward(internalbuffer, envelope, spare);

    vector_exp(envelope, halfsize + 1); // do exp because the previous computation of envelop is based on the log-spectrum
    vector_divide(mag, envelope, halfsize + 1); // whiten the magnitude spectrum by the spectral envelope

    // scale the envelope according to m_pitchScale
    // float envScale = env_pitch_scale;
    // float envScale = m_pitchScale * 0.85; // try 0.8 - 0.9, and use a pitchshift of an octave (12)
    // float envScale = m_pitchScale * env_pitch_scale;

    if (env_comp > 1.0) {
        // scaling up, we want a new envelope that is lower by the pitch factor
        for (int target = 0; target <= halfsize; ++target) {
            int source = lrint(target * env_comp);
            if (source > halfsize) {
                envelope[target] = 0.0;
                // phi[target] = 0.0;
            } else {
                envelope[target] = envelope[source];
                // phi[target] = phi[source];
            }
        }
    } else {
        // scaling down, we want a new envelope that is higher by the pitch factor
        for (int target = halfsize; target > 0; ) {
            --target;
            int source = lrint(target * env_comp);
            envelope[target] = envelope[source];
            // phi[target] = phi[source];
        }
    }

    // modify the whitened magnitude with the new envelope
    // FIXME: shouldn't the same enveloped be applied at the end of the flow in order to preserve formant?
    vector_multiply(mag, envelope, halfsize+1);

    // FIXME: phase will become inconsistent again after this process, here we should modify the phase accordingly

}

void
phasevocodercore::Impl::synthesiseSlice(size_t channel,
                                           size_t shiftIncrement)
{
    /* post process region */
    if ((m_options & OptionFormantPreserved) && (m_pitchScale != 1.0)) {
        formantPreserveSlice(channel);
    }

    if ((m_options & OptionGenderChange) && (m_pitchScale != 1.0)) {
        if (m_pitchScale > 1) {
            maleToFemale(channel); // better use this with an octave pitch shift
        } else {
            femaleToMale(channel);
        }
    } else if (m_options & OptionGenderChange) {
        // formantShiftSlice(channel, 0.5); // male to female default
        freqCompSlice(channel, 0.8);
        // freqCompSlice(channel, 0.5);
        // freqCompSlice(channel, 1.7);
        // formantShiftSlice(channel, 2);
    }

    /* synthesis region */
    channelinfo &ad = *m_audioData[channel];

    data_type *const R__ internalbuffer = ad.internalbuffer;
    float *const R__ interfacebuffer = ad.interfacebuffer;
    float *const R__ outputAccumulator = ad.outputAccumulator;
    float *const R__ windowAccumulator = ad.windowAccumulator;
    
    const int halfsize = m_fftSize / 2;

    // Our FFTs produced unscaled results. Scale before inverse
    // transform rather than after, to avoid overflow if using a
    // fixed-point FFT.
    float factor = 1.f / m_fftSize;
    vector_gain(ad.mag, factor, halfsize + 1);

    ad.fft->inversePolar(ad.mag, ad.phase, ad.internalbuffer); // tangkk:: ifft here

    ifftshift(interfacebuffer, m_fftSize, internalbuffer, m_synthesisWindow);
    // vector_assign(interfacebuffer, internalbuffer + halfsize, halfsize);
    // vector_assign(interfacebuffer + halfsize, internalbuffer, halfsize);
    // m_synthesisWindow->apply(interfacebuffer); // apply the synthesis window to the ifft outcome

    /* 
    Collect the outputAccumulator information - the overlapped synthesized output information
    now outputAccumulator contains the windowed ifft output
    plus the previous "shifted version" of output
    thus the OLA happens here
    the shifting happens in "WriteSlice":
    vector_move(outputAccumulator, outputAccumulator + shiftIncrement, m_synthesisWindowSize - shiftIncrement);
    vector_zeros(outputAccumulator + m_synthesisWindowSize - shiftIncrement, shiftIncrement);
    outputAccumulator and windowAccumulator are in-sync
    */
    vector_add(outputAccumulator, interfacebuffer, m_synthesisWindowSize);
    // ad.accumulatorFill = wsz;

    /*
    Collect windowAccumulator information - the overlapped window information
    to be used in the OLA process
    now windowAccumulator contains the window vector + m_analyzeWindow->GetArea()
    plus the previous "shifted version" of the previous iteration - thus achieve the overlap and add windows
    the shifting happens in "WriteSlice":
    vector_move(windowAccumulator, windowAccumulator + shiftIncrement, m_synthesisWindowSize - shiftIncrement);
    vector_zeros(windowAccumulator + m_synthesisWindowSize - shiftIncrement, shiftIncrement);
    windowAccumulator and outputAccumulator are in-sync
    */
   // add awindow's area (scalar) to swindow (vector) and copy the values to windowAccumulator
   // area*1.5 because we have to keep the original volume of the input audio
   // FIXME: can we get rid of the 1.5 magic number here?
   m_synthesisWindow->add2dst(windowAccumulator, m_analyzeWindow->GetArea() * 1.5);

}

void
phasevocodercore::Impl::synthesiseSliceCarrier(size_t channel,size_t shiftIncrement) {
    channelinfo &caad = *m_carrierData[channel];

    data_type *const R__ internalbuffer = caad.internalbuffer;
    float *const R__ interfacebuffer = caad.interfacebuffer;
    float *const R__ outputAccumulator = caad.outputAccumulator;
    float *const R__ windowAccumulator = caad.windowAccumulator;
    
    const int halfsize = m_fftSize / 2;

    float factor = 1.f / m_fftSize;
    vector_gain(caad.mag, factor, halfsize + 1);

    caad.fft->inversePolar(caad.mag, caad.phase, caad.internalbuffer); // tangkk:: ifft here

    ifftshift(interfacebuffer, m_fftSize, internalbuffer, m_synthesisWindow);
    // vector_assign(interfacebuffer, internalbuffer + halfsize, halfsize);
    // vector_assign(interfacebuffer + halfsize, internalbuffer, halfsize);
    // m_synthesisWindow->apply(interfacebuffer);

    vector_add(outputAccumulator, interfacebuffer, m_synthesisWindowSize);
    // caad.accumulatorFill = wsz;

    // we assume m_synthesisWindowSize always = fftsize
    // add awindow's area (scalar) to swindow (vector) and copy the values to windowAccumulator
    // FIXME: can we get rid of the 1.5 magic number here?
    m_synthesisWindow->add2dst(windowAccumulator, m_analyzeWindow->GetArea() * 1.5);


}

// void
// phasevocodercore::Impl::synthesiseSliceCarrierRealImag(size_t channel,size_t shiftIncrement) {
//     channelinfo &caad = *m_carrierData[channel];

//     data_type *const R__ internalbuffer = caad.internalbuffer;
//     float *const R__ interfacebuffer = caad.interfacebuffer;
//     float *const R__ outputAccumulator = caad.outputAccumulator;
//     float *const R__ windowAccumulator = caad.windowAccumulator;
    
//     const int halfsize = m_fftSize / 2;

//     float factor = 1.f / m_fftSize;
    
//     vector_gain(caad.real, factor, halfsize + 1);
//     vector_gain(caad.imag, factor, halfsize + 1);

//     caad.fft->inverse(caad.real, caad.imag, caad.internalbuffer); // tangkk:: ifft here

//     ifftshift(interfacebuffer, m_fftSize, internalbuffer, m_synthesisWindow);
//     // vector_assign(interfacebuffer, internalbuffer + halfsize, halfsize);
//     // vector_assign(interfacebuffer + halfsize, internalbuffer, halfsize);
//     // m_synthesisWindow->apply(interfacebuffer);

//     vector_add(outputAccumulator, interfacebuffer, m_synthesisWindowSize);
//     // caad.accumulatorFill = wsz;

//     // we assume m_synthesisWindowSize always = fftsize
//     m_synthesisWindow->add2dst(windowAccumulator, m_analyzeWindow->GetArea()); // add awindow's area (scalar) to swindow (vector) and copy the values to windowAccumulator

// }

int
phasevocodercore::Impl::writeSlice(size_t channel, size_t shiftIncrement, bool last) {

    channelinfo &ad = *m_audioData[channel];
    
    float *const R__ outputAccumulator = ad.outputAccumulator;
    float *const R__ windowAccumulator = ad.windowAccumulator;

    if (m_logLevel > 2) {
        cerr << "writeSlice(" << channel << ", " << shiftIncrement << ", " << last << ")" << endl;
    }

    vector_divide(outputAccumulator, windowAccumulator, shiftIncrement); // factor outputAccumulator by windowAccumulator

    size_t outframes = 0;
    if (m_pitchScale != 1.0 && ad.res) {

        size_t reqSize = int(ceil(shiftIncrement / m_pitchScale));

        if (reqSize > resamplebufSize) {
            // This shouldn't normally happen -- the buffer is
            // supposed to be initialised with enough space in the
            // first place.  But we retain this check in case the
            // pitch scale has changed since then, or the stretch
            // calculator has gone mad, or something.
            cerr << "WARNING: phasevocodercore::Impl::writeSlice: resizing resampler buffer from "
                      << resamplebufSize << " to " << reqSize << endl;
            SetResampleBufSize(reqSize);
        }

        // tangkk: this is where the resampling happens (pitchshift is implemented as time-stretch + resampling)
        outframes = ad.res->doresample(&ad.outputAccumulator,
                                                  &resamplebuf,
                                                  shiftIncrement,
                                                  1.0 / m_pitchScale,
                                                  last);

        writeOutSimple(*ad.outbuf, resamplebuf, outframes, ad.outcnt);

    } else {
        outframes = shiftIncrement;
        writeOutSimple(*ad.outbuf, outputAccumulator, shiftIncrement, ad.outcnt);
    }

    // data has been written, shift the outputAccumulator (and zeros the tail), for OLA the result of the next slice
    vector_move(outputAccumulator, outputAccumulator + shiftIncrement, m_synthesisWindowSize - shiftIncrement);
    vector_zeros(outputAccumulator + m_synthesisWindowSize - shiftIncrement, shiftIncrement);
    
    // also do the same for windowAccumulator, because it is also involved in the OLA process
    vector_move(windowAccumulator, windowAccumulator + shiftIncrement, m_synthesisWindowSize - shiftIncrement);
    vector_zeros(windowAccumulator + m_synthesisWindowSize - shiftIncrement, shiftIncrement);


    return outframes;
}

int
phasevocodercore::Impl::writeSliceCarrier(size_t channel, size_t shiftIncrement) {
    channelinfo &caad = *m_carrierData[channel];
    
    float *const R__ outputAccumulator = caad.outputAccumulator;
    float *const R__ windowAccumulator = caad.windowAccumulator;

    int required = shiftIncrement;
    int ws = caad.outbuf->GetWriteSpace();
    if (ws < required) {
        // if the outbuf is initialized properly, this should not be entered
        // circularqueue<float> *oldbuf = caad.outbuf;
        // caad.outbuf = oldbuf->resized(oldbuf->GetSize() + (required - ws));
        // m_emergencyScavenger.claim(oldbuf);
        cerr << "Buffer overrun on output for channel, return 0 " << endl;
        return 0;
    }

    // size_t theoreticalOut = caad.inputSize;
    vector_divide(outputAccumulator, windowAccumulator, shiftIncrement);

    size_t outframes = 0;
    writeOutSimple(*caad.outbuf, outputAccumulator, shiftIncrement, caad.outcnt);
    outframes = shiftIncrement;

    // data has been written, shift the outputAccumulator (and zeros the tail), for OLA the result of the next slice
    vector_move(outputAccumulator, outputAccumulator + shiftIncrement, m_synthesisWindowSize - shiftIncrement);
    vector_zeros(outputAccumulator + m_synthesisWindowSize - shiftIncrement, shiftIncrement);
    
    // also do the same for windowAccumulator, because it is also involved in the OLA process
    vector_move(windowAccumulator, windowAccumulator + shiftIncrement, m_synthesisWindowSize - shiftIncrement);
    vector_zeros(windowAccumulator + m_synthesisWindowSize - shiftIncrement, shiftIncrement);

    return outframes;

}

void
phasevocodercore::Impl::writeOutSimple(circularqueue<float> &to, float *from, size_t qty, size_t &outcnt) {
    size_t written = to.write(from, qty);
    outcnt += written;
    return;
}

int
phasevocodercore::Impl::numsamples_available() const {
    size_t ret = 0;

    for (size_t i = 0; i < m_channels; ++i) {
        size_t availOut = m_audioData[i]->outbuf->GetReadSpace();
        if (i == 0 || availOut < ret) ret = availOut; // tangkk: this must enter once
    }

    return ret;

}

int
phasevocodercore::Impl::numsamples_availableCarrier() const {
    size_t ret = 0;

    for (size_t i = 0; i < m_channels; ++i) {
        size_t availOut = m_carrierData[i]->outbuf->GetReadSpace();
        if (i == 0 || availOut < ret) ret = availOut; // tangkk: this must enter once
    }

    return ret;

}

size_t
phasevocodercore::Impl::retrieve(float *const *output, size_t num_samples) const {
    size_t ret = num_samples;

    // tangkk: this is the normal case (normal left and right)
    for (size_t chan = 0; chan < m_channels; ++chan) {
        size_t thisret = m_audioData[chan]->outbuf->read(output[chan], ret);
        if (thisret < ret) {
            if (chan > 0) {
                if (m_logLevel > 0) {
                    cerr << "phasevocodercore::Impl::retrieve: WARNING: channel imbalance detected" << endl;
                }
            }
            ret = thisret; // balance them by the shorter one
        }
    }

    return ret;
}

size_t
phasevocodercore::Impl::retrieveCarrier(float *const *output, size_t num_samples) const {
    size_t ret = num_samples;

    // tangkk: this is the normal case (normal left and right)
    for (size_t chan = 0; chan < m_channels; ++chan) {
        size_t thisret = m_carrierData[chan]->outbuf->read(output[chan], ret);
        if (thisret < ret) {
            if (chan > 0) {
                if (m_logLevel > 0) {
                    cerr << "phasevocodercore::Impl::retrieve: WARNING: channel imbalance detected" << endl;
                }
            }
            ret = thisret; // balance them by the shorter one
        }
    }

    return ret;
}

} // end of namespace audiomod
