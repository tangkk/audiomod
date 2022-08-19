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

#include <cstddef>
#include <map>
#include <vector>

namespace audiomod {

class  phasevocodercore
{
public:
    
    enum Option {

        // OptionProcessOffline       = 0x00000000,
        // OptionProcessRealTime      = 0x00000001,

        // OptionStretchElastic       = 0x00000000,
        // OptionStretchPrecise       = 0x00000010,
    
        // OptionTransientsCrisp      = 0x00000000,
        // OptionTransientsMixed      = 0x00000100,
        // OptionTransientsSmooth     = 0x00000200,

        // OptionDetectorCompound     = 0x00000000,
        // OptionDetectorPercussive   = 0x00000400,
        // OptionDetectorSoft         = 0x00000800,

        // OptionPhaseLaminar         = 0x00000000,
        // OptionPhaseIndependent     = 0x00002000,
    
        // OptionThreadingAuto        = 0x00000000,
        // OptionThreadingNever       = 0x00010000,
        // OptionThreadingAlways      = 0x00020000,

        // OptionWindowStandard       = 0x00000000,
        // OptionWindowShort          = 0x00100000,
        // OptionWindowLong           = 0x00200000,

        // OptionSmoothingOff         = 0x00000000,
        // OptionSmoothingOn          = 0x00800000,

        OptionFormantShifted       = 0x00000000,
        OptionFormantPreserved     = 0x01000000,
        OptionGenderChange         = 0x02000000,
        OptionRobotic              = 0x04000000,
        OptionWhisper              = 0x08000000,

        // OptionPitchHighSpeed       = 0x00000000,
        // OptionPitchHighQuality     = 0x04000000,
        // OptionPitchHighConsistency = 0x08000000,

        // OptionChannelsApart        = 0x00000000,
        // OptionChannelsTogether     = 0x10000000,

        // n.b. Options is int, so we must stop before 0x80000000
    };

    typedef int Options;

    // enum PresetOption {
    //     DefaultOptions             = 0x00000000,
    //     PercussiveOptions          = 0x00102000
    // };


    phasevocodercore(size_t sampleRate,
                        size_t channels,
                        Options options = 0x00000000,
                        double initialTimeRatio = 1.0,
                        double initialPitchScale = 1.0);
    ~phasevocodercore();


    /**
     * Return the last time ratio value that was set (either on
     * construction or with setTimeRatio()).
     */
    double getTimeRatio() const;

    /**
     * Return the last pitch scaling ratio value that was set (either
     * on construction or with setPitchScale()).
     */
    double getPitchScale() const;

    // size_t getLatency() const;

    void processNormal(const float *const *input, size_t num_samples);
    
    /**
     * simplified version of process at a constant speed without stretch ratio. Used in realtime and single threaded
     */
    void processConstant(const float *const *input, size_t num_samples);

    /**
     * process input with a generated carrier to produce a vocoder effect
     */
    void processVocoder(const float *const *input, size_t num_samples, int carrierType);


    int numsamples_available() const;
    int numsamples_availableCarrier() const;


    size_t retrieve(float *const *output, size_t num_samples) const;
    size_t retrieveCarrier(float *const *output, size_t num_samples) const;

    /**
     * Return the value of internal frequency cutoff value n.
     *
     * This function is not for general use.
     */
    // float getFrequencyCutoff(int n) const;

    /** 
     * Set the value of internal frequency cutoff n to f Hz.
     *
     * This function is not for general use.
     */
    // void setFrequencyCutoff(int n, float f);
    
    /**
     * Retrieve the value of the internal input block hopsize value.
     *
     * This function is provided for diagnostic purposes only.
     */
    size_t getInputHopSize() const;


    /**
     * Return the number of channels this stretcher was constructed
     * with.
     */
    size_t getChannelCount() const;


    /**
     * Set the level of debug output.  The value may be from 0 (errors
     * only) to 3 (very verbose, with audible ticks in the output at
     * phase reset points).  The default is whatever has been set
     * using setDefaultDebugLevel, or 0 if that function has not been
     * called.
     */
    void setDebugLevel(int level);


    static void setDefaultDebugLevel(int level);

    static void setDefaultFftSize(int fftsize);
    static void setDefaultHopSize(int hopsize);
    static void setDefaultCoreMode(int coremode);

protected:
    class Impl;
    Impl *m_pv;
};

} // end of namespace audiomod
