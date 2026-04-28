/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    pYIN - A fundamental frequency estimator for monophonic audio
    Centre for Digital Music, Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "PYinVamp.h"
#include "MonoNote.h"
#include "MonoPitch.h"

#include "FFT.h"

#include <vector>
#include <algorithm>

#include <cstdio>
#include <cmath>
#include <complex>
#include <math.h>

using std::string;
using std::vector;
using Vamp::RealTime;

PYinVamp::PYinVamp(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_stepSize(256),
    m_blockSize(2048),
    m_fmin(50),		// m_fmin(40),// orig
    m_fmax(1000), //m_fmax(1600), // orig
    m_yin(2048, inputSampleRate, 0.0),
    m_oF0Candidates(0),
    m_oF0Probs(0),
    m_oVoicedProb(0),
    m_oCandidateSalience(0),
    m_oSmoothedPitchTrack(0),
    m_oNotes(0),
    m_threshDistr(2.0f),
    m_outputUnvoiced(0.0f),
    m_preciseTime(0.0f),
	m_lowAmp(0.1f),
    //m_lowAmp(0.1f), // 0.1 orig
    //m_onsetSensitivity(0.7f),// m_onsetSensitivity(0.7f),//orig, 0.2 better than 0.7
    //m_pruneThresh(0.1f),//m_pruneThresh(0.1f),// orig,0.1:100ms, 0.5:50ms
	m_onsetSensitivity(0.7f),// m_onsetSensitivity(0.7f),//orig, 0.2 better than 0.7
    m_pruneThresh(0.1f),//m_pruneThresh(0.1f),// orig,0.1:100ms, 0.5:50ms
    m_pitchProb(0),
    m_timestamp(0),
    m_level(0)
{
}

PYinVamp::~PYinVamp()
{
}

string
PYinVamp::getIdentifier() const
{
    return "pyin";
}

string
PYinVamp::getName() const
{
    return "pYin";
}

string
PYinVamp::getDescription() const
{
    return "Monophonic pitch and note tracking based on a probabilistic Yin extension.";
}

string
PYinVamp::getMaker() const
{
    return "Matthias Mauch";
}

int
PYinVamp::getPluginVersion() const
{
    // Increment this each time you release a version that behaves
    // differently from the previous one
    return 2;
}

string
PYinVamp::getCopyright() const
{
    return "GPL";
}

PYinVamp::InputDomain
PYinVamp::getInputDomain() const
{
    return TimeDomain;
}

size_t
PYinVamp::getPreferredBlockSize() const
{
    return 2048;
}

size_t 
PYinVamp::getPreferredStepSize() const
{
    return 256;
}

size_t
PYinVamp::getMinChannelCount() const
{
    return 1;
}

size_t
PYinVamp::getMaxChannelCount() const
{
    return 1;
}

PYinVamp::ParameterList
PYinVamp::getParameterDescriptors() const
{
    ParameterList list;
    
    ParameterDescriptor d;

    d.identifier = "threshdistr";
    d.name = "Yin threshold distribution";
    d.description = ".";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 7.0f;
    d.defaultValue = 2.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    d.valueNames.push_back("Uniform");
    d.valueNames.push_back("Beta (mean 0.10)");
    d.valueNames.push_back("Beta (mean 0.15)");
    d.valueNames.push_back("Beta (mean 0.20)");
    d.valueNames.push_back("Beta (mean 0.30)");
    d.valueNames.push_back("Single Value 0.10");
    d.valueNames.push_back("Single Value 0.15");
    d.valueNames.push_back("Single Value 0.20");
    list.push_back(d);

    d.identifier = "outputunvoiced";
    d.valueNames.clear();
    d.name = "Output estimates classified as unvoiced?";
    d.description = ".";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 2.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    d.valueNames.push_back("No");
    d.valueNames.push_back("Yes");
    d.valueNames.push_back("Yes, as negative frequencies");
    list.push_back(d);

    d.identifier = "precisetime";
    d.valueNames.clear();
    d.name = "Use non-standard precise YIN timing (slow).";
    d.description = ".";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    list.push_back(d);

    d.identifier = "lowampsuppression";
    d.valueNames.clear();
    d.name = "Suppress low amplitude pitch estimates.";
    d.description = ".";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.1f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "onsetsensitivity";
    d.valueNames.clear();
    d.name = "Onset sensitivity";
    d.description = "Adds additional note onsets when RMS increases.";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.7f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "prunethresh";
    d.valueNames.clear();
    d.name = "Duration pruning threshold.";
    d.description = "Prune notes that are shorter than this value.";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 0.2f;
    d.defaultValue = 0.1f;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float
PYinVamp::getParameter(string identifier) const
{
    if (identifier == "threshdistr") {
            return m_threshDistr;
    }
    if (identifier == "outputunvoiced") {
            return m_outputUnvoiced;
    }
    if (identifier == "precisetime") {
            return m_preciseTime;
    }
    if (identifier == "lowampsuppression") {
            return m_lowAmp;
    }
    if (identifier == "onsetsensitivity") {
            return m_onsetSensitivity;
    }
    if (identifier == "prunethresh") {
            return m_pruneThresh;
    }
    return 0.f;
}

void
PYinVamp::setParameter(string identifier, float value) 
{
    if (identifier == "threshdistr")
    {
        m_threshDistr = value;
    }
    if (identifier == "outputunvoiced")
    {
        m_outputUnvoiced = value;
    }
    if (identifier == "precisetime")
    {
        m_preciseTime = value;
    }
    if (identifier == "lowampsuppression")
    {
        m_lowAmp = value;
    }
    if (identifier == "onsetsensitivity")
    {
        m_onsetSensitivity = value;
    }
    if (identifier == "prunethresh")
    {
        m_pruneThresh = value;
    }
}

PYinVamp::ProgramList
PYinVamp::getPrograms() const
{
    ProgramList list;
    return list;
}

string
PYinVamp::getCurrentProgram() const
{
    return ""; // no programs
}

void
PYinVamp::selectProgram(string name)
{
}

PYinVamp::OutputList
PYinVamp::getOutputDescriptors() const
{
    OutputList outputs;

    OutputDescriptor d;
    
    int outputNumber = 0;

    d.identifier = "f0candidates";
    d.name = "F0 Candidates";
    d.description = "Estimated fundamental frequency candidates.";
    d.unit = "Hz";
    d.hasFixedBinCount = false;
    // d.binCount = 1;
    d.hasKnownExtents = true;
    d.minValue = m_fmin;
    d.maxValue = 500;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oF0Candidates = outputNumber++;

    d.identifier = "f0probs";
    d.name = "Candidate Probabilities";
    d.description = "Probabilities  of estimated fundamental frequency candidates.";
    d.unit = "";
    d.hasFixedBinCount = false;
    // d.binCount = 1;
    d.hasKnownExtents = true;
    d.minValue = 0;
    d.maxValue = 1;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oF0Probs = outputNumber++;
    
    d.identifier = "voicedprob";
    d.name = "Voiced Probability";
    d.description = "Probability that the signal is voiced according to Probabilistic Yin.";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = true;
    d.minValue = 0;
    d.maxValue = 1;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oVoicedProb = outputNumber++;

    d.identifier = "candidatesalience";
    d.name = "Candidate Salience";
    d.description = "Candidate Salience";
    d.hasFixedBinCount = true;
    d.binCount = m_blockSize / 2;
    d.hasKnownExtents = true;
    d.minValue = 0;
    d.maxValue = 1;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oCandidateSalience = outputNumber++;
    
    d.identifier = "smoothedpitchtrack";
    d.name = "Smoothed Pitch Track";
    d.description = ".";
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    // d.minValue = 0;
    // d.maxValue = 1;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oSmoothedPitchTrack = outputNumber++;

    d.identifier = "notes";
    d.name = "Notes";
    d.description = "Derived fixed-pitch note frequencies";
    // d.unit = "MIDI unit";
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = true;
    outputs.push_back(d);
    m_oNotes = outputNumber++;

    return outputs;
}

bool
PYinVamp::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() || channels > getMaxChannelCount()) return false;

/*
    std::cerr << "PYinVamp::initialise: channels = " << channels
          << ", stepSize = " << stepSize << ", blockSize = " << blockSize
          << std::endl;
*/
    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;
    
    reset();

    return true;
}

void
PYinVamp::reset()
{    
    m_yin.setThresholdDistr(m_threshDistr);
    m_yin.setFrameSize(m_blockSize);
    m_yin.setFast(!m_preciseTime);
    
    m_pitchProb.clear();
    m_timestamp.clear();
    m_level.clear();
/*    
    std::cerr << "PYinVamp::reset"
          << ", blockSize = " << m_blockSize
          << std::endl;
*/
}

PYinVamp::FeatureSet
PYinVamp::process(const float *inputBuffers, RealTime timestamp)	//const float *const *inputBuffers
{
#if 0
    int offset = m_preciseTime == 1.0 ? m_blockSize/2 : m_blockSize/4;
    //timestamp = timestamp + Vamp::RealTime::frame2RealTime(offset, (int)m_inputSampleRate); //lrintf(m_inputSampleRate)
	RealTime tmpTime = Vamp::RealTime::frame2RealTime(offset, (int)m_inputSampleRate); 
	timestamp.nsec += tmpTime.nsec;

	//printf("%d, %d, %d, %d, %f\n", tmpTime.nsec, timestamp.nsec/1e6, offset, m_blockSize, m_inputSampleRate);
#endif
    FeatureSet fs;
    
    float rms = 0;
    
	//printf("point 0, m_blockSize = %d\n", m_blockSize);

    float *dInputBuffers = new float[m_blockSize];
    for (size_t i = 0; i < m_blockSize; ++i) {
        dInputBuffers[i] = inputBuffers[i];
        rms += inputBuffers[i] * inputBuffers[i];
    }
    rms /= m_blockSize;
    rms = sqrt(rms);
    
    bool isLowAmplitude = (rms < m_lowAmp);

	//printf("point 1\n");
	//printf("m_pitchProb.size = %d, rms = %f %f\n", m_pitchProb.size(), rms, m_lowAmp);
    
    Yin::YinOutput yo = m_yin.processProbabilisticYin(dInputBuffers);
    delete [] dInputBuffers;

	//printf("point 2\n");
    m_level.push_back(yo.rms);

    // First, get the things out of the way that we don't want to output 
    // immediately, but instead save for later.
    vector<pair<float, float> > tempPitchProb;
    for (size_t iCandidate = 0; iCandidate < yo.freqProb.size(); ++iCandidate)
    {
        float tempPitch = 12 * std::log(yo.freqProb[iCandidate].first/440)/std::log(2.) + 69;
        if (!isLowAmplitude)
        {
            tempPitchProb.push_back(pair<float, float>
                (tempPitch, yo.freqProb[iCandidate].second));
        } else {
            float factor = ((rms+0.01*m_lowAmp)/(1.01*m_lowAmp));
            tempPitchProb.push_back(pair<float, float>
                (tempPitch, yo.freqProb[iCandidate].second*factor));
        }
    }
    m_pitchProb.push_back(tempPitchProb);
    m_timestamp.push_back(timestamp);

#if 0
    // F0 CANDIDATES
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = timestamp;
    for (size_t i = 0; i < yo.freqProb.size(); ++i)
    {
        f.values.push_back(yo.freqProb[i].first);
    }
    fs[m_oF0Candidates].push_back(f);
    
    // CANDIDATES PROBABILITIES
    f.values.clear();
    float voicedProb = 0;
    for (size_t i = 0; i < yo.freqProb.size(); ++i)
    {
        f.values.push_back(yo.freqProb[i].second);
        voicedProb += yo.freqProb[i].second;
    }
    fs[m_oF0Probs].push_back(f);
    
	// VOICEDPROB
	// f.values.clear();
    f.values.push_back(voicedProb);
    fs[m_oVoicedProb].push_back(f);

    // SALIENCE -- maybe this should eventually disappear
    f.values.clear();
    float salienceSum = 0;
    for (size_t iBin = 0; iBin < yo.salience.size(); ++iBin)
    {
        f.values.push_back(yo.salience[iBin]);
        salienceSum += yo.salience[iBin];
    }
    fs[m_oCandidateSalience].push_back(f);
#endif
    return fs;
}

#if 1
PYinVamp::FeatureSet
PYinVamp::getRemainingFeatures()
{
    FeatureSet fs;
    Feature f;
    f.hasTimestamp = true;
    f.hasDuration = false;
    
    if (m_pitchProb.empty()) {
        return fs;
    }

    // MONO-PITCH STUFF
    MonoPitch mp;
    vector<float> mpOut = mp.process(m_pitchProb);
	for (size_t iFrame = 0; iFrame < mpOut.size(); ++iFrame)
	{
		if (mpOut[iFrame] < 0 && (m_outputUnvoiced==0)) continue;
		f.timestamp = m_timestamp[iFrame];
		f.values.clear();
        if (m_outputUnvoiced == 1)
        {
            f.values.push_back(fabs(mpOut[iFrame]));
        } else {
            f.values.push_back(mpOut[iFrame]);
        }
        
        fs[m_oSmoothedPitchTrack].push_back(f);
    }
    
	// printf("m_pitchProb mpOut = %d, %d, %d\n", m_pitchProb.size(), mpOut.size(), m_oSmoothedPitchTrack);

    // MONO-NOTE STUFF
	// std::cerr << "Mono Note Stuff" << std::endl;
	MonoNote mn;
	std::vector<std::vector<std::pair<float, float> > > smoothedPitch;
    for (size_t iFrame = 0; iFrame < mpOut.size(); ++iFrame) {
        std::vector<std::pair<float, float> > temp;
        if (mpOut[iFrame] > 0)
        {
            float tempPitch = 12 * std::log(mpOut[iFrame]/440)/std::log(2.) + 69;
            temp.push_back(std::pair<float,float>(tempPitch, .9));
        }
        smoothedPitch.push_back(temp);
    }
    // vector<MonoNote::FrameOutput> mnOut = mn.process(m_pitchProb);
    vector<MonoNote::FrameOutput> mnOut = mn.process(smoothedPitch);

	//printf("m_pitchProb mpOut mnOut = %d, %d, %d, %d, %d\n", m_pitchProb.size(), mpOut.size(), m_oSmoothedPitchTrack, mnOut.size(), m_oNotes);
    
    // turning feature into a note feature
    f.hasTimestamp = true;
    f.hasDuration = true;
    f.values.clear();
        
    int onsetFrame = 0;
    bool isVoiced = 0;
    bool oldIsVoiced = 0;
    size_t nFrame = m_pitchProb.size();

    float minNoteFrames = (m_inputSampleRate*m_pruneThresh) / m_stepSize;
    
    std::vector<float> notePitchTrack; // collects pitches for one note at a time
    for (size_t iFrame = 0; iFrame < nFrame; ++iFrame)
    {
        isVoiced = mnOut[iFrame].noteState < 3
                   && smoothedPitch[iFrame].size() > 0
                   && (iFrame >= nFrame-2
				    //|| ((m_level[iFrame]/m_level[iFrame+3]) > m_onsetSensitivity));		// rivsed code
                    || ((m_level[iFrame]/m_level[iFrame+2]) > m_onsetSensitivity)); // orig code
        // std::cerr << m_level[iFrame]/m_level[iFrame-1] << " " << isVoiced << std::endl;
        

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
                    f.values.clear();
                    f.values.push_back(medianFreq);
                    f.timestamp = m_timestamp[onsetFrame];
                    f.duration = m_timestamp[iFrame] - m_timestamp[onsetFrame];
					//printf("%d, %d, %d\n", m_oNotes, f.timestamp.nsec, f.duration.nsec);
                    fs[m_oNotes].push_back(f);
                }
                notePitchTrack.clear();
            }
        }
        oldIsVoiced = isVoiced;
	}

	return fs;
}
#else
PYinVamp::FeatureSet
PYinVamp::getRemainingFeatures()
{
	vector<Feature> noteRecordLess, noteRecordMore;
	vector<int> noteFrameLess, noteFrameMore;
	FeatureSet fs;
    Feature f;
    f.hasTimestamp = true;
    f.hasDuration = false;
    
    if (m_pitchProb.empty()) {
        return fs;
    }

    // MONO-PITCH STUFF
    MonoPitch mp;
    vector<float> mpOut = mp.process(m_pitchProb);
	for (size_t iFrame = 0; iFrame < mpOut.size(); ++iFrame)
	{
		if (mpOut[iFrame] < 0 && (m_outputUnvoiced==0)) continue;
		f.timestamp = m_timestamp[iFrame];
		f.values.clear();
        if (m_outputUnvoiced == 1)
        {
            f.values.push_back(fabs(mpOut[iFrame]));
        } else {
            f.values.push_back(mpOut[iFrame]);
        }
        
        fs[m_oSmoothedPitchTrack].push_back(f);
    }
    
    // MONO-NOTE STUFF
	// std::cerr << "Mono Note Stuff" << std::endl;
	MonoNote mn;
	std::vector<std::vector<std::pair<float, float> > > smoothedPitch;
    for (size_t iFrame = 0; iFrame < mpOut.size(); ++iFrame) {
        std::vector<std::pair<float, float> > temp;
        if (mpOut[iFrame] > 0)
        {
            float tempPitch = 12 * std::log(mpOut[iFrame]/440)/std::log(2.) + 69;
            temp.push_back(std::pair<float,float>(tempPitch, .9));
        }
        smoothedPitch.push_back(temp);
    }
    // vector<MonoNote::FrameOutput> mnOut = mn.process(m_pitchProb);
    vector<MonoNote::FrameOutput> mnOut = mn.process(smoothedPitch);
    
    // turning feature into a note feature
    f.hasTimestamp = true;
    f.hasDuration = true;
    f.values.clear();
        
    int onsetFrame = 0;
    bool isVoiced = 0;
    bool oldIsVoiced = 0;
    size_t nFrame = m_pitchProb.size();

    float minNoteFrames = (m_inputSampleRate*m_pruneThresh) / m_stepSize;
    
    std::vector<float> notePitchTrack; // collects pitches for one note at a time
    for (size_t iFrame = 0; iFrame < nFrame; ++iFrame)
    {
        isVoiced = mnOut[iFrame].noteState < 3
                   && smoothedPitch[iFrame].size() > 0
                   && (iFrame >= nFrame-2
				    //|| ((m_level[iFrame]/m_level[iFrame+3]) > m_onsetSensitivity));		// rivsed code
                    || ((m_level[iFrame]/m_level[iFrame+2]) > m_onsetSensitivity)); // orig code
        // std::cerr << m_level[iFrame]/m_level[iFrame-1] << " " << isVoiced << std::endl;
        

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
				int notePitchTrack_size = notePitchTrack.size();
                if (notePitchTrack_size >= minNoteFrames)
                {
                    std::sort(notePitchTrack.begin(), notePitchTrack.end());
                    float medianPitch = notePitchTrack[notePitchTrack.size()/2];
                    float medianFreq = std::pow(2,(medianPitch - 69) / 12) * 440;
                    f.values.clear();
                    f.values.push_back(medianFreq);
                    f.timestamp = m_timestamp[onsetFrame];
                    f.duration = m_timestamp[iFrame] - m_timestamp[onsetFrame];

					noteRecordMore.push_back(f);
					noteFrameMore.push_back(notePitchTrack_size);
					if (notePitchTrack_size >= minNoteFrames*2)
					{
						noteRecordLess.push_back(f);
						noteFrameLess.push_back(notePitchTrack_size);
					}
                    //fs[m_oNotes].push_back(f);
                }
                notePitchTrack.clear();
            }
        }
        oldIsVoiced = isVoiced;
	}

	// select thresold according median(note_duaration), select small thresold
	 std::sort(noteFrameMore.begin(), noteFrameMore.end());
     int medianDuration;
	 if (noteFrameMore.size() >= 1)
		medianDuration = noteFrameMore[noteFrameMore.size()/2];
	 else
		 medianDuration = minNoteFrames*5;

	 if ( medianDuration < minNoteFrames*3 )
	 //if ( medianDuration < minNoteFrames*0 )
	 {
		 for(int i = 0; i < noteRecordMore.size(); i++)
			 fs[m_oNotes].push_back(noteRecordMore[i]);
	 }
	 else
	 {
		 for(int i = 0; i < noteRecordLess.size(); i++)
			 fs[m_oNotes].push_back(noteRecordLess[i]);
	 }

	return fs;
}
#endif

PYinVamp::FeatureSet
PYinVamp::getRemainingFeatures0(vector<int> onsetMsAll)
	{
    FeatureSet fs;
    Feature f;
    f.hasTimestamp = true;
    f.hasDuration = false;
    
    if (m_pitchProb.empty()) {
        return fs;
    }

    // MONO-PITCH STUFF
    MonoPitch mp;
    vector<float> mpOut = mp.process(m_pitchProb);
	for (size_t iFrame = 0; iFrame < mpOut.size(); ++iFrame)
	{
		if (mpOut[iFrame] < 0 && (m_outputUnvoiced==0)) continue;
		f.timestamp = m_timestamp[iFrame];
		f.values.clear();
        if (m_outputUnvoiced == 1)
        {
            f.values.push_back(fabs(mpOut[iFrame]));
        } else {
            f.values.push_back(mpOut[iFrame]);
        }
        
        fs[m_oSmoothedPitchTrack].push_back(f);
    }
    
    // MONO-NOTE STUFF
	// std::cerr << "Mono Note Stuff" << std::endl;
	MonoNote mn;
	std::vector<std::vector<std::pair<float, float> > > smoothedPitch;
    for (size_t iFrame = 0; iFrame < mpOut.size(); ++iFrame) {
        std::vector<std::pair<float, float> > temp;
        if (mpOut[iFrame] > 0)
        {
            float tempPitch = 12 * std::log(mpOut[iFrame]/440)/std::log(2.) + 69;
            temp.push_back(std::pair<float,float>(tempPitch, .9));
        }
        smoothedPitch.push_back(temp);
    }
    // vector<MonoNote::FrameOutput> mnOut = mn.process(m_pitchProb);
    vector<MonoNote::FrameOutput> mnOut = mn.process(smoothedPitch);
    
    // turning feature into a note feature
    f.hasTimestamp = true;
    f.hasDuration = true;
    f.values.clear();
        
    int onsetFrame = 0;
    bool isVoiced = 0;
    bool oldIsVoiced = 0;
    size_t nFrame = m_pitchProb.size();

    float minNoteFrames = (m_inputSampleRate*m_pruneThresh) / m_stepSize;
    
#if 1
	std::vector<Feature> noteRecord;
    std::vector<float> notePitchTrack; // collects pitches for one note at a time
    for (size_t iFrame = 0; iFrame < nFrame; ++iFrame)
    {
        isVoiced = mnOut[iFrame].noteState < 3
                   && smoothedPitch[iFrame].size() > 0
                   && (iFrame >= nFrame-2
                       || ((m_level[iFrame]/m_level[iFrame+2]) > m_onsetSensitivity));
        // std::cerr << m_level[iFrame]/m_level[iFrame-1] << " " << isVoiced << std::endl;
        

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
                    f.values.clear();
                    f.values.push_back(medianFreq);
                    f.timestamp = m_timestamp[onsetFrame];
                    f.duration = m_timestamp[iFrame] - m_timestamp[onsetFrame];
					//printf("%d, %d, %d\n", m_oNotes, f.timestamp.nsec, f.duration.nsec);
                    //fs[m_oNotes].push_back(f); // orig
					if (f.duration.nsec > 50)
						noteRecord.push_back(f);
                }
                notePitchTrack.clear();
            }
        }
        oldIsVoiced = isVoiced;
	}
	
#if 1
	// note segment using onset detection results
	std::vector<Feature> noteRecordNew;
	std::vector<bool> noteInd;
	for(int i = 0; i < noteRecord.size()-1; i++)
	{
		// note < 150ms, merge with other note
#if 0
		if ( (noteRecord[i].duration.nsec < 50) && 
			( (noteRecord[i+1].timestamp.nsec-noteRecord[i].timestamp.nsec-noteRecord[i].duration.nsec) < 20) )
			noteInd.push_back(false);
		else
#endif
			noteInd.push_back(true);
	}
	noteInd.push_back(true);

	// merger note
	bool mergerFlag = false;
	for(int i = 0; i < noteInd.size()-1; )
	{
		if (true == noteInd[i])
		{
			noteRecordNew.push_back( noteRecord[i] );
			mergerFlag = false;
			i++;
		}
		else
		{
			// merge nerghboring note
			if (false == mergerFlag)
			{
				Feature tmp;
				tmp.values.push_back( (noteRecord[i].values[0] + noteRecord[i].values[0])/2 );
				tmp.timestamp = noteRecord[i].timestamp;
				tmp.duration = ( noteRecord[i].duration + noteRecord[i+1].duration );
				mergerFlag = true;
				i += 2;
			}
			// avoid re-merge note
			else
			{
				noteRecordNew.push_back( noteRecord[i] );
				mergerFlag = false;
				i++;
			}
		}
	}
	if (false == mergerFlag)
		noteRecordNew.push_back( noteRecord[noteRecord.size()-1] );

	for (int i = 0; i < noteRecordNew.size(); i++)
	{
		// left&right shrink 50ms
		int t1 = noteRecordNew[i].timestamp.nsec + 200;
		int t2 = noteRecordNew[i].timestamp.nsec - 200 + noteRecordNew[i].duration.nsec;

		// note re-segment or-not
		vector<int> segmentOnset;
		for (int j = 0; (j < onsetMsAll.size()) && (t1 < t2); j++)
		{
			if ( (onsetMsAll[j] > t1 ) && (onsetMsAll[j] < t2 ) )
			{
				printf("t1, t2, onset = %d, %d, %d\n", t1, t2, onsetMsAll[j]);
				segmentOnset.push_back( onsetMsAll[j] );	// onsets between two note
			}
		}

		if (0 == segmentOnset.size())
			fs[m_oNotes].push_back(noteRecordNew[i]);
		else
		{
			// head
			Feature tmp;
			tmp.values.push_back( noteRecordNew[i].values[0] );			// note
			tmp.timestamp.nsec = noteRecordNew[i].timestamp.nsec;		// onset
			tmp.duration.nsec = segmentOnset[0] - tmp.timestamp.nsec;	// duration
			fs[m_oNotes].push_back(tmp);

			printf("0 onset = %d\n", tmp.timestamp.nsec);

			// middle
			for (int j = 0; j < segmentOnset.size()-1; j++)
			{
				tmp.values.clear();
				tmp.values.push_back( noteRecordNew[i].values[0] );		// note
				tmp.timestamp.nsec = segmentOnset[j];					// onset
				tmp.duration.nsec = segmentOnset[j+1] -segmentOnset[j];	// duration
				fs[m_oNotes].push_back(tmp);
				printf("1 onset = %d\n", tmp.timestamp.nsec);
			}

			// tail
			tmp.values.clear();
			int tmpLen = segmentOnset.size()-1;
			tmp.values.push_back( noteRecordNew[i].values[0] );			// note
			tmp.timestamp.nsec = segmentOnset[tmpLen];					// onset
			tmp.duration.nsec = noteRecordNew[i].duration.nsec -
				(segmentOnset[tmpLen] - noteRecordNew[i].timestamp.nsec);	// duration
			fs[m_oNotes].push_back(tmp);
			printf("2 onset = %d\n", tmp.timestamp.nsec);
		}
	}

#endif

#else
	std::vector<int> onsetFrameAll;
	int offsetFrame = 0;
	// onsetMs to onsetFrame
	for (size_t i = 0; i < onsetMsAll.size(); i++)
	{
		
		int timeStepMs = m_stepSize * 1000 / m_inputSampleRate;
		int tmpOnset = onsetMsAll[i] / (float)timeStepMs + 0.5;
		onsetFrameAll.push_back(tmpOnset);
		//printf("onsetFrameAll[%d] = %d, %d\n", i, tmpOnset, timeStepMs);
	}
	
	std::vector<float> notePitchTrack; // collects pitches for one note at a time
	for (size_t i = 0; i < onsetMsAll.size()-1; i++)
	{
		if ( onsetFrameAll[i+1] < mpOut.size())
		{
			onsetFrame = onsetFrameAll[i];
			offsetFrame = onsetFrameAll[i+1];
			//printf("<%d>onsetFrame/offsetFrame = %d, %d, %d\n", i, onsetFrame, offsetFrame, offsetFrame-onsetFrame );

			// offset frame detection
			//for (size_t iFrame = (onsetFrame+offsetFrame)/2; (iFrame < offsetFrame+5) && (offsetFrame+4 < mpOut.size()); ++iFrame)
			for (size_t iFrame = (onsetFrame+offsetFrame)/2; iFrame <= offsetFrame; ++iFrame)
			{
				if (mpOut[iFrame] <= 0)
				{
					offsetFrame = iFrame;
					break;
				}
			}
			// onset frame detection
			/*
			for (size_t iFrame = (onsetFrame+offsetFrame)/2; (iFrame > onsetFrame-5) && (onsetFrame-4 >= 0); --iFrame)
			{
				if (mpOut[iFrame] <= 0)
				{
					onsetFrame = iFrame;
					break;
				}
			}*/

			//printf("<%d>onsetFrame/offsetFrame = %d, %d, %d\n", i, onsetFrame, offsetFrame, offsetFrame-onsetFrame );

			// pitchs in one note
			for (size_t iFrame = onsetFrame; iFrame <= offsetFrame; ++iFrame)
			{
				 if (mpOut[iFrame] > 0)
				 {
					 float tempPitch = 12 * std::log(mpOut[iFrame]/440)/std::log(2.) + 69;
					 notePitchTrack.push_back(tempPitch); // add to the note's pitch track
				 }
			}

			// midian pitch in one note
			//if (notePitchTrack.size() >= minNoteFrames)
			if (notePitchTrack.size() >= 2)
			{
				std::sort(notePitchTrack.begin(), notePitchTrack.end());
				float medianPitch = notePitchTrack[notePitchTrack.size()/2];
				float medianFreq = std::pow(2,(medianPitch - 69) / 12) * 440;
				if (medianFreq <= 90)
					medianFreq *=2.0;

				f.values.clear();
				f.values.push_back(medianFreq);
				f.timestamp = m_timestamp[onsetFrame];
				f.duration = m_timestamp[offsetFrame] - m_timestamp[onsetFrame];
				//printf("%d, %d, %d\n", m_oNotes, f.timestamp.nsec, f.duration.nsec);
				fs[m_oNotes].push_back(f);
			}
			notePitchTrack.clear();
		}
	}
#endif

    return fs;
}
