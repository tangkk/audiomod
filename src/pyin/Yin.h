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

#ifndef _YIN_H_
#define _YIN_H_

#include "FFT.h"
#include "MeanFilter.h"

#include <cmath>

#include <iostream>
#include <vector>
#include <exception>

using std::vector;
using std::pair;



class Yin
{
public:
    Yin(size_t frameSize, size_t inputSampleRate, float thresh = 0.2, bool fast = true);
    virtual ~Yin();

    struct YinOutput {
        float f0;
        float periodicity;
        float rms;
        vector<float> salience;
        vector<pair<float, float> > freqProb;
        YinOutput() :  f0(0), periodicity(0), rms(0), 
            salience(vector<float>(0)), freqProb(vector<pair<float, float> >(0)) { }
        YinOutput(float _f, float _p, float _r) :
            f0(_f), periodicity(_p), rms(_r), 
            salience(vector<float>(0)), freqProb(vector<pair<float, float> >(0)) { }
        YinOutput(float _f, float _p, float _r, vector<float> _salience) :
            f0(_f), periodicity(_p), rms(_r), salience(_salience), 
            freqProb(vector<pair<float, float> >(0)) { }
    };
    
    int setThreshold(float parameter);
    int setThresholdDistr(float parameter);
    int setFrameSize(size_t frameSize);
    int setFast(bool fast);
    // int setRemoveUnvoiced(bool frameSize);
    YinOutput process(const float *in) const;
    YinOutput processProbabilisticYin(const float *in) const;

private:
    mutable size_t m_frameSize;
    mutable size_t m_inputSampleRate;
    mutable float m_thresh;
    mutable size_t m_threshDistr;
    mutable size_t m_yinBufferSize;
    mutable bool   m_fast;
    // mutable bool m_removeUnvoiced;
};

#endif
