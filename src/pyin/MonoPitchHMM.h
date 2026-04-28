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

#ifndef _MONOPITCHHMM_H_
#define _MONOPITCHHMM_H_

#include "SparseHMM.h"

#include <boost/math/distributions.hpp>

#include <vector>
#include <cstdio>

using std::vector;

class MonoPitchHMM : public SparseHMM
{
public:
    MonoPitchHMM();
    const std::vector<float> calculateObsProb(const vector<pair<float, float> >);
    // float getMidiPitch(size_t index);
    // float getFrequency(size_t index);
    void build();
    float m_minFreq; // 82.40689f/2
    size_t m_nBPS;
    size_t m_nPitch;
    size_t m_transitionWidth;
    float m_selfTrans;
    float m_yinTrust;
    vector<float> m_freqs;
};

#endif
