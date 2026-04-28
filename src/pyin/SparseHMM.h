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

#ifndef _SPARSEHMM_H_
#define _SPARSEHMM_H_

#include <vector>
#include <cstdio>

using std::vector;
using std::pair;

class SparseHMM
{
public:
    virtual const std::vector<float> calculateObsProb(const vector<pair<float, float> >);
    const std::vector<int> decodeViterbi(std::vector<vector<float> > obs, 
                                   vector<float> *scale);
    vector<float> init;
    vector<size_t> from;
    vector<size_t> to;
    vector<float> transProb;
};

#endif
