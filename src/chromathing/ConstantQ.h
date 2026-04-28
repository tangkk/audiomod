/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    QM DSP Library

    Centre for Digital Music, Queen Mary, University of London.
    This file 2005-2006 Christian Landone.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#ifndef QM_DSP_CONSTANTQ_H
#define QM_DSP_CONSTANTQ_H

#include <vector>
#include "MathAliases.h"
#include "MathUtilities.h"

struct CQConfig {
    float FS;         // samplerate
    float min;        // minimum frequency
    float max;        // maximum frequency
    int BPO;           // bins per octave
    float CQThresh;   // threshold
};

class ConstantQ
{
public:
    ConstantQ(CQConfig config);
    ~ConstantQ();

    void process(const float* FFTRe, const float* FFTIm,
                 float* CQRe, float* CQIm);

    float* process(const float* FFTData);

    void sparsekernel();

    float getQ() { return m_dQ; }
    int getK() { return m_uK; }
    int getFFTLength() { return m_FFTLength; }
    int getHop() { return m_hop; }

private:
    void initialise(CQConfig config);
    void deInitialise();
        
    float* m_CQdata;
    float m_FS;
    float m_FMin;
    float m_FMax;
    float m_dQ;
    float m_CQThresh;
    int m_hop;
    int m_BPO;
    int m_FFTLength;
    int m_uK;

    struct SparseKernel {
        std::vector<int> is;
        std::vector<int> js;
        std::vector<float> imag;
        std::vector<float> real;
    };

    SparseKernel *m_sparseKernel;
};


#endif//CONSTANTQ_H

