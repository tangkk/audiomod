#ifndef QBSH_FPRINTER_GEN__H
#define QBSH_FPRINTER_GEN__H

#include "PYinVamp.h"

struct pYINInfo
{
	float	iStart;
	float	iDuration;
	float	iNote;
};

class pYINAnalyzer
{
public:
	pYINAnalyzer(int sr = 44100, int blockSize = 1024, int stepSize = 256, float onsetSens = 0.7);
	~pYINAnalyzer();

	void InitpYINAnalyzer(int sr);

	// process the whole piece of audio (all samples)
	std::vector<pYINInfo> pYINProcess(const float samples[], int nSamples, bool flushFlag);

private:
	// process one block of samples in m_blockBuf
	void pYINProcessBlock(bool flushFlag, int sr);
private:
	PYinVamp *m_pyin;
	float* m_blockBuf;

	int m_LeftSampleNum;
	int m_frame_count;

	std::vector<float> m_LeftSamples;
	std::vector<pYINInfo> m_pYINNotes;

	int m_SampleRate;
	int m_BlockSize;
	int m_StepSize;
	float m_OnsetSens;
};

#endif
