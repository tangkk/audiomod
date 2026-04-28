#include "PYinVamp.h"
#include "pYINAnalyzer.h"

pYINAnalyzer::pYINAnalyzer(int sr, int blockSize, int stepSize, float onsetSens)
{
	m_BlockSize = blockSize;
	m_StepSize = stepSize;
	m_SampleRate = sr;
	m_OnsetSens = onsetSens;
	std::cerr << "m_BlockSize:" << m_BlockSize << ", m_StepSize:" << m_StepSize << ", m_SampleRate:" << m_SampleRate << ", m_OnsetSens:" << m_OnsetSens << std::endl;
	m_LeftSamples.resize(m_BlockSize);
	InitpYINAnalyzer(sr);
}

pYINAnalyzer::~pYINAnalyzer()
{
	delete[] m_blockBuf;
	m_blockBuf = NULL;

	delete m_pyin;
}

void pYINAnalyzer::InitpYINAnalyzer(int sr)
{
	m_LeftSampleNum = 0;
	m_frame_count = 0;
	m_blockBuf = new float[m_BlockSize];
    m_pyin = new PYinVamp(sr);
	m_pyin->initialise(1, m_StepSize, m_BlockSize);
	m_pyin->getOutputDescriptors();
	m_pyin->setParameter("onsetsensitivity", m_OnsetSens);
}



std::vector<pYINInfo> pYINAnalyzer::pYINProcess(const float samples[], int nSamples, bool flushFlag){
	int totalSamplesNum = m_LeftSampleNum + nSamples;
	if (totalSamplesNum < m_BlockSize)
	{
		for (int i = 0; i < nSamples; ++i)
			m_LeftSamples[m_LeftSampleNum + i] = samples[i];

		m_LeftSampleNum = totalSamplesNum;
	}
	else{
		vector<float> vSamples(totalSamplesNum);
		std::cerr << "m_LeftSampleNum:" << m_LeftSampleNum << std::endl;
		for (int i = 0; i < m_LeftSampleNum; ++i)
			vSamples[i] = m_LeftSamples[i];

        for (int i = 0; i < nSamples; ++i)
			vSamples[m_LeftSampleNum + i] = samples[i];

		int FrameNum = (totalSamplesNum - m_BlockSize) / m_StepSize + 1;
		for (int i = 0; i < FrameNum; i++)
		{
			for (int j = 0; j < m_BlockSize; j++) {
                m_blockBuf[j] = vSamples[i* m_StepSize + j];
				//m_blockBuf[j] = vSamples[i* m_StepSize + j] / 32768.0;
				// printf("%f", m_blockBuf[j]);
			}
			pYINProcessBlock(flushFlag, m_SampleRate);
			m_frame_count++;

		}

		m_LeftSampleNum = totalSamplesNum - FrameNum*m_StepSize;
		for (int i = 0, j = FrameNum*m_StepSize; i < m_LeftSampleNum; ++i, ++j)
			m_LeftSamples[i] = vSamples[j];
	}

	if (flushFlag)
	{
		PYinVamp::FeatureSet fs = m_pyin->getRemainingFeatures();
		std::cerr << "fs[5].size():" << fs[5].size() << std::endl;

		for (int i = 0; i < fs[5].size(); i++)
		{
			pYINInfo tmpNote;
			tmpNote.iNote = fs[5][i].values[0];
			//tmpNote.iNote = 12 * std::log(fs[5][i].values[0] / 440) / std::log(2.) + 69.5;
			tmpNote.iStart = fs[5][i].timestamp.nsec / 1000.0;
			tmpNote.iDuration = fs[5][i].duration.nsec / 1000.0;
			m_pYINNotes.push_back(tmpNote);
		}
	}

	return m_pYINNotes;
}

void pYINAnalyzer::pYINProcessBlock(bool flushFlag, int sr){

	Vamp::RealTime tm;
    tm.nsec = (int)(((float(m_frame_count*m_StepSize + m_BlockSize / 2) * 1000) / sr));

	m_pyin->process(m_blockBuf, tm);
}
