// Reverb model implementation
//
// Written by Jezar at Dreampoint, June 2000
// http://www.dreampoint.co.uk
// This code is public domain

#include "revmodel.h"

revmodel::revmodel(float sampleRate)
{
	// Tie the components to their buffers
    for(int i = 0; i < numcombs; i++) {
        int sizeL = combtuningL[i] * sampleRate / 44100;
        bufcombLs[i] = new float[sizeL];
        combL[i].setbuffer(bufcombLs[i], sizeL);
        int sizeR = (combtuningL[i] + stereospread) * sampleRate / 44100;
        bufcombRs[i] = new float[sizeR];
        combR[i].setbuffer(bufcombRs[i], sizeR);
    }
    for(int i = 0; i < numallpasses; i++) {
        int sizeL = allpasstuningL[i] * sampleRate / 44100;
        bufallpassLs[i] = new float[sizeL];
        allpassL[i].setbuffer(bufallpassLs[i], sizeL);
        int sizeR = (allpasstuningL[i] + stereospread) * sampleRate / 44100;
        bufallpassRs[i] = new float[sizeR];
        allpassR[i].setbuffer(bufallpassRs[i], sizeR);
    }

	// Set default values
	allpassL[0].setfeedback(0.5f);
	allpassR[0].setfeedback(0.5f);
	allpassL[1].setfeedback(0.5f);
	allpassR[1].setfeedback(0.5f);
	allpassL[2].setfeedback(0.5f);
	allpassR[2].setfeedback(0.5f);
	allpassL[3].setfeedback(0.5f);
	allpassR[3].setfeedback(0.5f);
	setmode(initialmode);
	
	setroomsize(initialroom);
	setdamp(initialdamp);
	setwidth(initialwidth);
	setdry(initialdry);
	setwet(initialwet);

	// Buffer will be full of rubbish - so we MUST mute them
	mute();
}

revmodel::~revmodel() {
	for(int i = 0; i < numcombs; i++) {
		delete[] bufcombLs[i];
		delete[] bufcombRs[i];
	}

	for(int i = 0; i < numallpasses; i++) {
		delete[] bufallpassLs[i];
		delete[] bufallpassRs[i];
	}

}

void revmodel::mute()
{
	if (getmode() >= freezemode)
		return;

	for (int i=0;i<numcombs;i++)
	{
		combL[i].mute();
		combR[i].mute();
	}
	for (int i=0;i<numallpasses;i++)
	{
		allpassL[i].mute();
		allpassR[i].mute();
	}
}

void revmodel::processreplace(float *inputL, float *inputR, float *outputL, float *outputR, long numsamples, int skip)
{
	float outL,outR,input;

	while(numsamples-- > 0)
	{
		outL = outR = 0;

		if (inputR != nullptr) {
			input = (*inputL + *inputR) * gain;
		} else {
			input = (*inputL) * gain;
		}

		// Accumulate comb filters in parallel
		for(int i=0; i<numcombs; i++)
		{
			outL += combL[i].process(input);
			outR += combR[i].process(input);
		}

		// Feed through allpasses in series
		for(int i=0; i<numallpasses; i++)
		{
			outL = allpassL[i].process(outL);
			outR = allpassR[i].process(outR);
		}

		// Calculate output REPLACING anything already there
		if (inputR != nullptr) {
			*outputL = outL*wet1 + outR*wet2 + *inputL*dry;
			*outputR = outR*wet1 + outL*wet2 + *inputR*dry;
		} else {
			*outputL = outL*wet1 + outR*wet2 + *inputL*dry;
		}

		// Increment sample pointers, allowing for interleave (if any)
		if (inputR != nullptr) {
			inputL += skip;
			inputR += skip;
			outputL += skip;
			outputR += skip;
		} else {
			inputL += skip;
			outputL += skip;
		}
	}
}

void revmodel::processmix(float *inputL, float *inputR, float *outputL, float *outputR, long numsamples, int skip)
{
	float outL,outR,input;

	while(numsamples-- > 0)
	{
		outL = outR = 0;
		input = (*inputL + *inputR) * gain;

		// Accumulate comb filters in parallel
		for(int i=0; i<numcombs; i++)
		{
			outL += combL[i].process(input);
			outR += combR[i].process(input);
		}

		// Feed through allpasses in series
		for(int i=0; i<numallpasses; i++)
		{
			outL = allpassL[i].process(outL);
			outR = allpassR[i].process(outR);
		}

		// Calculate output MIXING with anything already there
		*outputL += outL*wet1 + outR*wet2 + *inputL*dry;
		*outputR += outR*wet1 + outL*wet2 + *inputR*dry;

		// Increment sample pointers, allowing for interleave (if any)
		inputL += skip;
		inputR += skip;
		outputL += skip;
		outputR += skip;
	}
}

void revmodel::update()
{
// Recalculate internal values after parameter change

	int i;

	wet1 = wet*(width/2 + 0.5f);
	wet2 = wet*((1-width)/2);

	if (mode >= freezemode)
	{
		roomsize1 = 1;
		damp1 = 0;
		gain = muted;
	}
	else
	{
		roomsize1 = roomsize;
		damp1 = damp;
		gain = fixedgain;
	}

	for(i=0; i<numcombs; i++)
	{
		combL[i].setfeedback(roomsize1);
		combR[i].setfeedback(roomsize1);
	}

	for(i=0; i<numcombs; i++)
	{
		combL[i].setdamp(damp1);
		combR[i].setdamp(damp1);
	}
}

// The following get/set functions are not inlined, because
// speed is never an issue when calling them, and also
// because as you develop the reverb model, you may
// wish to take dynamic action when they are called.

void revmodel::setroomsize(float value)
{
	roomsize = (value*scaleroom) + offsetroom;
	update();
}

float revmodel::getroomsize()
{
	return (roomsize-offsetroom)/scaleroom;
}

void revmodel::setdamp(float value)
{
	damp = value*scaledamp;
	update();
}

float revmodel::getdamp()
{
	return damp/scaledamp;
}

void revmodel::setwet(float value)
{
	wet = value*scalewet;
	update();
}

float revmodel::getwet()
{
	return wet/scalewet;
}

void revmodel::setdry(float value)
{
	dry = value*scaledry;
}

float revmodel::getdry()
{
	return dry/scaledry;
}

void revmodel::setwidth(float value)
{
	width = value;
	update();
}

float revmodel::getwidth()
{
	return width;
}

void revmodel::setmode(float value)
{
	mode = value;
	update();
}

float revmodel::getmode()
{
	if (mode >= freezemode)
		return 1;
	else
		return 0;
}

//ends
