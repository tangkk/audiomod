// Reverb model declaration
//
// Written by Jezar at Dreampoint, June 2000
// http://www.dreampoint.co.uk
// This code is public domain

#pragma once

#include <iostream>
#include <memory>
#include "../common/filters/comb.h"
#include "../common/filters/allpass.h"
#include "../common/filters/tuning.h"

class revmodel
{
public:
					revmodel(float sampleRate);
					~revmodel();
			void	mute();
			// non-interleave: skip = 1; interleave skip = 2
			void	processmix(float *inputL, float *inputR, float *outputL, float *outputR, long numsamples, int skip = 1);
			// non-interleave: skip = 1; interleave skip = 2
			void	processreplace(float *inputL, float *inputR, float *outputL, float *outputR, long numsamples, int skip = 1);
			void	setroomsize(float value);
			float	getroomsize();
			void	setdamp(float value);
			float	getdamp();
			void	setwet(float value);
			float	getwet();
			void	setdry(float value);
			float	getdry();
			void	setwidth(float value);
			float	getwidth();
			void	setmode(float value);
			float	getmode();
private:
			void	update();
private:
	float	gain;
	float	roomsize,roomsize1;
	float	damp,damp1;
	float	wet,wet1,wet2;
	float	dry;
	float	width;
	float	mode;

	// The following are all declared inline 
	// to remove the need for dynamic allocation
	// with its subsequent error-checking messiness

	// Comb filters
	comb	combL[numcombs];
	comb	combR[numcombs];

	// Allpass filters
	allpass	allpassL[numallpasses];
	allpass	allpassR[numallpasses];

	// Buffers for the combs
    float *bufcombLs[numcombs];
	float *bufcombRs[numcombs];

	// Buffers for the allpasses
    float *bufallpassLs[numallpasses];
	float *bufallpassRs[numallpasses];
};

//ends
