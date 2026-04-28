/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp

    An API for audio analysis and feature extraction plugins.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright 2006-2012 Chris Cannam and QMUL.
  
    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
    ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Except as contained in this notice, the names of the Centre for
    Digital Music; Queen Mary, University of London; and Chris Cannam
    shall not be used in advertising or otherwise to promote the sale,
    use or other dealings in this Software without prior written
    authorization.
*/

//#include <malloc.h>
#include "FFT.h"
#include <math.h>

namespace Vamp {

static void
	fft(unsigned int n, bool inverse,
	const float *ri, const float *ii,
	float *ro, float *io)
{
	if (!ri || !ro || !io) return;

	unsigned int bits;
	unsigned int i, j, k, m;
	unsigned int blockSize, blockEnd;

	float tr, ti;

	if (n < 2) return;
	if (n & (n-1)) return;

	float angle = 2.0 * M_PI;
	if (inverse) angle = -angle;

	for (i = 0; ; ++i) {
		if (n & (1 << i)) {
			bits = i;
			break;
		}
	}

	int *table = (int *)malloc(n * sizeof(int));

	for (i = 0; i < n; ++i) {
		m = i;
		for (j = k = 0; j < bits; ++j) {
			k = (k << 1) | (m & 1);
			m >>= 1;
		}
		table[i] = k;
	}

	if (ii) {
		for (i = 0; i < n; ++i) {
			ro[table[i]] = ri[i];
			io[table[i]] = ii[i];
		}
	} else {
		for (i = 0; i < n; ++i) {
			ro[table[i]] = ri[i];
			io[table[i]] = 0.0;
		}
	}

	blockEnd = 1;

	for (blockSize = 2; blockSize <= n; blockSize <<= 1) {

		float delta = angle / (float)blockSize;
		float sm2 = -sin(-2 * delta);
		float sm1 = -sin(-delta);
		float cm2 = cos(-2 * delta);
		float cm1 = cos(-delta);
		float w = 2 * cm1;
		float ar[3], ai[3];

		for (i = 0; i < n; i += blockSize) {

			ar[2] = cm2;
			ar[1] = cm1;

			ai[2] = sm2;
			ai[1] = sm1;

			for (j = i, m = 0; m < blockEnd; j++, m++) {

				ar[0] = w * ar[1] - ar[2];
				ar[2] = ar[1];
				ar[1] = ar[0];

				ai[0] = w * ai[1] - ai[2];
				ai[2] = ai[1];
				ai[1] = ai[0];

				k = j + blockEnd;
				tr = ar[0] * ro[k] - ai[0] * io[k];
				ti = ar[0] * io[k] + ai[0] * ro[k];

				ro[k] = ro[j] - tr;
				io[k] = io[j] - ti;

				ro[j] += tr;
				io[j] += ti;
			}
		}

		blockEnd = blockSize;
	}

	if (inverse) {

		float denom = (float)n;

		for (i = 0; i < n; i++) {
			ro[i] /= denom;
			io[i] /= denom;
		}
	}

	free(table);
}

void
FFT::forward(unsigned int n,
	     const float *ri, const float *ii,
	     float *ro, float *io)
{
    fft(n, false, ri, ii, ro, io);
}

void
FFT::inverse(unsigned int n,
	     const float *ri, const float *ii,
	     float *ro, float *io)
{
    fft(n, true, ri, ii, ro, io);
}

}


