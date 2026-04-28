#include <math.h>

#ifndef SRFFT_H
#define SRFFT_H

#ifndef  M_PI
#define M_PI  (3.14159265358979)
#endif

typedef struct
{
	int real;
	int imag;
} COMPLEX;

typedef struct
{
	float real;
	float imag;
} fCOMPLEX;

class SRFFT
{
public :
	SRFFT(int fft_len);
	
	~SRFFT();
	
	void Split_radix(COMPLEX *x);
	
	void invert_FFT(COMPLEX *x);
	
	void Split_radix(COMPLEX *x, COMPLEX *xx);
	
	void invert_FFT(COMPLEX *x, COMPLEX *xx);
	
private :
	int *cosTab;
	int *sinTab;
	int *reorder;
	int n;
	COMPLEX *temp;
};

#endif
