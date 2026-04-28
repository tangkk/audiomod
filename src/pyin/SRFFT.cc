#include "SRFFT.h"

SRFFT::SRFFT(int fft_len)  
{
	int i, j, k;

	n = fft_len;

	cosTab = new int[n];
	sinTab = new int[n];
	temp   = new COMPLEX[n];

	for(i=0; i<(n); i++)
	{
		cosTab[i]=cos( 2 * M_PI * i * 1.0/n) * 1073741824;  // 2^30
		sinTab[i]=sin( 2 * M_PI * i * 1.0/n) * 1073741824;  // 2^30
	}
    
	reorder = new int[n];

    reorder[0]=0;
    for(j=0,i=0;i<(n-1); i++)
    {   
        k=n/2;
        while(k<(j+1))
        {  
           j=j-k;
           k=k/2;
        }
        j=j+k;
        
        reorder[i+1]=j;
    }

}

SRFFT::~SRFFT()
{
	delete [] cosTab;
	delete [] sinTab;
	delete [] reorder;
	delete [] temp;

}

void SRFFT::Split_radix(COMPLEX *x)
{
    int e, a, a3, i,j,k,m,i1,i2,i3,n1,n2,n4,id,is;
    int r1,r2;
    int s1,s2,s3;
    int cc1,cc3,ss1,ss3;

    for (j=1,i=1;i<n;i++)
    {
        m=i;                 
        j=2*j;
        if(j==n) break;
    }

    n2=2*n;                  

    for(k=1;k<m;k++)         
    {
        n2=n2/2;             
        n4=n2/4;             
        e=n/n2;              
        a=0;                 

        for(j=0;j<n4;j++)    
        {
            a3=3*a;          
            cc1=cosTab[a];   
            ss1=sinTab[a];   
            cc3=cosTab[a3];  
            ss3=sinTab[a3];  
            a=(j+1)*e;       
            is=j;            
            id=2*n2;         
            do{
                  for ( i=is;i<(n-1);i=i+id)  
                  {
                      i1=i+n4;             
                      i2=i1+n4;            
                      i3=i2+n4;            
                      
                      r1=x[i].real-x[i2].real;
                      x[i].real=x[i].real+x[i2].real;
					  r2=x[i1].real-x[i3].real;
					  x[i1].real=x[i1].real+x[i3].real;   
                      
					  s1=x[i].imag-x[i2].imag;       
					  x[i].imag=x[i].imag+x[i2].imag;     
					  s2=x[i1].imag-x[i3].imag;      
					  x[i1].imag=x[i1].imag+x[i3].imag;  
                      
                      s3=r1-s2;            
                      r1=r1+s2;            
                      s2=r2-s1;            
                      r2=r2+s1;            
                      
					  x[i2].real = (((long long)( r1)*cc1)>>30) - (((long long)s2*ss1)>>30); 
					  x[i2].imag = (((long long)(-s2)*cc1)>>30) - (((long long)r1*ss1)>>30); 
					  x[i3].real = (((long long)( s3)*cc3)>>30) + (((long long)r2*ss3)>>30); 
					  x[i3].imag = (((long long)( r2)*cc3)>>30) - (((long long)s3*ss3)>>30); 
                  }
                  is=2*id-n2+j;   
                  id=4*id;        
              }   while(is<(n-1));
            }
        }

    is=0;  
    id=4;  
    do
    {
        for(i=is;i<n;i=i+id)
        {   
            i1=i+1;
			r1=x[i].real;
			r2=x[i].imag;
			x[i].real=r1+x[i1].real;
			x[i].imag=r2+x[i1].imag;
			x[i1].real=r1-x[i1].real;
			x[i1].imag=r2-x[i1].imag;
        }
        is=2*id-2;
        id=4*id;
    }while(is<(n-1));
    
    n1=n-1;
    
    for(j=0,i=0;i<n1;i++)
    {   
        j=reorder[i];
        if(i<j)
        {   
			r1=x[j].real;
			s1=x[j].imag;
			x[j].real=x[i].real;
			x[j].imag=x[i].imag;
			x[i].real=r1;
			x[i].imag=s1;
        }
    }

}

void SRFFT::invert_FFT(COMPLEX *x)
{
    int i;
	int shift_bit;

	switch(n)
	{
	case 2048: shift_bit = 11; break;
	case 1024: shift_bit = 10; break;
	case  512: shift_bit =  9; break;
	case  256: shift_bit =  8; break;
	case  128: shift_bit =  7; break;
	case   64: shift_bit =  6; break;
	case   32: shift_bit =  5; break;
	case   16: shift_bit =  4; break;
	case    8: shift_bit =  3; break;
	case    4: shift_bit =  2; break;
	default  : shift_bit =  1;
	}



    for(i=0; i<n; i++)
    {
		x[i].imag *= (-1);
    }
    
	Split_radix(x);
    
    for(i=0; i<n; i++)
    {
		x[i].imag *= (-1);

		x[i].real >>= shift_bit;
		x[i].imag >>= shift_bit;
		
    }

}



void SRFFT::Split_radix(COMPLEX *x, COMPLEX *xx)
{

	int m = n>>1;
	for(int i=0;i<n;i++ )
	{
		temp[i].real = x[i].real;
		temp[i].imag = xx[i].real;
	}
	Split_radix(temp);

	x[0].real  = temp[0].real;
	x[0].imag  = 0;
	xx[0].real = temp[0].imag;
	xx[0].imag = 0;
	x[m].real  = temp[m].real;
	x[m].imag  = 0;
	xx[m].real = temp[m].imag;
	xx[m].imag = 0;

	for(int i=1;i<m;i++)
	{
		x[i].real = (temp[i].real + temp[n-i].real)>>1;
		x[i].imag = (temp[i].imag - temp[n-i].imag)>>1;
		xx[i].real =(temp[i].imag + temp[n-i].imag)>>1;
		xx[i].imag =(temp[n-i].real-temp[i].real)>>1;

		x[n-i].real = x[i].real;
		x[n-i].imag = -x[i].imag;
		xx[n-i].real =xx[i].real;
		xx[n-i].imag =-xx[i].imag;
	}

}
	
void SRFFT::invert_FFT(COMPLEX *x, COMPLEX *xx)
{
	int i;
	for(i=0; i<n; i++)
	{
		temp[i].real  = x[i].real - xx[i].imag;
		temp[i].imag  = x[i].imag + xx[i].real;
	}

	invert_FFT(temp);

	for(int i=0;i<n;i++)
	{
		x[i].real = temp[i].real;
	//	x[i].imag = 0;
		xx[i].real =temp[i].imag;
	//	xx[i].imag =0;
	}

}