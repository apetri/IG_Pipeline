/*
 *  mathematics.c
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 5/23/07.
 *  Copyright 2007. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <time.h>

#include <complex.h>
#include <fftw3.h>

#include "mathematics.h"
#include "allocation.h"
#include "main.h"

struct matrix2x2 M2x2_add(struct matrix2x2 A, struct matrix2x2 B)
{
	struct matrix2x2 C;
	
	C.m11=A.m11+B.m11;
	C.m12=A.m12+B.m12;
	C.m21=A.m21+B.m21;
	C.m22=A.m22+B.m22;

	return C;
}


struct matrix2x2 M2x2_mult(struct matrix2x2 A, struct matrix2x2 B)
{
	struct matrix2x2 C;
	
	C.m11=A.m11*B.m11+A.m12*B.m21;
	C.m12=A.m11*B.m12+A.m12*B.m22;
	C.m21=A.m21*B.m11+A.m22*B.m21;
	C.m22=A.m21*B.m12+A.m22*B.m22;

	return C;
}


struct matrix2x2 M2x2_scalar(double scalar, struct matrix2x2 A) // multiplication with a scalar:
{
	A.m11 *= scalar;
    A.m12 *= scalar;
	A.m21 *= scalar;
	A.m22 *= scalar;
	
	return A;
}



#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(int *idum) // Random number generator from Numerical Recipes in C
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


void ensure_random(void)
{
  int i, presequence;
  double number;
  presequence=11731+187*parameters.process_number;
  for (i=0; i<presequence; i++)
  {
      number=ran2(&parameters.seed);
  }

}



struct scramble Scrambler(int scramble_mode, double ra1, double ra2, double ra3, double ra4, double ra5, double ra6)
{
	struct scramble scrambled;
	
	if (scramble_mode==0) // then no scrambling/turning of coordinates
	{
		scrambled.x1=0;
		scrambled.x2=1;
		scrambled.x3=2;
		scrambled.signx1=1;
		scrambled.signx2=1;
		scrambled.signx3=1;
		
		scrambled.centeroffset1=0.5*parameters.boxsize;
		scrambled.centeroffset2=0.5*parameters.boxsize;
		scrambled.centeroffset3=0.5*parameters.boxsize;
	}
	
	
	if (scramble_mode==2)
	{
		// Fully randomized scramble mode for rotating and mirroring. Crosshairs not implemented, this mode works only with the advanced FITS file plane generating mode.
		
		double random_number1, random_number2, random_number3, random_number4, random_number5, random_number6;

                random_number1=ra1;
                random_number2=ra2;
                random_number3=ra3;
                random_number4=ra4;
                random_number5=ra5;
                random_number6=ra6;

		if (feedback >3) printf("Random Numbers in fully random mode of Scrambler: %e %e %e %e %e %e \n", random_number1, random_number2, random_number3, random_number4, random_number5, random_number6);
		
		if (random_number1<0.333333333333333333333)
		{	
			scrambled.x1=0;
			if (random_number2<0.5)
			{ 
				scrambled.x2=1;
				scrambled.x3=2;
			}
			else
			{ 
				scrambled.x2=2;
				scrambled.x3=1;
			}
		}
		else if (random_number1>0.666666666666666666666)
		{
			scrambled.x1=2;
			if (random_number2<0.5)
			{ 
				scrambled.x2=0;
				scrambled.x3=1;
			}
			else
			{ 
				scrambled.x2=1;
				scrambled.x3=0;
			}
		}
		else 
		{
			scrambled.x1=1;
			if (random_number2<0.5)
			{ 
				scrambled.x2=0;
				scrambled.x3=2;
			}
			else
			{ 
				scrambled.x2=2;
				scrambled.x3=0;
			}
		}
		// Axes Mirroring:
		if (random_number4<0.5) scrambled.signx1=1;
		else scrambled.signx1=-1;
		if (random_number5<0.5) scrambled.signx2=1;
		else scrambled.signx2=-1;
		if (random_number6<0.5) scrambled.signx3=1;
		else scrambled.signx3=-1;
		
		if (feedback>3) printf("Randomized Coordinates (0, 1, 2, three signs): %d %d %d %d %d %d \n", scrambled.x1, scrambled.x2, scrambled.x3, scrambled.signx1, scrambled.signx2, scrambled.signx3);
				
	}
	
	
	return scrambled;
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// Bicubic interpolation routines from Numerical Recipes in C:

void bcucof(double y[], double y1[], double y2[], double y12[], double d1, double d2,
	    double **c)
{
  static int wt[16][16]=
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
      -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
      2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
      0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
      0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
      0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
      -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
      9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
      -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
      2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
      -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
      4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1};
  int l,k,j,i;
  double xx,d1d2,cl[16],x[16];
  
  d1d2=d1*d2;
  for (i=1;i<=4;i++) {
    x[i-1]=y[i];
    x[i+3]=y1[i]*d1;
    x[i+7]=y2[i]*d2;
    x[i+11]=y12[i]*d1d2;
  }
  for (i=0;i<=15;i++) {
    xx=0.0;
    for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
    cl[i]=xx;
  }
  l=0;
  for (i=1;i<=4;i++)
    for (j=1;j<=4;j++) c[i][j]=cl[l++];
}


//#define NRANSI
// #include "nrutil.h"

void bcuint(double y[], double y1[], double y2[], double y12[], double x1l, double x1u, double x2l, double x2u, double x1, double x2, double *ansy, double *ansy1, double *ansy2)
{
  void bcucof(double y[], double y1[], double y2[], double y12[], double d1, double d2, double **c);
  int i;
  double t,u,d1,d2,**c;

  
  // c=matrix(1,4,1,4);
  c=Matrix(4,4);
  d1=x1u-x1l;
  d2=x2u-x2l;
  bcucof(y,y1,y2,y12,d1,d2,c);
  if (x1u == x1l || x2u == x2l) 
    {
      printf("Bad input in routine bcuint: %e %e, %e %e\n", x1l, x1u, x2l, x2u);
      fflush(stdout);
      exit(1);
    }
  t=(x1-x1l)/d1;
  u=(x2-x2l)/d2;
  // *ansy=(*ansy2)=(*ansy1)=0.0;
  // printf("ANSI: %e\n", *ansy);
  // fflush(stdout);

  
  *ansy=0.0;
  *ansy1=0.0;
  *ansy2=0.0;
  
  for (i=4;i>=1;i--) {
    
    *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
    *ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
    *ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
    
  }
  
  *ansy1 /= d1;
  *ansy2 /= d2;
  
  // free_matrix(c,1,4,1,4);
  free_Matrix(c, 4);
  
}
//#undef NRANSI

