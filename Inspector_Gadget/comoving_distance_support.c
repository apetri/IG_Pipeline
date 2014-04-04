/*
 *  comoving_distance_support.c
 *  Comoving Distance Computation for Inspector Gadget and for V. Springel's Gadget-2
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "allocation.h"
#include "comoving_distance.h"
#include "comoving_distance_support.h"


///////////////////////////////////////////
// INTEGRATION:
///////////////////////////////////////////


//#define NRANSI
//#include "nrutil.h"
#define MAXSTP 1000000
// orignial setting from Numerical Recipes in C was MAXSTP 10000
#define TINY 1.0e-30

// Macros for definitions of simple functions that are called by the integration subroutines from Numerical Recipes in C.
// They have been here reverse engineered based on Numerical Recipes in FORTRAN.
#define FMAX(A, B) (((A)>(B)) ? (A) : (B))
#define FMIN(A, B) (((A)>(B)) ? (B) : (A))
#define SIGN(A, B) ((((B)>=0.0) ? 1.0 : -1.0) * fabs(A))



void odeint_c(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])))
{
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;

	yscal=Vector(nvar);
	y=Vector(nvar);
	dydx=Vector(nvar);

	//double yscal[nvar+1], y[nvar+1], dydx[nvar+1];
	//yscal=vector(1,nvar);
	//y=vector(1,nvar);
	//dydx=vector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);  // returns magnitude of h1 with the sign of x2-x1
	*nok = (*nbad) = kount_c = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax_c > 0) xsav=x-dxsav_c*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,y,dydx);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax_c > 0 && kount_c < kmax_c-1 && fabs(x-xsav) > fabs(dxsav_c)) {
			kount_c++;
			xp_c[kount_c]=x;
			for (i=1;i<=nvar;i++) yp_c[i][kount_c]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		if (hdid == h) (*nok)++; else (*nbad)++;
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax_c) {
				kount_c++;
				xp_c[kount_c]=x;
				for (i=1;i<=nvar;i++) yp_c[i][kount_c]=y[i];
			}
		//	free_vector(dydx,1,nvar);
		//	free_vector(y,1,nvar);
		//	free_vector(yscal,1,nvar);
		free_Vector(dydx);
		free_Vector(y);
		free_Vector(yscal);

			return;
		}
		if (fabs(hnext) <= hmin)
		{
			printf("Step size too small in odeint_c");
			exit(1);
		}
		h=hnext;
	}
	printf("Too many steps in routine odeint_c");
	exit(1);
}
#undef MAXSTP
#undef TINY
//#undef NRANSI





#undef FMAX
#undef FMIN
#undef SIGN

