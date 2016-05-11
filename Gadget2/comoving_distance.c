/*
 *  darkenergy.c
 *  Dark Energy Extension for V. Springel's Gadget-2
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "allvars.h"

#include "darkenergy.h"
#include "comoving_distance.h"

//#include <nrutil.h>
//#define NRANSI

 
// Global variables for integration (comoving distance)
extern double dxsav_c; // distance at which steps are to be saved.
extern double *xp_c; // array of "times" at which output is saved ("time" = the evolution parameter of the differential equation).
extern double **yp_c; // array containing the values of the funtions at the output times given in array xp. 
extern int kmax_c; // maximal number of intermediate steps saved (last step is always saved)
extern int kount_c; // kounts through saved steps up to kmax
extern int nrhs_c;   // counts function evaluations (increased by one each time derivs is called, which contains the ODE's) 

// double *y2_c;
extern double *ap_c;
extern double **DEp_c;


// This function specifies the ordinary differential equations.
// They must by first order ODE's, higher order equations must be rewritten as a system of ordinary first order equations.
// First derivative is on the left-hand side, the corresponding equation on the right-hand side.
// Functions are enumerated from 1 to neqs, where neqs specifies the number of equations and must be set properly within the (sub)routine that calls odeint. 
void derivs_c (double x, double y[], double dydx[])
{
	nrhs_c++; //counts function evaluations
	// dydx=f(y(x),x) first order ordinary differential equation:
	dydx[1]=1.0/(sqrt(All.Omega0*(1.0+x)*(1.0+x)*(1.0+x)+All.OmegaLambda*DarkEnergy(1.0/(1.0+x),&de_cosmo))); // arbitrary examples of first order differential equations being solved (this line and the next two).
}


double calculate_comoving_distance(double a) {
	
	int i;
	int neqs; // number of differential equations
	//double ystart[neqs+1];
	double *ystart; // initial conditions array
	double x1, x2; // starting and end point of integration
	double eps, h1, hmin; // Performance control parameters for numerical integrator odeint.
	int nok, nbad; // counts number of good and bad steps (is passed on as a pointer to subroutines called by odeint, so can be modified by those correctly).

	double speedoflight=2.99792458e5; // in km/s (exact value, meter is defined that way).
	double comoving_distance;
	comoving_distance=0.0;

	// Number of ordinary differential equations to be solved:
	neqs=1;
	// The Differential Equations are specified in the function derivs.

	// Performance and Output Control Parameters for numerical integrator:
	// Performance:
	eps=pow(10,-18); // Precision, maximal allowed error
	h1=0.01; // guess for first stepsize
	hmin=0; // minimal stepsize (can be zero)
	// Output (output stored in (xp, yp[])):
	kmax_c=1; //100000; // maximum number of intermediate steps stored (first one and last one are always stored, and count towards the total number of steps specified by kmax).
	dxsav_c=0.0001; // steps saved only in intervals larger than this value.

	// Allocate arrays for differential equation ("time" parameter if the equation is xp, the functions are enumerated by yp[1-neqs]): 
	xp_c=Vector(kmax_c); // Initializes vector, ready for NR-C (Numerical Recipes in C) component enumeration 1 through kmax.
	yp_c=Matrix(neqs,kmax_c); // Initializes neqs x kmax matrix with NR-C enumeration. 
	ystart=Vector(neqs); // Initial conditions (position) for functions solving the differential equations (only one in example here).
	// WARNING: NEVER call xp[0], yp[0][...], or ystart[0] !!! Count starts at 1. (Otherwise you will overwrite some other variables!)

	// Allocate dark energy array:
	ap_c=Vector(kmax_c);
	DEp_c=Matrix(neqs,kmax_c);
		
	//Initial conditions (for first oder equation example here, only one starting value, no derivative, needed):
	ystart[1]=0.0; // function value of first ODE at starting point is 0, because it's an integral.
	x1=0.0; // starting point of integration is at redshift 0.
	x2=(1.0/a)-1.0; // end point of integration at this redshift (redshift needs to be larger than begin of N-body simulation, make higher if necessary).
	
	// Call driver for numerical integrator with above parameters (the driver calls then further subroutines):
	odeint_c(ystart, neqs, x1, x2, eps, h1, hmin, &nok, &nbad, derivs_c, rkqs);

	// printf("Kount: %d.\n", kount_c);

	// Sample output to check that everything is o.k. and demonstrate how the integrator works:
	// Output should be only correct for writeouts from i=1 to i=kmax. The rest is included just as a reference.
	
	// printf("ystart, nok, nbad, nrhs: %e %d %d %d\n", ystart[1], nok, nbad, nrhs_c);

	// Before splining, replace redshift z by scale factor a (the name of the variable is xp), and integral by whole dark energy factor expression, reorder by ascending scale factor:
	for (i=1;i<=kount_c;i++)
	{
		ap_c[i]=1.0/(1.0+xp_c[kount_c+1-i]);
		DEp_c[1][i]=yp_c[1][kount_c+1-i];
		//		printf("Eq 1: i, xp, yp: %d --  %e %e\n", i, ap_c[i], DEp_c[1][i]);
	}

	// printf("Pre-Comoving distance: %e, scale factor %e.\n", DEp_c[1][kount_c], ap_c[kount_c]);
	comoving_distance=speedoflight*DEp_c[1][kount_c]/100.0*1000.0; // gives result in [h^-1 kpc].
	// Comoving distance = c * \int_{0}^{z} dz' 1/H(z'), where H(z) = H0 * sqrt(OM*(1+z)^3+OL*DarkEnergy(1/(z+1))).
	// Then the get h^-1 in, write the H0 in the denominator as 100 * h, multiply everything by h (such that one needs to divide by h to get kpc), which removes the h from the denominator to make the quantity h-free.  
	// The factor of 1000.0 converts to kpc (since H0=100*h is in km/s/Mpc).

	free_Vector(xp_c);
	free_Matrix(yp_c, neqs);
	free_Vector(ystart);
	
	// Do not free those until the very end of the whole Gadget run:
	//free_Vector(y2_c);
	free_Vector(ap_c);
	free_Matrix(DEp_c, neqs);
		
	//    printf("Finished calculating comoving distance, and dark energy parameters are: w0=%e, wa=%e.\n", All.w0, All.wa);
    return comoving_distance;
}

//#undef NRANSI
