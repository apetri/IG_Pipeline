/*
 *  darkenergy.c
 *  Dark Energy Extension for V. Springel's Gadget-2
 *  Also used in Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "main.h"
#include "allocation.h"
#include "darkenergy.h"
#include "darkenergy_support.h"



//#include <nrutil.h>
//#define NRANSI

// Global variables for integration:
double dxsav; // distance at which steps are to be saved.
double *xp; // array of "times" at which output is saved ("time" = the evolution parameter of the differential equation).
double **yp; // array containing the values of the funtions at the output times given in array xp. 
int kmax; // maximal number of intermediate steps saved (last step is always saved)
int kount; // kounts through saved steps up to kmax
int nrhs;   // counts function evaluations (increased by one each time derivs is called, which contains the ODE's) 

double *y2;
double *ap;
double **DEp;

// This function specifies the ordinary differential equations.
// They must by first order ODE's, higher order equations must be rewritten as a system of ordinary first order equations.
// First derivative is on the left-hand side, the corresponding equation on the right-hand side.
// Functions are enumerated from 1 to neqs, where neqs specifies the number of equations and must be set properly within the (sub)routine that calls odeint. 
void derivs (double x, double y[], double dydx[])
{
	nrhs++; //counts function evaluations
	// dydx=f(y(x),x) first order ordinary differential equation:
	dydx[1]=(1+w(x))/(1+x); // arbitrary examples of first order differential equations being solved (this line and the next two).
}

// Initializes w(z) profile; call this function before any first use of calculate_comoving_distance()
void initialize_darkenergy (void) {
	
	int i;
	int neqs; // number of differential equations
	//double ystart[neqs+1];
	double *ystart; // initial conditions array
	double x1, x2; // starting and end point of integration
	double eps, h1, hmin; // Performance control parameters for numerical integrator odeint.
	int nok, nbad; // counts number of good and bad steps (is passed on as a pointer to subroutines called by odeint, so can be modified by those correctly).

	if (parameters.darkenergy_initialized==1) return; // if dark energy is already initialized, don't initialize again, could lead to error during allocation of global variables used for integration. 

	// Number of ordinary differential equations to be solved:
	neqs=1;
	// The Differential Equations are specified in the function derivs.

	// Performance and Output Control Parameters for numerical integrator:
	// Performance:
	eps=pow(10,-18); // Precision, maximal allowed error
	h1=0.01; // guess for first stepsize
	hmin=0; // minimal stepsize (can be zero)
	// Output (output stored in (xp, yp[])):
	kmax=100000; // maximum number of intermediate steps stored (first one and last one are always stored, and count towards the total number of steps specified by kmax).
	dxsav=0.0001; // steps saved only in intervals larger than this value.

	// Allocate arrays for differential equation ("time" parameter if the equation is xp, the functions are enumerated by yp[1-neqs]): 
	xp=Vector(kmax); // Initializes vector, ready for NR-C (Numerical Recipes in C) component enumeration 1 through kmax.
	yp=Matrix(neqs,kmax); // Initializes neqs x kmax matrix with NR-C enumeration. 
	ystart=Vector(neqs); // Initial conditions (position) for functions solving the differential equations (only one in example here).
	// WARNING: NEVER call xp[0], yp[0][...], or ystart[0] !!! Count starts at 1. (Otherwise you will overwrite some other variables!)

	// Allocate dark energy array:
	ap=Vector(kmax);
	DEp=Matrix(neqs,kmax);
		
	//Initial conditions (for first oder equation example here, only one starting value, no derivative, needed):
	ystart[1]=0.0; // function value of first ODE at starting point is 0, because it's an integral.
	x1=0; // starting point of integration is at redshift 0.
	x2=100; // end point of integration at this redshift (redshift needs to be larger than begin of N-body simulation, make higher if necessary).
	
	// Call driver for numerical integrator with above parameters (the driver calls then further subroutines):
	odeint(ystart, neqs, x1, x2, eps, h1, hmin, &nok, &nbad, derivs, rkqs);

	//	printf("Kount: %d.\n", kount);

	// Sample output to check that everything is o.k. and demonstrate how the integrator works:
	// Output should be only correct for writeouts from i=1 to i=kmax. The rest is included just as a reference.
	
	// printf("ystart, nok, nbad, nrhs: %e %d %d %d\n", ystart[1], nok, nbad, nrhs);

	// Before splining, replace redshift z by scale factor a (the name of the variable is xp), and integral by whole dark energy factor expression, reorder by ascending scale factor:
	for (i=1;i<=kount;i++)
	{
		ap[i]=1.0/(1.0+xp[kount+1-i]);
		DEp[1][i]=exp(3*yp[1][kount+1-i]);
		// printf("Eq 1: i, xp, yp: %d --  %e %e\n", i, ap[i], DEp[1][i]);
	}
	//Now can spline this final expression as a function of scale factor:

	// Now do interpolation of above tabulated solutions:
	y2=Vector(kmax);
	// Initialize spline (need only do once):
	// Arguments for spline(): table of arguments of function (x), table of function values at those arguments (y(x)), number of points tabulated, first derivatives at first and last point, output: second derivatives of function at tabulated points).
	//spline(xp, yp[1], kmax, yp[1][1], yp[1][kmax], y2);
	spline(ap, DEp[1], kount, DEp[1][1], DEp[1][kount], y2);
	// Finding the first derivatives at start and end point is easy, because the ODE in derivs is always given as a first order equation, so it's always just the right-hand side of the equation, evaluated at [1] and [kmax] respectively.

	free_Vector(xp);
	free_Matrix(yp, neqs);
	free_Vector(ystart);
	
	// Do not free those until the very end of the whole Instector Gadget run:
	//free_Vector(y2);
	//free_Vector(ap);
	//free_Matrix(DEp, neqs);

	parameters.darkenergy_initialized=1;

    printf("Finished initializing dark energy.\n");
    return;
}

double w(double z)
{
        double ww;
	
	ww=parameters.w0+(z/(1.0+z))*parameters.wa; // example of a redshift-dependent dark energy equation of state parameter w(z).
	return ww;

}

double DarkEnergy (double a)
{
	double yy;
	if (parameters.w0>-1.0001 && parameters.w0<-0.9999 && parameters.wa>-0.0001 && parameters.wa<0.0001) yy=1.0; // in this case the model is so close to LCDM that we assume that's what was intended (avoids wiggling of spline, because for w=-1 the integrator makes only very few steps because the integrand is always zero).
	else splint(ap, DEp[1], y2, kount, a, &yy);
	
	return yy;
}


//#undef NRANSI