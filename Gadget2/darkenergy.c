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


void initialize_darkenergy (void) {
		
    printf("Finished initializing dark energy. (DO NOTHING!)\n");
    return;
}

double w(double z)
{
        double ww; // , w0, wa;
	
	// ww=-1.0;
	//	w0=-0.5; wa=-0.5;
	ww=All.w0+(z/(1+z))*All.wa; // example of a redshift-dependent dark energy equation of state parameter w(z).
	return ww;

}

//<AP>

double DarkEnergy (double a)
{
	
	if (All.w0>-1.0001 && All.w0<-0.9999 && All.wa>-0.0001 && All.wa<0.0001) return 1.0 ; 
	return pow(a,-3*(1+All.w0+All.wa))*exp(-3*All.wa*(1-a))
	
}

//</AP>

//#undef NRANSI
