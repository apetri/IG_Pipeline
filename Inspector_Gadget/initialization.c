/*
 *  initialization.c
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at the University of Miami on 8/30/11 (originally for a Limber approximation code).
 *  Copyright 2011. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "main.h"
#include "allocation.h"

#include "darkenergy_support.h" // necessary for spline() and splint() from Numerical Recipes in C, which are defined there.
#include "comoving_distance.h"

#include "initialization.h"






// Global variables for integration:
//double dxsav_chi; // distance at which steps are to be saved.
double *xp_chi; // array of "times" at which output is saved ("time" = the evolution parameter of the differential equation).
double *yp_chi; // array containing the values of the funtions at the output times given in array xp. 
int kmax_chi; // maximal number of intermediate steps saved (last step is always saved)
//int kount_chi; // kounts through saved steps up to kmax
//int nrhs_chi;   // counts function evaluations (increased by one each time derivs is called, which contains the ODE's) 
double *y2_chi; // contains derivatives for spline


// Initialize spline for chi(z) table:
// (spline will return comoving distance values in units of Mph/h.
void initialize_chi_in_Mpc(void)
{
	int i;
	double zmax, a1, amax, da;
	
    if (parameters.chi_initialized==1) return;
    
    
	// Settable parameters:
	kmax_chi=10000; // kmax_chi=10000 for good precision even at z=0.01 (and reasonable at z=0.001) if zmax is set to 100 (computation takes ca. 1 minute, this is good for production runs, but tedious for debugging). kmax_chi=1000 for reasonable precision at z=0.01 and fast computation (for debugging).
	zmax=100.0; // 100.0 for weak lensing N-body, 1500.0 for CMB lensing (not done by this code); // redshift out to which we want to spline.
	
	a1=1.0; // scale factor today (leave =1).
	amax=a1/(1.0+zmax); // scale factor at starting redshift (or as far as we want to spline here
	
	xp_chi=Vector(kmax_chi); // Arrays in NR in C start at 1.
	yp_chi=Vector(kmax_chi);
	
	// Steps logarithmically spaced in a (gives approximately equal spacing in comoving distance):
	da=fabs(log(a1)-log(amax))/((double) kmax_chi);
	
	// Compute chi(z) table:
	for (i=1; i<=kmax_chi; i++)
	{
		// xp_chi[i]=a1*exp(-(i-1.0)*da);
		// yp_chi[i]=calculate_comoving_distance((1.0/(1.0+xp_chi[i])));
		
		xp_chi[i]=1.0/(a1*exp(-(i-1.0)*da))-1.0;
		yp_chi[i]=calculate_comoving_distance((1.0/(1.0+xp_chi[i])))/1000.0; // divides by 1000 here, because calculate comoving distance in Inspector Gadget returns comoving distance in kpc/h.
		
		
		// printf("i, xp_chi, yp_chi: %d %e %e\n", i, xp_chi[i], yp_chi[i]);
	}
	
	
	y2_chi=Vector(kmax_chi);
	// Initialize spline (need only do once):
	// Arguments for spline(): table of arguments of function (x), table of function values at those arguments (y(x)), number of points tabulated, first derivatives at first and last point, output: second derivatives of function at tabulated points).
	spline(xp_chi, yp_chi, kmax_chi, yp_chi[1], yp_chi[kmax_chi], y2_chi);
	// Finding the first derivatives at start and end point is easy, because the ODE in derivs is always given as a first order equation, so it's always just the right-hand side of the equation, evaluated at [1] and [kmax] respectively.
	
    printf("Warning: chi (comoving distance in Mpc/h) -- to be obtained with function get_chi_in_Mpc(z) -- has been initialized only out to redshift z_max=%e. Do not call get_chi_in_Mpc(z) for redshifts z>z_max. (The function calculate_comoving_distance(scale factor) is still available for aribtrary redshifts; not that that function returns the comoving distance in kpc/h (not Mpc/h).\n", zmax);
    
    parameters.chi_initialized=1;
    
}


// This funtion evaluates spline from chi(z) table generated earlier by repetitive calls of function calculate_comoving_distance:
// (For speed reasons, interpolation is used rather than calculating the comoving distance accurately from scratch for each point needed.)
double get_chi_in_Mpc(double z) // returns comoving distance in Mpc/h (because initialized by initialize_chi_in_Mpc.
{
	double chi;
	// double xx, yy; // points at which we evaluate spline and result we get.
	// xx=1.0; // example value at which to evaluate.
	
	// Evaluate spline:  
	// Arguments for splint(): table of arguments of function (x), table of function values at those arguments (y(x)), output table from function spline above (second derivatives probably), number of tabulated points, point at which splined function is to be evaluated (x), output needs to be supplied as address and is the value of the splined function at that point (y(x)).
	splint(xp_chi, yp_chi, y2_chi, kmax_chi, z, &chi);
	// splint(xp_chi, yp_chi[1], y2_chi, kmax_chi, xx, &yy);
	
	// printf("Spline evaluated at x=%e is y=%e \n", z, chi);
	
	return chi;
	
}

