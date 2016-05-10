#ifndef __DARKENERGY_H
#define __DARKENERGY_H

typedef struct {

	double w0;
	double wa;

} DECosmo ;

// Global variables for integration (dark energy density)
double dxsav; // distance at which steps are to be saved.
double *xp; // array of "times" at which output is saved ("time" = the evolution parameter of the differential equation).
double **yp; // array containing the values of the funtions at the output times given in array xp. 
int kmax; // maximal number of intermediate steps saved (last step is always saved)
int kount; // kounts through saved steps up to kmax
int nrhs;   // counts function evaluations (increased by one each time derivs is called, which contains the ODE's) 

double *y2;
double *ap;
double **DEp;


// Global variables for integration (comoving distance)
double dxsav_c; // distance at which steps are to be saved.
double *xp_c; // array of "times" at which output is saved ("time" = the evolution parameter of the differential equation).
double **yp_c; // array containing the values of the funtions at the output times given in array xp. 
int kmax_c; // maximal number of intermediate steps saved (last step is always saved)
int kount_c; // kounts through saved steps up to kmax
int nrhs_c;   // counts function evaluations (increased by one each time derivs is called, which contains the ODE's) 

// double *y2_c;
double *ap_c;
double **DEp_c;

double DarkEnergy (double a, DECosmo *cosmo);

/*
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */

//////////////////////////////////
// INTEGRATION:
//////////////////////////////////

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])));

void odeint_c(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])));
	
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
	
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double []));


////////////////////////
// INTERPOLATION:
////////////////////////

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);


////////////////////////////////
// ALLOCATION:
////////////////////////////////

double* Vector(int m);
double** Matrix(int m, int n);
void free_Vector(double* vector);
void free_Matrix(double** matrix, int m);

#endif
