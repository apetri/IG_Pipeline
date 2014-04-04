/*
 *  darkenergy_support.h
 *  Comoving Distance Computation for Inspector Gadget and for V. Springel's Gadget-2
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */


// will be included at the beginning of main.c
// list (empty) declaration here of all subroutines that are defined in integration.h, so they appear at the beginning of the main.c file, due to the include.

//////////////////////////////////
// INTEGRATION:
//////////////////////////////////

void odeint_c(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])));



