/*
 *  darkenergy.h
 *  Dark Energy Extension for V. Springel's Gadget-2
 *  Also used in Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */


// Must be included in all .c files that are not main.c
// Contains all references to global variables and functions which are declared in main.c

extern int kmax, kount, nrhs;
extern double dxsav, *xp, **yp;
extern double *y2;
extern double *ap;
extern double **DEp;

void derivs (double x, double y[], double dydx[]);
void initialize_darkenergy (void);

double w(double z);
double DarkEnergy (double a);
