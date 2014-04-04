/*
 *  darkenergy.h
 *  Comoving Distance Computation for Inspector Gadget and for V. Springel's Gadget-2
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */


// Must be included in all .c files that are not main.c
// Contains all references to global variables and functions which are declared in main.c

extern int kmax_c, kount_c, nrhs_c;
extern double dxsav_c, *xp_c, **yp_c;
// extern double *y2_c;
extern double *ap_c;
extern double **DEp_c;

void derivs_c (double x, double y[], double dydx[]);
double calculate_comoving_distance(double time);
