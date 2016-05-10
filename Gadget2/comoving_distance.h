/*
 *  darkenergy.h
 *  Dark Energy Extension for V. Springel's Gadget-2
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */


// Must be included in all .c files that are not main.c
// Contains all references to global variables and functions which are declared in main.c

#ifndef __COMOVING_DISTANCE_H
#define __COMOVING_DISTANCE_H

#include "darkenergy.h"

void derivs_c (double x, double y[], double dydx[]);
double calculate_comoving_distance(double time);

#endif
