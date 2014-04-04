/*
 *  mathematics.h
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 5/23/07.
 *  Copyright 2007. All rights reserved.
 *
 */


struct matrix2x2
{
	double m11;
	double m12;
	double m21;
	double m22;
} ;


struct matrix2x2 M2x2_add(struct matrix2x2 A, struct matrix2x2 B);
struct matrix2x2 M2x2_mult(struct matrix2x2 A, struct matrix2x2 B);
struct matrix2x2 M2x2_scalar(double scale, struct matrix2x2 A);

double ran2(int *idum);

void ensure_random(void);

struct scramble Scrambler(int scramble_mode, double ra1, double ra2, double ra3, double ra4, double ra5, double ra6);


// Bicubic interpolation routines from Numerical Recipes in C:
void bcucof(double y[], double y1[], double y2[], double y12[], double d1, double d2, double **c);
void bcuint(double y[], double y1[], double y2[], double y12[], double x1l, double x1u, double x2l, double x2u, double x1, double x2, double *ansy, double *ansy1, double *ansy2);

