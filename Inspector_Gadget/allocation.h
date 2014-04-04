/*
 *  allocation.h
 *  Integrator-Test
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 4/30/07.
 *  Copyright 2007. All rights reserved.
 *
 */

double* Vector(int m);
double** Matrix(int m, int n);
double* Vector0(int m);
double** Matrix0(int m, int n);
double*** Matrix0_3D(int l, int m, int n);
double**** Matrix0_4D(int k, int l, int m, int n);
struct matrix2x2** ArrayMatrix(int m, int n);
void free_Vector(double* vector);
void free_Matrix(double** matrix, int m);
void free_Vector0(double* vector);
void free_Matrix0(double** matrix, int m);
void free_ArrayMatrix(struct matrix2x2 **A, int m);

void timing(void);

