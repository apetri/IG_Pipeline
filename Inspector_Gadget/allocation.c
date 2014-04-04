/*
 *  allocation.c
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 4/30/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <time.h>

#include "main.h"
#include "mathematics.h"
#include "allocation.h"



// Subroutines for convenient initialization of arrays (1D vectors and 2D matrices)
// complying with the enumeation convention of numerical recipes, staring to count at 1.

double* Vector(int m) // Initializes 1D array (declare array by "double *arrayname;")
// with m rows. Count starts at 1 and runs through m.
// This function was written to initialize 1D arrays for use with Numerical Recipes in C,
// which use this enumeration of components.
{
         int row;
		 double *vector = (double *) malloc(m * sizeof(double));
         assert(vector != NULL);
		 vector--;
         for (row = 1; row <= m; row++) {
			// Initialize new vector with zeroes:
               vector[row] = 0.0;
         } /* END for */
		return vector;
}

double** Matrix(int m, int n) // Initializes 2D array (declare array by "double **arrayname;")
// with m rows and n columns. Counts start at 1 and run through m and n, respectively.
// This function was written to initialize 2D arrays for use with Numerical Recipes in C,
// which use this enumeration of components.
{
         int row, col;
		 double **matrix = (double **) malloc(m * sizeof(double *));
         assert(matrix != NULL);
		 matrix--;
         for (row = 1; row <= m; row++) {
            matrix[row] = (double *) malloc(n * sizeof(double));
            assert(matrix[row] != NULL);
			matrix[row]--;
			// Initialize new matrix with zeroes:
            for (col = 1; col <= n; col++)
               matrix[row][col] = 0.0;
         } /* END for */
		return matrix;
}


// Memory freeing routines for above allocations of Vector and Matrix:

void free_Vector(double* vector)
{
	vector++; // compensate for offset "vector--;" given at allocation of vector, to conform with Numerical Recipes components enumeration, starting at 1 rather than 0.
	free(vector);
}

void free_Matrix(double** matrix, int m) // need to specify how many rows original matrix had in this routine! (could be improved) 
{
	int row;
	for (row = 1; row <= m; row++)
	{
	matrix[row]++; // compensate for offset "matrix[row]--;" at matrix allocation.
	free(matrix[row]);
	}
	matrix++;
	free(matrix);
}

double* Vector0(int m) // Initializes 1D array (declare array by "double *arrayname;")
// with m rows. Count starts at 0 and runs through m-1, as is standard in C.
{
         int row;
		 double *vector = (double *) malloc(m * sizeof(double));
         assert(vector != NULL);
		 
         for (row = 0; row < m; row++) {
			// Initialize new vector with zeroes:
               vector[row] = 0.0;
         } /* END for */
		return vector;
}

double** Matrix0(int m, int n) // Initializes 2D array (declare array by "double **arrayname;")
// with m rows and n columns. Counts start at 0 and run through m-1 and n-1, respectively, as is standard in C.
{
         int row, col;
		 double **matrix = (double **) malloc(m * sizeof(double *));
         assert(matrix != NULL);
		 
         for (row = 0; row < m; row++) {
            matrix[row] = (double *) malloc(n * sizeof(double));
            assert(matrix[row] != NULL);
			
			// Initialize new matrix with zeroes:
            for (col = 0; col < n; col++)
               matrix[row][col] = 0.0;
         } /* END for */
		 return matrix;
}

double*** Matrix0_3D(int l, int m, int n) // Initializes 3D array (declare array by "double ***arrayname;")
{
         int row;
		 double ***vector = (double ***) malloc(l * sizeof(double **));
         assert(vector != NULL);
		 
		 // Allocate and initialize new matrix with zeroes:
         for (row = 0; row < l; row++) {
			vector[row] = Matrix0(m,n);
         } /* END for */
		return vector;
}

double**** Matrix0_4D(int k, int l, int m, int n) // Initializes 4D array (declare array by "double ****arrayname;")
{
         int row, col;
		 double ****matrix = (double ****) malloc(k * sizeof(double ***));
         assert(matrix != NULL);
		 
         for (row = 0; row < k; row++) {
            matrix[row] = malloc(l * sizeof(double **));
            assert(matrix[row] != NULL);
			
			// Allocate and initialize new matrix with zeroes:
            for (col = 0; col < l; col++)
               matrix[row][col] = Matrix0(m,n);
         } /* END for */
		 return matrix;
}


void free_Vector0(double* vector)
{
	free(vector);
}

void free_Matrix0(double** matrix, int m) // need to specify how many rows original matrix had in this routine! (could be improved) 
{
	int row;
	for (row = 0; row < m; row++)
	{
		free(matrix[row]);
	}
	free(matrix);
}

// Should write freeing functions here for Matrix3D_0 and Matrix4D_0. Are presently only being allocated and automatically freed when program finishes.



/////////////// TIMING OF RUN: /////////////////////
void timing(void)
{
  // clock_t time_now;
  time_t time_now;

	div_t temp1, temp2;
	float time_ratio, uhr1, uhr2, runtime, timediff;
	int seconds, runtime_H, runtime_M, runtime_S;
	
	// time_now=clock(); // get time now.
	time_now=time(NULL); // get time now.
	
	runtime=((long) (time_now-parameters.start_time));
	seconds=((int) runtime);
	
	// Convert to hours:minutes:seconds :
	temp1=div(seconds, 60);
	temp2=div(temp1.quot, 60);
	runtime_S=temp1.rem;
	runtime_M=temp2.rem;
	runtime_H=temp2.quot; 
	
	timediff=runtime-parameters.runtime; // difference to previous call of timing function.
	
	printf("Process %d (superrank) %d (local process number): Time: %d:%d:%d (h:min:sec), %e sec, difference to previous timing step: %e sec.\n", superrank, parameters.process_number, runtime_H, runtime_M, runtime_S, runtime, timediff);
	parameters.runtime=runtime; // make new time now the old time (for next call of this routine).
	fflush(stdout);
}


