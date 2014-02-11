/*
 *  main.h
 *  3D Power Spectrum Calculator (and FFTW3-MPI Test)
 *
 *  Created by Jan Michael Kratochvil on 11/02/2012.
 *  Copyright 2012 Jan Michael Kratochvil. All rights reserved.
 *
 */


// Global variables for determination of endianness of binary Gadget-2 snapshot files:
extern int machine_endian;
extern int snapshot_endian;

extern int superrank, supersize;

// fftw_complex my_function(int i, int j);
void show(fftw_complex *data, ptrdiff_t local_0_start, ptrdiff_t local_n0, ptrdiff_t N1, int ThisTask, int NTasks);

// fftw_complex my_function3d(int i, int j, int k);
void show3d(fftw_complex *data, ptrdiff_t local_0_start_3d, ptrdiff_t local_n0_3d, ptrdiff_t M1, ptrdiff_t M2, int ThisTask, int NTasks);

void project_to_2D(fftw_complex *data3d, int N0_local, int N0_local_start, int N0, int N1, int N2, double *imagearray, int projection_axis, int min_cell, int max_cell, long *naxes);


