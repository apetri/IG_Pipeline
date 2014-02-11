/*
 *  power_spectrum.h
 *  FITS-Modifier
 *
 *  Created by Jan Kratochvil on 4/4/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

void Power_Spectrum3D_MPI(fftw_complex *data3d, ptrdiff_t M0, ptrdiff_t M1, ptrdiff_t M2, ptrdiff_t local_n0_3d, ptrdiff_t local_0_start_3d, double boxsize, int number_of_bins, char filename[], int ThisTask);
void Power_Spectrum2D(double* imagearray, int nx, int ny, double survey_angle, int number_of_bins, char filename[]);
