/*
 *  power_spectrum.c
 *  Power Spectrum
 *
 *  Created by Jan Kratochvil on 3/16/08 at Columbia University.
 *  Copyright 2008 Jan Kratochvil. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <complex.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

#include "main.h"

// NOTES:
// External variables needed:
// parameters.nx and parameters.ny: resolution of 2D input map (in pixels). Note that this map is input formally as a 1D-array.
// parameters.nxny=parameters.nx*parameters.ny
// parameters.survey_angle: angular size (in one dimension) of map, in degrees.



//void Power_Spectrum3D(double* imagearray, int nx, int ny, double survey_angle, int number_of_bins, char filename[]) // reads in image array, copies it into Fourier transformable array, Fourier transforms, and then overwrites original image array
void Power_Spectrum3D_MPI(fftw_complex *data3d, ptrdiff_t M0, ptrdiff_t M1, ptrdiff_t M2, ptrdiff_t local_n0_3d, ptrdiff_t local_0_start_3d, double boxsize, int number_of_bins, char filename[], int ThisTask) 
{
	// NOTE: ThisTask as argument for this function is only needed to identify master process for writing out the power spectrum file.
	
	int i, bin;
	double binsize;
	double comoving_grid_cell_size, norm;
	double *binbound_low, *binbound_high, *power_in_bin, *bin_averager, *power_in_bin_copy, *bin_averager_copy;
	double k_min, k_max, k_dist;
	int k0, k1, k2, kk0, kk1, kk2;
	int m, m_local;
	
	
	int N0, N1, N2;
	int N0_local, N0_local_start;
	N0=M0;
	N1=M1;
	N2=M2;
	N0_local=local_n0_3d;
	N0_local_start=local_0_start_3d;
	
	double N0d, N1d, N2d;
	N0d=N0;
	N1d=N1;
	N2d=N2;
	
	fftw_plan plan3d;
	
	printf("Task %d starting Power_Spectrum3D_MPI.\n", ThisTask);
	fflush(stdout);
	
	// fftw_mpi_init();
		
	// create plan for in-place forward DFT:
	plan3d=fftw_mpi_plan_dft_3d(M0, M1, M2, data3d, data3d, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	printf("Task %d starting FFTW.\n", ThisTask);
	fflush(stdout);
	
	// compute transforms, in-place, as many times as desired:
	fftw_execute(plan3d);
	
	fftw_destroy_plan(plan3d);

	norm=(boxsize*boxsize*boxsize)/(1000.0*1000.0*1000.0)/((N0d*N1d*N2d)*(N0d*N1d*N2d)); // FFTW normalization according to Donhui Jeong's Dissertation Thesis, Chapter 7.
	for (i=0; i<N0_local*N1*N2; i++) data3d[i]=data3d[i]*conj(data3d[i])*norm; // normalized power spectrum.
	
	
	k_max=N0/2.0; // this is maximal k-vector that fits into Brillouin zone (half the m\
	ax real length because of plus-minus directions).
	k_min=1.0;
	
	// Determines comoving scale of FFT grid:
	comoving_grid_cell_size=boxsize/N0; // WARNING: Works only for cubical grid.
	
	
	// Allocate histogram arrays:
	binbound_low=(double *) malloc(number_of_bins*sizeof(double));
	assert(binbound_low != NULL);

	binbound_high=(double *) malloc(number_of_bins*sizeof(double));
	assert(binbound_high != NULL);
	
	power_in_bin=(double *) malloc(number_of_bins*sizeof(double));
	assert(power_in_bin != NULL);

	power_in_bin_copy=(double *) malloc(number_of_bins*sizeof(double));
	assert(power_in_bin_copy != NULL);

	bin_averager=(double *) malloc(number_of_bins*sizeof(double));
	assert(bin_averager != NULL);

	bin_averager_copy=(double *) malloc(number_of_bins*sizeof(double));
	assert(bin_averager_copy != NULL);
	
	// Zero out power spectrum array:
	
	for (bin=0;bin<number_of_bins;bin++)
	{
		power_in_bin[bin]=0.0;
		power_in_bin_copy[bin]=0.0;
	}
	
	// Zero out bin averager:
	for (bin=0;bin<number_of_bins;bin++)
	{
		bin_averager[bin]=0.0;
		bin_averager_copy[bin]=0.0;
	}
	
	// **Linear bins**:
	// binsize=(k_max-k_min+0.001)/((double) number_of_bins);
	binsize=(k_max-k_min)/((double) number_of_bins);
	// **Log bins**:
	//binsize=(log(k_max)/log(10)-log(k_min)/log(10))/((double) parameters.number_of_bins);
	
	for (bin=0;bin<number_of_bins;bin++)
	{
		binbound_low[bin]=k_min+(bin)*binsize;
		binbound_low[bin]*=(2*M_PI*1000.0/boxsize); //(2*M_PI/comoving_grid_cell_size);
		binbound_high[bin]=k_min+(bin)*binsize+binsize;
		binbound_high[bin]*=(2*M_PI*1000.0/boxsize); // (2*M_PI/comoving_grid_cell_size);
        // WARNING: k should come out in units of h/kpc with that. Seems to be a factor of 1000 off...
    }
	
	
	//printf("Done Fourier transforming 2D image array.\n");
	
	printf("Task %d filling power spectrum bins...\n", ThisTask);
	fflush(stdout);
	
	for (k0=0;k0<N0_local;k0++)
	{
		
		//printf("Task %d: k0=%d\n", ThisTask, k0);
		//fflush(stdout);
		
		
		for (k1=0;k1<N1;k1++)
		{
			
			//printf("Task %d: k1=%d\n", ThisTask, k1);
			//fflush(stdout);
			
			for (k2=0;k2<N2;k2++)
			{
				
				//printf("Task %d: k2=%d\n", ThisTask, k2);
				//fflush(stdout);
				
				
				m=(k0+N0_local_start)*(N1*N2)+k1*N2+k2;
				m_local=k0*(N1*N2)+k1*N2+k2;
				
				if (m!=0) // don't include zero mode in power spectrum
				{
					// Periodic shortening of vectors:
					if (k0+N0_local_start>N0/2) kk0=k0+N0_local_start-N0;
					else kk0=k0+N0_local_start;
					if (k1>N1/2) kk1=k1-N1;
					else kk1=k1;
					if (k2>N2/2) kk2=k2-N2;
					else kk2=k2;
					
					k_dist=sqrt(kk0*kk0+kk1*kk1+kk2*kk2);
				
					// **Linear bins**:
					bin=floor((k_dist-k_min)/binsize);
					// **Log bins**:
					//bin=floor((log(k_dist)/log(10)-log(k_min)/log(10))/binsize);
					if (k_dist<k_max) // if prevents long k-vectors to be written that do not go full circle in square 3D zone (and are thus undercounted).
					{
						power_in_bin_copy[bin]+=data3d[m_local]; // imagearray[m]; // average over all bins with same |k|
						// power_in_bin[i][bin]+=support_active_image[m];
						bin_averager_copy[bin]++; // counts how many spectrum points inserted in bin.
					}
				}
			}
		}
	}
	
	printf("Task %d done filling power spectrum bins.\n", ThisTask);
	fflush(stdout);
		
	MPI_Allreduce(power_in_bin_copy, power_in_bin, number_of_bins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(bin_averager_copy, bin_averager, number_of_bins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	printf("Task %d done with MPI_Allreduce.\n", ThisTask);
	fflush(stdout);
	
	
	// Finish averaging procedure:
	for (bin=0;bin<number_of_bins;bin++) 
	{	
			if (bin_averager[bin]>0) power_in_bin[bin]/=bin_averager[bin]; // average over all k's with same |k|, now divided by number of spectrum points inserted in each bin.
			else power_in_bin[bin]=0.0;
		//power_in_bin[i][bin]/=(4.0*pi*bin*bin); // k-space averaging
	}
	
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Print power spectrum to file (can change this to return the arrays "power_in_bin" and "binvalue" to main program):
		if (ThisTask==0)
		{
			FILE *outfile1;
			outfile1=fopen(filename, "w");
		
			for (bin=0; bin<number_of_bins; bin++)
			{
				fprintf(outfile1, "%e %e %e\n", binbound_low[bin], binbound_high[bin], power_in_bin[bin]);
				// fprintf(outfile1, "%e %e\n", binvalue[bin], (binvalue[bin]*(binvalue[bin]+1.0)/(2*M_PI))*power_in_bin[bin]);
			}
		
			fclose(outfile1);

		}
		
	printf("Task %d done with Power_Spectrum3D_MPI function.\n", ThisTask);
	fflush(stdout);	
	
		free(binbound_low);
		free(binbound_high);
		free(power_in_bin);
		free(power_in_bin_copy);
		free(bin_averager);
		free(bin_averager_copy);
	
	
}




void Power_Spectrum2D(double* imagearray, int nx, int ny, double survey_angle, int number_of_bins, char filename[]) // reads in image array, copies it into Fourier transformable array, Fourier transforms, and then overwrites original image array
{
	int k;
	int nxny;
	nxny=nx*ny;
	
	fftw_complex *array, *ft_array;
	fftw_plan plan_ft;
	
	array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);
	ft_array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);
	assert(array!=NULL);
	assert(ft_array!=NULL);
	
	plan_ft = fftw_plan_dft_2d(ny, nx, array, ft_array, FFTW_FORWARD, FFTW_ESTIMATE);
	
	
	double binsize; // width of a k-bin;
	double* binvalue; // k-value (x-axis) of k-bin with a certain bin number.
	double* power_in_bin; // power in k-bin, averaged over all k-vectors of one realization with \
	same length |k|.
	double* bin_averager;
	double k_max, k_min;
	
	int bin, m, l, k1, k2;
	double k_dist;
	
	k_max=nx/2.0; // this is maximal k-vector that fits into Brillouin zone (half the m\
	ax real length because of plus-minus directions).
	k_min=1.0;

	// Allocate histogram arrays:
	binvalue=(double *) malloc(number_of_bins*sizeof(double));
	assert(binvalue != NULL);
	power_in_bin=(double *) malloc(number_of_bins*sizeof(double));
	assert(power_in_bin != NULL);
	bin_averager=(double *) malloc(number_of_bins*sizeof(double));
	assert(bin_averager != NULL);
	
	// Zero out power spectrum array:

	for (bin=0;bin<number_of_bins;bin++)
	{
		power_in_bin[bin]=0.0;
	}
	
	// Zero out bin averager:
	for (bin=0;bin<number_of_bins;bin++)
	{
		bin_averager[bin]=0.0;
	}
	
	// **Linear bins**:
	binsize=(k_max-k_min+0.001)/((double) number_of_bins);
	// **Log bins**:
	//binsize=(log(k_max)/log(10)-log(k_min)/log(10))/((double) parameters.number_of_bins);
	
	for (bin=0;bin<number_of_bins;bin++)
	{
		binvalue[bin]=k_min+(bin)*binsize;
		binvalue[bin]*=360.0/survey_angle; // unit conversion. 3.6 to be su\
		bstituted by survey angle in degrees. //360.0/(3.6*2.0*pi);//2.0*pi/((3.0/360.0*2.0*pi)/2048.0); // un\
		it conversion from integer bins to k=2*pi/theta.  1=3 deg/2048, need in radians
		// above unit conversion comes from: k_max=parameters.nx/2 in program units; c\
		ompared to k_max=survey_angle/(parameters.nx/2) in physical units to be converted to rad and by k=2*pi\
		/theta, theta angle in deg (half boxside vector in Fourier space corresponds to double spacing in real\
		space.
	}
	//printf("Done writing power bin centers.\n");
	
	
	
	
	printf("Processing file for power spectrum...\n");	
	
	// Read imagearray into Fourier-able array:
	for (k=0;k<nxny;k++) array[k]=imagearray[k];
	
	// Execute forward Fourier transform with FFTW:
	fftw_execute(plan_ft);
	
	// Normalize forward Fourier transform (use this convention here, helps with getting proper normalization of WL power spectrum):
	for (k=0;k<nxny;k++)
	{
		ft_array[k]/=((double) nxny); // performing FT normalization 1/L^2 on forward Fourier transformation, where L is length of map side.
	}
		
	// Compute F(k) x F*(k) = |F(k)|^2 (needed for power spectrum):
	for (k=0;k<nxny;k++) ft_array[k]=ft_array[k]*conj(ft_array[k]);
	
	// Overwrite original real imagearray with Fourier transformed squared array:
	for (k=0;k<nxny;k++)
	{
		imagearray[k]=((double *) ft_array)[2*k]; // reading only real part from complex array "ft_array" into real array "imagearray".
		//imagearray[k]/=((double) parameters.nxny*parameters.nxny); // old line, already done above // performing FT normalization 1/L^2 on forward Fourier transformation, do it squared, because |F|^2 contains two F's in array.
	}
	//or: return ft_array; (then function: fftw_complex* Fourier_forward{...}).
	
	// Free objects needed for FFTW (clean up before exiting function):
	fftw_destroy_plan(plan_ft);
	fftw_free(array);
	fftw_free(ft_array);

	fftw_cleanup(); // free any residual information kept by FFTW after destroying plans and freeing arrays.

	// Now normalize Power spectrum according to convention (takes into account normalization factor between continuous definition of power spectrum an discrete transform:
	// Power Spectrum P(k) normalized in the continuous case as <d(k)d(k)>=(2*Pi)^2 P(k).
	double angle_norm;
	angle_norm=(2*M_PI*survey_angle/360.0)*(2*M_PI*survey_angle/360.0);
	for (m=0; m<nxny; m++)
	{
		imagearray[m]*=angle_norm;
	}

	//printf("Done Fourier transforming 2D image array.\n");
	
	for (k=0;k<ny;k++)
	{
		for (l=0;l<nx;l++)
		{
			m=k*nx+l;
			
			if (m!=0) // don't include zero mode in power spectrum
			{
				// Periodic shortening of vectors:
				if (l>nx/2) k1=l-nx;
				else k1=l;
				if (k>ny/2) k2=k-ny;
				else k2=k;
				
				k_dist=sqrt(k1*k1+k2*k2);
				// **Linear bins**:
				bin=floor((k_dist-k_min)/binsize);
				// **Log bins**:
				//bin=floor((log(k_dist)/log(10)-log(k_min)/log(10))/binsize);
				if (k_dist<=k_max) // if prevents long k-vectors to be written that do not go full circle in square 2D zone (and are thus undercounted).
				{
					power_in_bin[bin]+=imagearray[m]; // average over all bins with same |k|
					// power_in_bin[i][bin]+=support_active_image[m];
					bin_averager[bin]++; // counts how many spectrum points inserted in bin.
				}
			}
		}
	}
	// Finish averaging procedure:
	for (bin=0;bin<number_of_bins;bin++) power_in_bin[bin]/=bin_averager[bin]; // average over all k's with same |k|, now divided by number of spectrum points inserted in each bin.
	//power_in_bin[i][bin]/=(4.0*pi*bin*bin); // k-space averaging
	
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Print power spectrum to file (can change this to return the arrays "power_in_bin" and "binvalue" to main program):
	FILE *outfile1;
	outfile1=fopen(filename, "w");
	
	for (bin=0; bin<number_of_bins; bin++)
	{
		fprintf(outfile1, "%e %e\n", binvalue[bin], power_in_bin[bin]);
		// fprintf(outfile1, "%e %e\n", binvalue[bin], (binvalue[bin]*(binvalue[bin]+1.0)/(2*M_PI))*power_in_bin[bin]);
	}
	
	fclose(outfile1);
	

	free(binvalue);
	free(power_in_bin);
	free(bin_averager);
	
	
}



