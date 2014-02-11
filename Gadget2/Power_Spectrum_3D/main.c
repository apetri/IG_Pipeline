/*
 *  main.c
 *  3D Power Spectrum Calculator (and FFTW3-MPI Test)
 *
 *  Created by Jan Michael Kratochvil on 11/02/2012.
 *  Copyright 2012 Jan Michael Kratochvil. All rights reserved.
 *
 */




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <complex.h>
//#include <fftw3.h>
#include <fftw3-mpi.h>
#include <mpi.h>

#include "main.h"
#include "snapshot_read.h"
#include "endianness.h"
#include "particle_processing.h"
#include "power_spectrum.h"

#include "fitsio.h"
#include "fits.h"

#include "ini.h"
#include "options.h"


int machine_endian;
int snapshot_endian;

int superrank, supersize;


int main (int argc, char ** argv) {

	if(argc<4){
		fprintf(stderr,"Usage %s <ini_options_file> <snapshot_number> <number_of_files_per_snapshot>\n",*argv);
		exit(1);
	}
	
	char series_foldername[1000],snapshot_path[1000], snapshot_filenamebase[1000], snapshot_filename[1000], power_spectrum_path[1000], power_spectrum_filenamebase[1000], power_spectrum_filename[1000];
	
	//////////////////////////////////////////////////////////////////////////////
	// MANUAL SECTION:
	// NOTE: Need to read quantities below from a parameter file. eventually:
	//AP: Done
	//////////////////////////////////////////////////////////////////////////////

	//Parse options from file
	sys_options *options = malloc(sizeof(sys_options));
	if(ini_parse(argv[1],handler,options)<0){
		fprintf(stderr,"ini options file %s not found\n",argv[1]);
		exit(1);
	}


	int snapshot_number=atoi(argv[2]);
	int number_of_files_per_snapshot=atoi(argv[3]);
	int particle_side=options->num_particles_side; // number of particles in N-body simulation along one dimension. 
	
	// FFT Grid:
	int FFT_grid=options->FFT_grid_size;
	
	// Number of power spectrum bins"
	int number_of_bins=options->number_of_bins;
	
	// Length of particle position buffer (in number of particles; each particle will have three coordinates):
	int particle_buffer_length=options->particle_buffer_length;

	// Gadget-2 snapshot path and filename:
	
	sprintf(series_foldername,"%s/Storage/sims/snapshots/%s-series",options->mass_storage_path,options->series_name);
	sprintf(snapshot_path,"%s/%s-%db%d_%s_ic%d",series_foldername,options->series_name,particle_side,options->box_size_snapshot_kpc,options->model_basename,options->realization_number);
	
	sprintf(snapshot_filenamebase, "snapshot");
	sprintf(power_spectrum_path, "%s",options->power_spectrum_savepath);
	sprintf(power_spectrum_filenamebase, "%s", options->power_spectrum_filebase);
	
	// Project to 2D (just as output check):
	char projection_filename[1000];
	long *naxes;
	naxes=malloc(2*sizeof(long));
	sprintf(projection_filename, "%s/%s",options->projection_savepath,options->projection_filebase);
	int projection_axis=2;
	int min_cell=0; // Determines start of slice in projected direction; in FFT grid cell units (integer), not in comoving distance.
	int max_cell=FFT_grid; // Determines end of slice in projected direction.
	//////////////////////////////////////////////////////////////////////////////

	//Don't need options anymore
	free(options);
	
	int i, file_number;
	double boxsize;
	
	int ThisTask, NTasks;
	
	double *imagearray;
	
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
	MPI_Comm_size( MPI_COMM_WORLD, &NTasks );
	
	superrank=ThisTask;
	supersize=NTasks;
	
	if (ThisTask==0) printf("Starting 3D matter power spectrum calculation code...\n");
	fflush(stdout);
	
	
	// Declare and allocate some variables and arrays necessary for the FFTW3-MPI calculation of the power spectrum:
	
	// const 
	ptrdiff_t M0, M1, M2;
	
	M0=FFT_grid;
	M1=FFT_grid;
	M2=FFT_grid;
	
	fftw_plan plan3d; //, plan3d_back;
	fftw_complex *data3d;
	ptrdiff_t alloc_local_3d, local_n0_3d, local_0_start_3d, ii, jj, kk;
	
	fftw_mpi_init();
	
	// get local data size and allocate:
	alloc_local_3d=fftw_mpi_local_size_3d(M0, M1, M2, MPI_COMM_WORLD, &local_n0_3d, &local_0_start_3d);
	data3d=fftw_alloc_complex(alloc_local_3d);
	
	// create plan for in-place forward DFT:
	//// plan3d=fftw_mpi_plan_dft_3d(M0, M1, M2, data3d, data3d, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
	// plan3d_back=fftw_mpi_plan_dft_3d(M0, M1, M2, data3d, data3d, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	int N0, N1, N2;
	int N0_local, N0_local_start;
	
	N0=M0;
	N1=M1;
	N2=M2;
	N0_local=local_n0_3d;
	N0_local_start=local_0_start_3d;
	
	// Process zero collects all boundaries:
	// MPI_Gather(N0_locals, MPI_INT, 1, N0_local, MPI_INT, 1, 0, MPI_COMM_WORLD);
	// MPI_Gather(N0_local_starts, MPI_INT, 1, N0_local_start, MPI_INT, 1, 0, MPI_COMM_WORLD);
	// perhaps send the collections to every process if needed.
	
	// Zero out Fourier grid array before filling it:
	for (i=0; i<N0_local*N1*N2; i++) data3d[i]=0.0;
		
	// Loop over all files comprising a adget-2 snapshot at a fixed redshift:
	for (file_number=0; file_number<number_of_files_per_snapshot; file_number++)
	{
		
		// Construct Gadget-2 snapshot filename:
		if(number_of_files_per_snapshot>1){
			sprintf(snapshot_filename, "%s/%s_%03d.%d", snapshot_path, snapshot_filenamebase, snapshot_number, file_number);
		} else{
			sprintf(snapshot_filename, "%s/%s_%03d", snapshot_path, snapshot_filenamebase, snapshot_number);
		}
		
		// Read particles from Gadget-2 snapshot file and insert them into FFT grid:
		if (ThisTask==0) boxsize=read_particles_master(snapshot_filename, particle_buffer_length, data3d, N0_local, N0_local_start, N0, N1, N2);
		else boxsize=read_particles_slave(snapshot_filename, particle_buffer_length, data3d, N0_local, N0_local_start, N0, N1, N2);
	}
	
	// Convert number of particles per grid cell into density contrast:
	convert_particles_to_density_contrast(data3d, N0_local, N0, N1, N2, particle_side);

	 
	
	// Print data array on screen before any transformations:
	if (ThisTask==0) printf("3D Data array before FT transformations:\n");
	//show3d(data3d, local_0_start_3d, local_n0_3d, M1, M2, ThisTask, NTasks);
	// Do the same for the MPI collected array:
	//MPI_Gather(data, (int) local_n0*N1, MPI_COMPLEX, alldata, (int) local_n0*N1, MPI_COMPLEX, 0, MPI_COMM_WORLD);

	if (projection_axis==0) imagearray=malloc(N1*N2*sizeof(double));
	else if (projection_axis==1) imagearray=malloc(N0*N2*sizeof(double));
	else if (projection_axis==2) imagearray=malloc(N0*N1*sizeof(double));
	assert(imagearray!=NULL);
	
	project_to_2D(data3d, N0_local, N0_local_start, N0, N1, N2, imagearray, projection_axis, min_cell, max_cell, naxes);
	if (ThisTask==0) writeFITSimage_f(projection_filename, 2, naxes, imagearray);
	free(imagearray);
	
	double sum_local=0;
	double sum=0;
	// Sum over all grid cells to see if equal to the number of particles inserted:
	for (i=0; i<N0_local*N1*N2; i++) sum_local+=data3d[i];
	
	MPI_Reduce(&sum_local, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if (ThisTask==0) printf("Sum over all grid cells: %e\n", sum);
	
	MPI_Barrier(MPI_COMM_WORLD);
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	sprintf(power_spectrum_filename, "%s/%s_%03d.txt", power_spectrum_path, power_spectrum_filenamebase, snapshot_number);
	Power_Spectrum3D_MPI(data3d, M0, M1, M2, local_n0_3d, local_0_start_3d, boxsize, number_of_bins, power_spectrum_filename, ThisTask);
	
	
	// compute transforms, in-place, as many times as desired:
	// fftw_execute(plan3d);
	
	// fftw_destroy_plan(plan3d);
	
	
	// MPI_Gather(alldata,);
	// Print data array on screen after forward Fourier transformations:
	// if (ThisTask==0) printf("3D Data array after forward FT transformations:\n");
	// show3d(data3d, local_0_start_3d, local_n0_3d, M1, M2, ThisTask, NTasks);
	
	
	/*
	fftw_execute(plan3d_back);
	
	fftw_destroy_plan(plan3d_back);
	
	// Normalize Fourier Transfrom:
	for (ii=0; ii<local_n0_3d; ii++)
	{
		for (jj=0; jj<M1; jj++)
		{
			for (kk=0; kk<M2; kk++)
			{
				data3d[ii*M1*M2+jj*M2+kk]/=(M0*M1*M2);
			}
		}
	}
	
	// Print data array on screen after backward Fourier transformations:
	if (ThisTask==0) printf("3D Data array after backward FT transformations:\n(Should equal the original array, up to a normalization constant (FFTW needs to be normalized by hand.)\n");
	show3d(data3d, local_0_start_3d, local_n0_3d, M1, M2, ThisTask, NTasks);
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	*/
	
	
	
	fftw_mpi_cleanup();
	
	MPI_Finalize();

    if (ThisTask==0) printf("Done with FFTW3-MPI Test. Program ran successfully.\n");
    return 0;
}


/*
// Arbitrary function, filling array to be Fourier transformed with initial values:
fftw_complex my_function(int i, int j)
{
	fftw_complex x;
	
	x=i*100+j;
	
	return x;
}

fftw_complex my_function3d(int i, int j, int k)
{
	fftw_complex x;
	
	// x=i*10000+j*100+k;
	x=cos(i+j+k);
	
	return x;
}
 */

// Print array on screen, one process at a time (without combining on one process):
void show(fftw_complex *data, ptrdiff_t local_0_start, ptrdiff_t local_n0, ptrdiff_t N1, int ThisTask, int NTasks)
{
	int k, l, m, localindex, index_0, index_1;
	double x, y;

	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (k=0;k<NTasks;k++)
	{
		
		if (ThisTask==k)
		{
			for (l=0; l<local_n0; l++)
			{
				for (m=0; m<N1; m++)
				{
					index_0=(local_0_start+l);
					index_1=m;
					localindex=l*N1+m;
					x=creal(data[localindex]);
					y=cimag(data[localindex]);
							
					printf("data[%d,%d]= %e + I * %e   ", index_0, index_1, x, y); 
				}
				printf("\n");
			}
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
	}
}



void show3d(fftw_complex *data, ptrdiff_t local_0_start_3d, ptrdiff_t local_n0_3d, ptrdiff_t M1, ptrdiff_t M2, int ThisTask, int NTasks)
{
	int k, l, m, n, localindex, index_0, index_1, index_2;
	double x, y;
	
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (k=0;k<NTasks;k++)
	{
		
		if (ThisTask==k)
		{
			for (l=0; l<local_n0_3d; l++)
			{
				for (m=0; m<M1; m++)
				{
					for (n=0; n<M2; n++)
					{
						index_0=(local_0_start_3d+l);
						index_1=m;
						index_2=n;
						localindex=l*(M1*M2)+m*M2+n;
						x=creal(data[localindex]);
						y=cimag(data[localindex]);
					
					printf("data[%d,%d,%d]= %e + I * %e   ", index_0, index_1, index_2, x, y); 
	
					}
					printf("\n");
				}
				printf("\n-----------------\n");
			}
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
	}
}


// Projects FFT 3D grid onto 2D grid along one coordinate axis:
void project_to_2D(fftw_complex *data3d, int N0_local, int N0_local_start, int N0, int N1, int N2, double *imagearray, int projection_axis, int min_cell, int max_cell, long *naxes)
{
	int i, j, k, index, m, insert, imagearray_size=0;
	int insertparticle=0;

//<AP>
	double *imagearray_copy;
//</AP>
	
	if (projection_axis==0)
	{	
		imagearray_size=N1*N2;
		naxes[0]=N2; // WARNING: May have to change the order of these two.
		naxes[1]=N1;
	}
	else if (projection_axis==1)
	{	
		imagearray_size=N0*N2;
		naxes[0]=N2;
		naxes[1]=N0;	
	}
	else if (projection_axis==2)
	{	
		imagearray_size=N0*N1;
		naxes[0]=N1;
		naxes[1]=N0;
	}

//<AP>
	imagearray_copy = malloc(sizeof(double)*imagearray_size);
//</AP>

	for (i=0; i<imagearray_size; i++){ 
		imagearray[i]=0.0;
		imagearray_copy[i]=0.0;
	}
	
	printf("projaxis: %d, min_cell %d, max_cell %d, N0_local %d, N0_start_local %d, N0 %d, N1 %d, N2 %d.\n", projection_axis, min_cell, max_cell, N0_local, N0_local_start, N0, N1, N2);
	
	for (i=0; i<N0_local; i++)
	{
		for (j=0; j<N1; j++)
		{
			for (k=0; k<N2; k++)
			{
				insert=1;
				
				//printf("projaxis: %d, min_cell %d, max_cell %d, N0_local %d, N0_start_local %d, N0 %d, N1 %d, N2 %d, i, j, k: (%d, %d, %d).\n", projection_axis, min_cell, max_cell, N0_local, N0_local_start, N0, N1, N2, i, j, k);
				
				if (projection_axis==0 && min_cell<=i+N0_local_start && max_cell>i+N0_local_start)
				{
					m=j*N2+k;
					//printf("Went through here 0.\n");
				}
				else if (projection_axis==1 && min_cell<=j && max_cell>j)
				{
				 	m=(i+N0_local_start)*N2+k;
					//printf("Went through here 1.\n");
				}
				else if (projection_axis==2 && min_cell<=k && max_cell>k)
				{
				 	m=(i+N0_local_start)*N1+j;
					//printf("Went through here 2.\n");
				}
				else 
				{	
					//printf("Went through here: none.\n");
					insert=0;
				}
					
				if (insert!=0)
				{	
					
					index=i*N1*N2+j*N2+k;	
					imagearray_copy[m]+=creal(data3d[index]);
					insertparticle++;
					if (insertparticle<10) printf("ThisTask: %d -- i, j, k: (%d, %d, %d), m: %d, index: %d, value: %e\n", superrank, i, j, k, m, index, creal(data3d[index])); 
				}
			}
		}
	}
	
	printf("Inserted %d cells\n", insertparticle);
	
	printf("Before MPI_Reduce: ThisTask: %d, imagearray: %e %e %e %e %e\n", superrank, imagearray_copy[0], imagearray_copy[1], imagearray_copy[2], imagearray_copy[3], imagearray_copy[4]);
	
	// MPI_Reduce(imagearray, imagearray, imagearray_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	// WRONG: Send and receive buffer must be different in this case.
	// AP: EXTREMELY WRONG! MPI_Allreduce with same buffer addresses causes crash at runtime!!! 
	MPI_Allreduce(imagearray_copy, imagearray, imagearray_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//<AP>
	free(imagearray_copy);
//</AP>
	
	printf("After MPI_Reduce: ThisTask: %d, imagearray: %e %e %e %e %e\n", superrank, imagearray[0], imagearray[1], imagearray[2], imagearray[3], imagearray[4]);
}


