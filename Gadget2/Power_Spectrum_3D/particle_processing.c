/*
 *  particle_processing.c
 *  Gadget-2 Snapshot File Reader (low memory version with buffer)
 *
 *  Created by Jan Michael Kratochvil - on 11/03/2012.
 *  Copyright 2012 Jan Michael Kratochvil. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <time.h>

#include <complex.h>
#include <mpi.h>
#include <fftw3-mpi.h>

#include "main.h"
#include "snapshot_read.h"
#include "endianness.h"
#include "particle_processing.h"






//int main(int argc, char ** argv)

double read_particles_master(char snapshot_filename[], int particle_buffer_length, fftw_complex *data3d, int N0_local, int N0_local_start, int N0, int N1, int N2)
{
	
	//if (counter==1) snapshot_endian=snapshotEndianness(checkfilename, file_number, snapshot_number);
	//load_snapshot_multi(checkfilename, file_number, snapshot_number);	
	
	FILE *snapshot_file;
	// char snapshot_path[1000], snapshot_filenamebase[200], snapshot_filename[1000];
	// int particle_buffer_length; 
	int particles_read_from_file, particles_read_this_time, particles_in_file, allparticles_read, snapshot_number, file_number;
	int particle_type, alltype_particles_read_from_file; // variables spanning several particle types.
	struct io_header_1 header1;
	float *particle_position_buffer;
	int readcycle, number_of_readcycles;
	double boxsize;
	double grid_cell_size[3], low_x, high_x;
	
	// Determine machine endianness (do this only once at the begin of a run):
	machine_endian=machineEndianness();
	printf("Endianness of Machine (0: Little-Endian, 1: Big-Endian): %d\n", machine_endian);
	fflush(stdout);
	
	// Determine snapshot endianness:
	snapshot_endian=snapshotEndianness(snapshot_filename);
	// Open snapshot file and read header (keeps snapshot file open):
	open_snapshot_multi(&snapshot_file, snapshot_filename, &header1);
	
	// Print header info on screen:
	show_header(header1);
	
	// Transfer header information and communicate relevant quantities:
	boxsize=header1.BoxSize;
	MPI_Bcast(&boxsize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// ... (need to implement this for all header variables used by code).
	// Compute derived quantities:
	get_grid_cell_size(grid_cell_size, boxsize, N0, N1, N2);
	get_local_comoving_bounds(N0_local, N0_local_start, grid_cell_size[0], &low_x, &high_x);
	
	// Allocate particle position buffer:
	particle_position_buffer=malloc(3*particle_buffer_length*sizeof(float));
	assert(particle_position_buffer!=NULL);
	
	// No particles read yet:
	alltype_particles_read_from_file=0;
	particles_read_this_time=0;
	allparticles_read=0;

	// Read all particles from file, one buffer filling at a time:
	
	read_skip(&snapshot_file); // skipping Fortran block dummy.
	
	// Loop over all particle types (gas, CDM, etc.). For pure CDM, only particle_type=1 gets executed).
	for(particle_type=0;particle_type<6;particle_type++)
	{
		particles_in_file=header1.npart[particle_type]; // how many particles are in the file of one type.
		particles_read_from_file=0; // particles of one type read from file.
		
		MPI_Bcast(&particles_in_file, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		number_of_readcycles=particles_in_file/particle_buffer_length;
		if (particles_in_file%particle_buffer_length!=0) number_of_readcycles++;
		
		//if (particle_type==1)
		{
			printf("particles_in_file: %d, particle_buffer_length %d, number_of_cycles %d\n", particles_in_file, particle_buffer_length, number_of_readcycles);
		}
		
		// Manual:
		//number_of_readcycles=1;
		
		MPI_Bcast(&number_of_readcycles, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		for (readcycle=0; readcycle<number_of_readcycles; readcycle++)
		// while (particles_in_file>0 && allparticles_read==0) // check if all particles of one type read.
		{
			
			particles_read_this_time=0; // zero out just as a precaution.
			// Fill particle position buffer with new particles from snapshot file: 
			allparticles_read=read_particle_positions(&snapshot_file, particle_position_buffer, &particles_read_from_file, &particles_read_this_time, particles_in_file, particle_buffer_length);
			MPI_Bcast(&particles_read_this_time, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
			// printf("particles_read_this_time: %d\n", particles_read_this_time);
			
			// Check integrity of read-in buffer:
			check_snapshot_multi(particle_position_buffer, particle_buffer_length, particles_read_this_time, header1.BoxSize);
		
			// show_particle_position_buffer(particle_position_buffer, particle_buffer_length);
		
			// Do some work here with each particle buffer filling (may distribute onto multiple processors):
			if (particle_type==0) // gas particles.
			{
				// nothing implemented for gas particles yet.
			}
			else if (particle_type==1) // CDM particles.
			{
				// here code something up what will do with 
				// process_CDM_particles(particle_position_buffer, particle_buffer_length, particles_read_this_time);
				
				// Send particles to all processors (they will do their filtering themselves):
				MPI_Bcast(particle_position_buffer, 3*particle_buffer_length, MPI_FLOAT, 0, MPI_COMM_WORLD);
				insert_particles_in_grid(particle_position_buffer, particles_read_this_time, data3d, N1, N2, grid_cell_size, low_x, high_x);
				// This assumes that the communications overhead is much smaller than the calculation overhead which may not be true.
			}
			// ...
				
		}
		alltype_particles_read_from_file+=particles_read_from_file;
		
	}
	
	
	read_skip(&snapshot_file); // skipping Fortran block dummy.
									
	// Check if all particles from snapshot read:
	check_snapshot_complete(particles_read_from_file, particles_in_file);
									
	
	close_snapshot_multi(&snapshot_file); // close snapshot file after reading particle positions (never read velocities, ID, etc., here - saves I/O bandwidth).
	
	free(particle_position_buffer);
	
	printf("In total, %d particles processed by this process.\n", alltype_particles_read_from_file);

	printf("Task: %d got to here.\n", superrank);
	fflush(stdout);
	
	// Do some more work after all the particles were read in and distributed/processed after this function exits.
	
	return boxsize;
	
}

double read_particles_slave(char snapshot_filename[], int particle_buffer_length, fftw_complex *data3d, int N0_local, int N0_local_start, int N0, int N1, int N2)
{
	float *particle_position_buffer;
	int particle_type, particles_in_file, number_of_readcycles, readcycle, particles_read_this_time;
	double boxsize, grid_cell_size[3], low_x, high_x;
	
	// Allocate particle position buffer:
	particle_position_buffer=malloc(3*particle_buffer_length*sizeof(float));
	assert(particle_position_buffer!=NULL);
	
	// Receive all relevant header variables:
	MPI_Bcast(&boxsize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	get_grid_cell_size(grid_cell_size, boxsize, N0, N1, N2);
	get_local_comoving_bounds(N0_local, N0_local_start, grid_cell_size[0], &low_x, &high_x);
	
	// Receive particle buffer contents and execute calculation:
	for(particle_type=0;particle_type<6;particle_type++)
	{
		MPI_Bcast(&particles_in_file, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&number_of_readcycles, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		for (readcycle=0; readcycle<number_of_readcycles; readcycle++)
		{
			MPI_Bcast(&particles_read_this_time, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
			if (particle_type==0) // gas particles.
			{
				// nothing implemented for gas particles yet.
			}
			else if (particle_type==1) // CDM particles.
			{
				MPI_Bcast(particle_position_buffer, 3*particle_buffer_length, MPI_FLOAT, 0, MPI_COMM_WORLD);
				insert_particles_in_grid(particle_position_buffer, particles_read_this_time, data3d, N1, N2, grid_cell_size, low_x, high_x);
			}
		}
	}
	
	free(particle_position_buffer);
	
	return boxsize;
}
	


void get_grid_cell_size(double *grid_cell_size, double boxsize, int N0, int N1, int N2)
{
	
	grid_cell_size[0]=boxsize/N0;
	grid_cell_size[1]=boxsize/N1;	
	grid_cell_size[2]=boxsize/N2;
	
}

void get_local_comoving_bounds(int N0_local, int N0_local_start, double grid_cell_size, double *low_x, double *high_x)
{
	*low_x=N0_local_start*grid_cell_size;
	*high_x=(N0_local_start+N0_local)*grid_cell_size;
}


void insert_particles_in_grid(float *buffer, int particles_read_this_time, fftw_complex *data3d, int N1, int N2, double *grid_cell_size, double low_x, double high_x)
{
	int i, index;
	int a, b, c; // integer coordinates on grid.
	int inserted_particles=0;
	
	// printf("Task %d: grid_cell_size %e, low_x %e, high_x %e.\n", superrank, grid_cell_size[0], low_x, high_x);
	
	for (i=0; i<particles_read_this_time; i++)
	{
		index=3*i;
		if (buffer[index]>=low_x && buffer[index]<high_x) // if this is true, particle is in proper range for this process to fill it into its grid.
		{
			a=((int) floor((buffer[index]-low_x)/grid_cell_size[0]));
			b=((int) floor((buffer[index+1])/grid_cell_size[1]));
			c=((int) floor((buffer[index+2])/grid_cell_size[2]));
			
			data3d[a*(N1*N2)+b*N2+c]+=1.0;
			inserted_particles++;
		}
		
	}
	// printf("Task %d inserted %d particles (of %d in buffer).\n", superrank, inserted_particles, particles_read_this_time);
	
}


// This function converts number of particles to density contrast:
void convert_particles_to_density_contrast(fftw_complex *data3d, int N0_local, int N0, int N1, int N2, int particle_side)
{
	int i;
	double N0d, N1d, N2d, particle_side_d, factor;
	
	N0d=(double) N0;
	N1d=(double) N1;
	N2d=(double) N2;
	particle_side_d=(double) particle_side;
	
	factor=(N0d/particle_side_d)*(N1d/particle_side_d)*(N2d/particle_side_d);
	for (i=0; i<N0_local*N1*N2; i++) data3d[i]=(data3d[i]*factor)-1.0;
	
}
