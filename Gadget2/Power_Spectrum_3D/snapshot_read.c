/*
 *  snapshot_read.c
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

// #define SKIP fread(&dummy, sizeof(dummy), 1, fd);
// #define SKIP fread(&dummy, sizeof(dummy), 1, *fd);



// Open Gadget-2 snapshot file (in Gadget's default binary file format) and read header:
//int load_snapshot_multi(char *fname, int file_number, int snapshot_number, int plane_number)
int open_snapshot_multi(FILE **fd, char *fname, struct io_header_1 *header1)
{
	char   buf[1000];
	int    k,dummy,ntot_withmasses;
	
  	int NumPart;
	int Ngas;
	int Nhalo;
	double Time;
	double Redshift;

	int NumPartThisFile;
	int NgasThisFile;
	int NhaloThisFile;
	
	
		sprintf(buf,"%s",fname);
		
		if(!(*fd=fopen(buf,"r")))
		{
			printf("can't open file `%s`\n",buf);
			exit(0);
		}
		
		printf("reading `%s' ...\n",buf); fflush(stdout);
		
		fread(&dummy, sizeof(dummy), 1, *fd);
		printf("First Dummy: %d\n", dummy);
		// header_fread(&header1, sizeof(header1), 1, fd);
		header_fread(header1, sizeof(*header1), 1, *fd);
		fread(&dummy, sizeof(dummy), 1, *fd);
		printf("Second Dummy: %d\n", dummy);
		printf("Size of header1 %d\n", ((int) sizeof(*header1)));
		
		printf("w0: %e \n", header1->w0);
		printf("wa: %e \n", header1->wa);
		printf("comoving distance: %e \n", header1->comoving_distance);
		fflush(stdout);
		
		int test_int;
		//int_fread(&test_int, sizeof(int), 1, fd);
		//int_fread(&test_int, sizeof(int), 1, fd);
		//int_fread(&test_int, sizeof(int), 1, fd);
		//test_int=IntSwap(test_int);
		//printf("Test int: %d\n", test_int);
		
		printf("Header read: %d, %e, %e, %e\n", header1->npart[1], header1->mass[1], header1->time, header1->redshift);
		fflush(stdout);
		
		/*      if(files==1)
		 {
		 for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
		 NumPart+= header1.npart[k];
		 Ngas= header1.npart[0];
		 Nhalo= header1.npart[1];
		 }
		 else
		 */	{
			 for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
				 NumPart+= header1->npartTotal[k];
			 Ngas= header1->npartTotal[0];
			 Nhalo= header1->npartTotal[1];
			 
			 for(k=0, NumPartThisFile=0, ntot_withmasses=0; k<5; k++)
				 NumPartThisFile+= header1->npart[k];
			 NgasThisFile= header1->npart[0];
			 NhaloThisFile= header1->npart[1];
		 }
		
		for(k=0, ntot_withmasses=0; k<5; k++)
		{
			if(header1->mass[k]==0)
				ntot_withmasses+= header1->npart[k];
		}
		
		
		// Done reading header. Transfer of header variables into general code parameters will happen in separate function. Particle positions will also be read by different function.
		
}

// Close Gadget-2 snapshot file:
void close_snapshot_multi(FILE **fd)
{
	fclose(*fd); // closes reading of snapshot file.
	// can add here other functions as well, like freeing of arrays, etc.
}

		
// This function fills particle buffer by reading from Gadget-2 snapshot file (which must have been already previously opened):
// (This function must be called consecutively multiple times until all particles from a snapshot file are read (the particle position buffer is typically smaller than the number of partcles in the file).)
// Structure of particle listing in Gadget-2 defaut binary file format:
// First positions, then velocities, then IDs, then species masses
// within those, k from 0 to 5 cycles through particle types,
// with number of particles for every species header1.npart[k] in that particular file.
// The three coordinates of a particle are always listed adjacently.
int read_particle_positions(FILE **fd, float *particle_position_buffer, int *particles_read_from_file, int *particles_read_this_time, int particles_in_file, int particle_buffer_length) 		
{
	// Here read particle positiond into particle buffer:
	
	int dummy, i, all_particles_read;
	
	if (*particles_read_from_file+particle_buffer_length<particles_in_file)
	{
		*particles_read_this_time=particle_buffer_length;
		all_particles_read=0;
	}
	else
	{
		*particles_read_this_time=particles_in_file-*particles_read_from_file; // if not enough particles left, don't fill whole particle read buffer.
		all_particles_read=1;
	}	
	
	// Read particle positions:
	float_fread(particle_position_buffer, sizeof(float), *particles_read_this_time*3, *fd);
	
	*particles_read_from_file+=*particles_read_this_time;
	
	
	// Zero out unused rest of buffer just as a precaution:
	for (i=*particles_read_this_time; i<particle_buffer_length; i++)
	{	
		particle_position_buffer[3*i]=0.0;
		particle_position_buffer[3*i+1]=0.0;
		particle_position_buffer[3*i+2]=0.0;
	}
	
	
	return all_particles_read; // returns 0 if still more particles left to read, or 1 if all particles in the file have been read.
	
}	

void read_skip(FILE **fd)
{
	int dummy;
	fread(&dummy, sizeof(dummy), 1, *fd); // skipping Fortran block dummy.
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Some functions to verify read-in of Gadget-2 snapshot file:

// Check if values read into buffer are consistent with basic parameters (checks for buffer corruption):
void check_snapshot_multi(float *buffer, int buffer_length, int particles_read_this_time, double boxsize)
{	
	int i, k, particle;
	int overstep, understep, bufferstep, error, error_buffer;
	double maxpos, minpos;
	
	overstep=0; understep=0; error=0; error_buffer=0;
	maxpos=-pow(10,300); minpos=pow(10,300);
	
	
	for (particle=0; particle<buffer_length; particle++) // loop over particles in buffer.
	{
		
		for (k=0; k<3; k++) // loop over the three spatial coordinates
		{
			i=3*particle+k;
			
			if (particle<particles_read_this_time)
			{
				if (buffer[i]>maxpos) maxpos=buffer[i]; // determines maximal particle coordinate value.
				if (buffer[i]<minpos) minpos=buffer[i]; // determines minimal particle coordinate value.
				
				
				if (buffer[i]>boxsize)
				{
					if (error==0) printf("ERROR: A particle lies outsize the box. Comoving particle coordinate %d of particle %d in buffer too big: %e .\n", k, particle, buffer[i]);
					error=1;
					overstep++;
				}
				if (buffer[i]<0.0)
				{
					if (error==0) printf("ERROR: A particle lies outsize the box. Comoving particle coordinate %d of particle %d in buffer too small: %e .\n", k, particle, buffer[i]);
					error=1;
					understep++;
				}				
				
				// Check for suspicious particles which lie precisely on box boundary (not impossible to happen, but should have zero likelihood measure of happening):
				if (buffer[i]==0.0 || buffer[i]==boxsize) printf("Suspicious particle: buffer index, particle, position, suspicious coordinate direction: %d, %d, (%e, %e, %e), %d.\n", i, particle, buffer[3*particle], buffer[3*particle+1], buffer[3*particle+2], k);
				
				
			}
			else // check that part of buffer without particles is really zero:
			{
				if (fabs(buffer[i])>1e-15)
				{
					if (error_buffer==0) printf("ERROR: Buffer not zero where it does not contain any particle anymore: buffer entry %d (particle %d, coord. %d) is %e instead of zero.\n", i, particle, k, buffer[i]);
					error_buffer=1;
					bufferstep++;
					
				}
				
			}
			
		} // end loop over k.
		
	}	// end loop over particle.
	
	if (error==1) printf("The following number of particle coordinate violations occurred (too large, too small): %d %d .\n", overstep, understep);
	// else printf("Particle coordinates all within BoxSize. No particle coordinate box boundary violations occurred.\n");
	
	// printf("Minimal and maximal particle position coordinate values were (allowed box boundaries): %e, %e, (%e, %e).\n", minpos, maxpos, 0.0, boxsize);
	
	if (error>0 || error_buffer>0)
	{
		printf("ERROR: Check check_snapshot_multi failed. Something is wrong with the particle positions read from the snapshot file (either a read-in error or the snapshot file is corrupted)\n");
		exit(1);
	}
}


// Check if number of particles read into buffer is equal to expected number of particles in that one snapshot file:
void check_snapshot_complete(int particles_read_from_file, int particles_in_file)
{
	if (particles_read_from_file!=particles_in_file)
	{
		printf("Error: Done with reading Gadget-2 snapshot file, but number of particles inserted into particle position buffer does not match expected number of particles for that file.\nExpected: %d. Actually read in: %d.\nAborting.\n", particles_in_file, particles_read_from_file);
		exit(1);
	}
	// Can add more sophisticated checks here.
	
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Some printing utilities for the screen:


// Print header content on screen:
void show_header(struct io_header_1 header1)
{
	
	printf("Header info:\n");
	
	printf("npart[6]: %d %d %d %d %d %d .\n", header1.npart[0], header1.npart[1], header1.npart[2], header1.npart[3], header1.npart[4], header1.npart[5]);
	printf("mass[6]: %e %e %e %e %e %e .\n", header1.mass[0], header1.mass[1], header1.mass[2], header1.mass[3], header1.mass[4], header1.mass[5]);
	printf("time: %e .\n", header1.time);
	printf("redshift: %e .\n", header1.redshift);
	printf("flag_sfr: %d .\n", header1.flag_sfr);
	printf("flag_feedback: %d .\n", header1.flag_feedback);
	printf("npartTotal[6]: %d %d %d %d %d %d .\n", header1.npartTotal[0], header1.npartTotal[1], header1.npartTotal[2], header1.npartTotal[3], header1.npartTotal[4], header1.npartTotal[5]);
	printf("flag_cooling: %d .\n", header1.flag_cooling);
	printf("num_files: %d .\n", header1.num_files);
	printf("BoxSize: %e .\n", header1.BoxSize);
	printf("Omega0: %e .\n", header1.Omega0);
	printf("OmegaLambda: %e .\n", header1.OmegaLambda);
	printf("HubbleParam: %e .\n", header1.HubbleParam);
	// Additional stuff added into header by JMK:
	printf("Dark Energy w0: %e .\n" , header1.w0);
	printf("Dark Energy wa: %e .\n" , header1.wa);
	printf("Comoving Distance from Observer (in kpc/h): %e .\n" , header1.comoving_distance);
	
	fflush(stdout);
}

// Print particle position buffer on screen:
void show_particle_position_buffer(float *particle_position_buffer, int number_of_particles)
{
	int i;
	for (i=0; i<number_of_particles; i++) printf("Position coordinates of Particle %d: %e %e %e\n", i, particle_position_buffer[3*i], particle_position_buffer[3*i+1], particle_position_buffer[3*i+2]);
	fflush(stdout);
}

