/*
 *  main.c
 *  Inspector Gadget, Version 5.3
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 5/8/2007.
 *  Later improved at the University of Miami and the University of KwaZulu-Natal.
 *  Last updated on April 4, 2014.
 *  Copyright 2007-14 Jan Michael Kratochvil. All rights reserved.
 *  This is not a public release.
 *
 */




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <time.h>

#include <complex.h>
#include <fftw3.h>

#include "main.h"
#include "io_routines.h"
#include "mathematics.h"
#include "allocation.h"
#include "2D-plane_multi.h"
#include "fits.h"
#include "preload_planes.h"
#include "weak_lensing_multi.h"
#include "mpi_support.h"
#include "endianness.h"

// Header files below optional?  (Introduced here with comoving distance splining, but probably not necessary in this file.)
#include "darkenergy.h"
#include "comoving_distance.h"
#include "initialization.h"




#if defined(MPI_COMPILE)
#include <mpi.h>
MPI_Comm sim_comm;
MPI_Comm snapshot_comm; // USED??
MPI_Comm map_comm;
#endif
int superrank;
int supersize;


double UnitMass_in_g;


struct io_header_1 header1;

struct snapshot_control *Snapshot;


struct analysis_parameters parameters;

FILE *input_file_ran; // global file pointer for random number file with snapshot box rotations during potential plane generation. Needs to be global in this implementation, so that it continuously advances lines on each use, and does not reset when the funtion is called (this could be programmed better).

int feedback;

int **cosmology_sampler;


// Inspector Gadget, main program:

int main(int argc, char **argv)
{

  

	int i, j;
	int start_seed;
	double original_survey_angle;
    int realization;
    int plane_number;
    
    
    // For preloading potential planes:
    int preload_planes;
    double *plane_storage; // array containing preloaded planes.
    MPI_Win plane_storage_window; // MPI RMA Window on the above array, allowing one-sided MPI RMA communication.
    int number_of_plane_realizations, number_of_sim_ics, first_sim_ic;

    // Copy from global parameters (note that the values of all these are set in function prepare_MPI() via global variables):
    int ray_tracing, plane_assignment_mode, raypoint_averaging, convergence_direct, galaxy_catalogue_type;
    int nx, ny, NbinsX, NbinsY, number_of_planes; // not sure these local versions are safe to use, as the global ones get updated during reading of the header of potential planes; so currently using these only where safe, in parallel to global ones.
    // Further variables to be localized: survey_angle, survey_angle_in_rad.
    
    
    
	printf("--------------------------------------------------\n");
	printf("Inspector Gadget, Version X-5.0\n");
	printf("--------------------------------------------------\n");
	fflush(stdout);

    // Start timing Inspector Gadget run:
	// parameters.start_time=clock();
	parameters.start_time=time(NULL);

	//Read in arguments:
	printf("Called with number of arguments (count: 1 means no additional argument): %d\n", argc);
	//printf("Argument content: %s\n", argv[0]);
	//printf("Argument 1 content: %s\n", argv[1]);
	fflush(stdout);


	// Set extension of map file type:
#ifdef HAVE_FITS
        sprintf(parameters.extension, "fit"); // FITS file                                        
#else
    exit(1); // Alternatives currently not implemented, must use FITS format (cfitsio library).
#endif



#if defined(MPI_COMPILE)
	MPI_Init( &argc, &argv );
	prepare_MPI(argc, argv);
#else
    
    printf("ERROR: Inspector Gadget does not run anymore without MPI. Must have MPI and compile with an MPI compiler like mpicc to run this program. Aborting.\n");
    printf("Note: You do not have to have OpenMP necessarily; that is optional. But MPI must be present.\n");
    fflush(stdout);
    exit(1);
#endif

        printf("Read-in completed.\n");
	fflush(stdout);
	// Read-in completed.

	// Initial settings and hardcoded parameters:
	parameters.byteswap=0;
	parameters.endianness_set=0;
	parameters.darkenergy_initialized=0;
    parameters.chi_initialized=0;
	parameters.snapshot_allocated=0;
	parameters.diagnostic=0;
	parameters.plane_assignment_mode=0; // this puts particles on planes with periodic completion on boundaries (all other modes are obsolete, since always full plane files produces, may one day change this back of partial box rotations allowed).

    
    // New parameter copying to eliminate global variables:
    preload_planes=parameters.preload_planes;
    convergence_direct=parameters.convergence_direct;
    ray_tracing=parameters.ray_tracing;
    plane_assignment_mode=parameters.plane_assignment_mode;
    raypoint_averaging=parameters.raypoint_averaging;
    galaxy_catalogue_type=parameters.galaxy_catalogue_type;
    nx=parameters.nx;
    ny=parameters.ny;
    NbinsX=parameters.NbinsX;
    NbinsY=parameters.NbinsY;
    // Used in ray-tracing with preloading of planes:
    number_of_plane_realizations=parameters.number_of_plane_realizations;
    number_of_sim_ics=parameters.number_of_sim_ics;
    first_sim_ic=parameters.first_sim_ic;


    
    // Calculating derived quantities from read in parameters (may get overwritten later by code automatically, e.g. after reading in planes):
	parameters.survey_angle_in_rad=(parameters.survey_angle/360.0)*2.0*M_PI;
	parameters.nxny = parameters.nx * parameters.ny;
    
    
    

	if (parameters.generate_planes==1 && parameters.generate_maps==0)
	{
        parameters.number_of_planes=(parameters.global_last_snapshot-parameters.global_first_snapshot)/parameters.snapskip+1;
		if (parameters.plane_before_source >= parameters.number_of_planes) 
		{
			printf("ERROR: plane_before_source=%d while number_of_planes only =%d.\n", parameters.plane_before_source, parameters.number_of_planes);
			parameters.number_of_planes=parameters.plane_before_source+1; 
			printf("Some source planes (and thus the planes before sources) are detached from contiguous plane stack calculated from observer. Ray-tracing will give wrong results if these far out lying last plane before source are used as terminal planes, because some planes in between will be missing.\n");
			fflush(stdout);
			exit(1);
		} 
	}
	else
	  {
	    parameters.number_of_planes=parameters.plane_before_source+1; // plus one for zero plane. This is sufficient, will utilize exactly as many planes as needed to reach source plane, and snapshot specifics don't need to be specified properly in the Condor Script (although it is recommended to always leave them there at the values used for plane generation). Will need to change this line, once adopt strategy to generate maps for multiple redshifts simultaneously.
	    printf("Number of planes set automatically to plane_before_source+1, and is: %d.\n", parameters.number_of_planes);
	  }


    // Further global-to-local variable transfers:
    number_of_planes=parameters.number_of_planes;
    
    
    
	feedback=parameters.feedback;
	
	parameters.seed+=parameters.process_number; // add process number to default seed.
	

	// Outputting on screen to double check values:
	printf("---------------------------\n");
	if (argc>1)
	{
		printf("This node runs with (parameters passed by script):\n");
		printf("Simulation Codename / Modelfolder (for output filenames and folders): %s \n", parameters.simulation_codename);
		printf("Generate planes (yes if =1): %d \n", parameters.generate_planes);
		printf("Generate maps (yes if =1): %d \n", parameters.generate_maps);
		printf("Process Number: %d \n", parameters.process_number);
		printf("Total number of parallel processes: %d \n", parameters.parallel);
		printf("File number of first snapshot processed by this process: %d \n", parameters.first_snapshot);
		printf("File number of last snapshot processed by this process: %d \n", parameters.last_snapshot);
		printf("File number of global first snapshot: %d \n", parameters.global_first_snapshot);
		printf("File number of global last snapshot: %d \n", parameters.global_last_snapshot);
		printf("File number of first realization processed by this process: %d \n", parameters.first_realization);
		printf("File number of last realization processed by this process: %d \n", parameters.last_realization);

		printf("Snapskip: %d\n", parameters.snapskip);
		printf("Source Redshift: %e\n", parameters.source_redshift);
		printf("Number of last FITS plane before source plane: %d\n", parameters.plane_before_source);
		printf("Source Comoving Distance: %e\n", parameters.source_comoving_distance);
		printf("Simulation Model Name (tag in header of FITS file): %s\n", parameters.modelname);


		printf("Snapshot basename: %s \n", parameters.snapshot_name);
		printf("Potential basename: %s \n", parameters.potential_basename);
		printf("Convergence: %s \n", parameters.convergence_basename);
		printf("Shear 1: %s \n", parameters.shear1_basename);
		printf("Shear 2: %s \n", parameters.shear2_basename);
	

		if (parameters.first_snapshot>parameters.last_snapshot || parameters.first_snapshot>parameters.global_last_snapshot)
		  {
		    printf("First snapshot larger than last snapshot: %d %d. Superfluous node. Exiting this node.\n", parameters.first_snapshot, parameters.last_snapshot);
		    exit(0);
		  }

		if (parameters.first_realization>parameters.last_realization)
		  {
		    printf("First realization > last realization: %d %d. Superfluous node. Exiting this node.\n", parameters.first_realization, parameters.last_realization);
		    exit(0);
		  }

	}

	printf("---------------------------\n");
	fflush(stdout);
	
	printf("These paths are used to read and write for run of code:\n");
	printf("Snapshot Path (N-body simulation snapshot boxes are read in from this location): %s\n", parameters.snapshot_path);
	printf("Parameter Path (Parameter files are read in from this location): %s\n", parameters.parameter_path);
	printf("Plane Path (Plane FITS files are read in from to this location for ray-tracing): %s\n", parameters.plane_path);
	printf("Plane Output Path (Plane Output FITS files are written to this location): %s\n", parameters.plane_output_path);
	printf("Map Output Path (Map Output FITS files are written to this location): %s\n", parameters.map_output_path);
				
	printf("---------------------------\n");

	printf("These parameter values are used for run of code:\n");
	printf("Particle Species plotted: %d\n", parameters.species);
	printf("NbinsX (number of bins for weak lensing maps): %d\n", parameters.NbinsX);
	printf("NbinsY (number of bins for weak lensing maps): %d\n", parameters.NbinsY);
	printf("nx (number of grid points in 2D lens planes): %d\n", parameters.nx);
	printf("ny (number of grid points in 2D lens planes): %d\n", parameters.ny);
	printf("Survey Angle (if zero, later set automatically to maximum allowed by N-body simulation): %e\n", parameters.survey_angle);
	printf("Plane Assignment Mode (Transverse Snapshot Box size if zero, fitting (padded) survey angle otherwise): %d\n", parameters.plane_assignment_mode);
	printf("Scramble Mode (Snapshot Boxes Orientation): %d\n", parameters.scramble_mode);
	printf("Ray Tracing (on if !=0): %d\n", parameters.ray_tracing);
	printf("(Global) Number of Snapshots: %d\n", parameters.snapshots);
	printf("Number of planes: %d\n", parameters.number_of_planes);
	printf("Feedback Level (0 is lowest): %d\n", parameters.feedback);
	printf("Generate planes? (If 0, already made planes are read in): %d \n", parameters.generate_planes);
	printf("Generate weak lensing maps (form existing or just generated planes)? (Yes if =1): %d \n", parameters.generate_maps);
	printf("Plane shift (additional randomization of plane center during WL map generation) (yes if =1): %d\n", parameters.plane_shift);
	printf("Initial Random Number Seed (for plane shift): %d\n", parameters.seed);
	printf("Cell Embedding (1: cloud-in-cell, 2: TSC (preferred)): %d\n", parameters.cell_embedding);
	printf("Raypoint Averaging (how potential derivatives are calculated at point where lightray pierces lens plane (1: linear (from NR), 2: TSC): %d\n", parameters.raypoint_averaging);
	printf("How many parallel processes expected (each to be executed manually or by IG Master Script): %d\n", parameters.parallel);
	printf("Process Number (among parallel processes): %d\n", parameters.process_number);
	printf("Convergence direct (do regular weak lensing ray-tracing only if =0) : %d\n", parameters.convergence_direct);
	printf("---------------------------\n");
	fflush(stdout);

	

	// Allocating snapshot_control:
       if(!(Snapshot=malloc(parameters.number_of_boxcenters*sizeof(struct snapshot_control))))
       {
           fprintf(stderr,"failed to allocate memory.\n");
           exit(0);
       } 

	
       printf("Number of planes before allocating Plane array: %d\n", parameters.number_of_planes);

	// Allocating Plane pointer array:
	Plane= (struct plane_2D *) malloc(parameters.number_of_planes * sizeof(struct plane_2D));
	assert(Plane != NULL);
	for (plane_number=0;plane_number<parameters.number_of_planes;plane_number++)
	{
		Plane[plane_number].particles_written= (double *) malloc((parameters.seed_block+1) * sizeof(double));
		assert(Plane[plane_number].particles_written != NULL);

		for (j=0; j<3; j++)
		  {
		    Plane[plane_number].rRot[j]= (double *) malloc((parameters.seed_block+1) * sizeof(double));
		    assert(Plane[plane_number].rRot[j] != NULL);
		    Plane[plane_number].rShift[j]=(double *) malloc((parameters.seed_block+1) * sizeof(double));
		    assert(Plane[plane_number].rShift[j] != NULL);
		  }
	}
	// Plane count starts at Plane 0, that's the first plane away from (and closest to) observer.
    
	// Allocating FITSheader pointer array:
	FITSheader= (struct fitsheader *) malloc(parameters.number_of_planes * sizeof(struct fitsheader));
	assert(FITSheader != NULL);
	

        // fastforward in random number file:
        char ran_line[1000];
        int ff, ffstop, ran_line_length;
        ff=0;
	// Fast forward making space for 1000 realizations:
	if (parameters.plane_before_source<0) ffstop=1000*parameters.process_number+(parameters.first_realization-1); // in the regular case when plane in the middle of the plane stack.
	else ffstop=1000*parameters.plane_before_source+(parameters.first_realization-1); // in the case of final plane before source plane (ensures same random numbers (snapshot box orientation) as for the same plane when calculated as a regular plane for source planes at higher redshifts.)

    // New IG 5.0:
    if (parameters.mode==1) ffstop=0;
    /*
    while ((ff < ffstop ) && ((ran_line_length=fgetline(input_file_ran, ran_line, sizeof(ran_line))) > 0))
        {
	        ff++;
        }
        printf("fastforward check: %d lines skipped in random number file.\n", ff);
     */
    
    
	timing();

        ensure_random();


		// Remember parameter values that have special meaning for map creation cycle (to set automatic mode):
		original_survey_angle=parameters.survey_angle;

		start_seed=parameters.seed;
		

		if (parameters.generate_planes!=0)
		{
		
		  for (plane_number=0; plane_number<parameters.number_of_planes; plane_number++)
			  {
			    match_planes(Plane, plane_number, parameters.global_last_snapshot, parameters.snapskip);
			    Plane[plane_number].set==0; // sets flag that plane header parameters other than Plane[i].snapshot (which was set by match_planes(i)) have not been set properly yet and setup_plane() will need to be called later.

			    // Set range which planes get processes by this process, based on which snapshots were uploaded for this process:
			    if (Plane[plane_number].snapshot==parameters.last_snapshot) parameters.first_plane=plane_number;
			    if (Plane[plane_number].snapshot==parameters.first_snapshot) parameters.last_plane=plane_number;
			  }

		
			for (plane_number=parameters.first_plane; plane_number<=parameters.last_plane; plane_number++) // process numbers start at 1.
			{
                input_file_ran = fopen(parameters.snapshot_rotation_randomizer_file, "r"); // open file with random numbers.

			    for (realization=parameters.first_realization; realization<=parameters.last_realization; realization+=parameters.seed_block) compute_plane(Plane, plane_number, parameters.nx, parameters.ny, parameters.seed_block, realization, parameters.last_realization, parameters.number_of_planes, parameters.source_redshift, parameters.plane_before_source, parameters.scramble_mode, parameters.species, parameters.cell_embedding, parameters.plane_assignment_mode);
                fclose(input_file_ran);
			}
			
			timing();
			printf("Done generating planes.\n");

		}
		
		parameters.seed=start_seed; // resets random number seed to original value in parameterfile (plus process_number addition), such that same ray-tracing is run with the same random number seed regardless if plane generation is done in the same program run or separately beforehand.
		
		if (parameters.generate_maps!=0)
		{

            // Determine number (and first) of plane realizations and N-body sim initial conditions by quickly going through random number file; then broadcast to all processes:
            /*
            if (superrank==0) cosmology_sampler_scan(&number_of_plane_realizations, &number_of_sim_ics, &first_plane_realization, &first_sim_ic);
            MPI_Bcast(&number_of_plane_realizations, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&number_of_sim_ics, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&first_plane_realization, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&first_sim_ic, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            */
            // The above is for a future code version with automatic determination of number of plane realizations and N-body ICs.
            
            if (preload_planes!=0)
            {

	        printf("Superrank %d starting to preload potential planes...\n", superrank);
	        fflush(stdout);

                plane_storage=read_in_all_planes(sim_comm, &plane_storage_window, number_of_planes, number_of_plane_realizations, number_of_sim_ics, first_sim_ic, nx, ny, convergence_direct);

		MPI_Barrier(sim_comm); // wait for each process to complete preloading planes before continuing.
		printf("Superrank %d successfully passed through MPI_Barrier after preloading of planes completed for all processes in the sim_comm communicator.\n", superrank);
		fflush(stdout);
            }

		  for (parameters.galaxy_subfield=parameters.first_galaxy_subfield; parameters.galaxy_subfield<=parameters.last_galaxy_subfield; parameters.galaxy_subfield++)
		    {
			
		      // Initialize Cosmology Sampler (for N-body initial conditions, plane realizations; for plane rotation, mirroring and shifts, etc.)
		      cosmology_sampler_initialization();

			for (realization=parameters.first_realization; realization<=parameters.last_realization; realization++)
			{
                parameters.survey_angle=original_survey_angle;
				    readpotentialplanes_header(FITSheader, Plane, parameters.number_of_planes, parameters.plane_shift, parameters.convergence_direct, nx, realization, parameters.source_redshift, &(parameters.source_comoving_distance)); // read headers of all potential planes into memory.
                // Note: the above line updates many of the parameters variables from the potential plane files, so don't use local versions for insertion below. However, the function below does not modify these.
                weak_lensing2(parameters.NbinsX, parameters.NbinsY, parameters.nx, parameters.ny, ray_tracing, plane_assignment_mode, raypoint_averaging, Plane, convergence_direct, galaxy_catalogue_type, realization, parameters.survey_angle_in_rad, parameters.source_comoving_distance, parameters.number_of_planes, parameters.plane_before_source, parameters.source_redshift, parameters.plane_shift, preload_planes, number_of_plane_realizations, first_sim_ic, sim_comm, &plane_storage_window); // perform ray-tracing over these potential planes (individual planes are read in here as well).

			} // end loop over parameters.realization.
			
			timing();
			printf("Superrank %d: Done generating WL quantities for Galaxy Subfield %d.\n",superrank, parameters.galaxy_subfield);


			// Free cosmology sampler:
			printf("Superrank %d now freeing cosmology sampler...", superrank);
			fflush(stdout);
			    int r;
			    for (r=0; r<parameters.global_last_realization; r++)
			      {
				free(cosmology_sampler[r]);
			      }
			    free(cosmology_sampler);

			printf("\nSuperrank %d: Got past cosmology_sampler freeing.\n", superrank);
			fflush(stdout);

                    } // end loop over parameters.galaxy_subfield.

			timing();
                        printf("Superrank %d: Done generating WL maps.\n", superrank);
                        fflush(stdout);

		} // end of if clause parameters.generate_maps.

		/////////////////////////////////////////////////////////////


		
		parameters.seed=start_seed; // resets randon number seed to default value again, for any future applications that may be implemented here (like plane generation and map generation above).
		
		///// Implement further applications here (on the fundamental level of creating planes or maps): /////////////
								
		//////////// End of original multi section. ///////////////////////
	

    if(parameters.generate_planes!=0) free_all_snapshots(); // Free memory allocated for particle arrays and snapshot_control before ending program.


    if (parameters.generate_maps!=0 && preload_planes!=0)
    {
        MPI_Win_free(&plane_storage_window);
        free(plane_storage);
        // Both of the above were allocated by read_in_all_planes().
    }

    fflush(stdout);
    
#if defined(MPI_COMPILE)
    MPI_Finalize();
#endif

    if (superrank==0) printf("Finished lens plane generation or ray tracing successfully. Inspector Gadget exited normally.\n");
    fflush(stdout);
    exit(0);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Functions taken from Gadget-2 analysis example codes and modified:
// These include mostly loading a multi-file snapshot, (de)allocating memory for snapshots, as well as unit conversions (the latter are not used in the code).
// (Other input/output functions can be found in io_routines.c, and allocation function in allocation.c.)



/* this routine allocates the memory for the 
 * particle data.
 */
void allocate_snapshot_memory(int NumPart, int snapshot_number)
{
  if (feedback>2) printf("Allocating memory for Snapshot %d...\n", snapshot_number);

  int bytes;

  if(!(Snapshot[snapshot_number].P=malloc(bytes = NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"Process %d (superrank) failed to allocate memory (%g MB instead of %g MB required).\n", superrank, bytes/(1024.0*1024.0),0);
      printf("FATAL ERROR: Process %d (superrank) failed to allocate memory (%g MB instead of %g MB required).\n", superrank, bytes/(1024.0*1024.0), 0);
#if defined(MPI_COMPILE)
      MPI_Finalize();
#endif
      exit(1);
    }
  
  Snapshot[snapshot_number].P--;   /* start with offset 1 */  // NEVER WRITE TO P[0] OR YOU MAY DELETE SOME OTHER VARIABLE!!
  
  
  if(!(Snapshot[snapshot_number].Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"Process %d (superrank) failed to allocate memory.\n", superrank);
      printf("FATAL ERROR: Process %d (superrank) failed to allocate memory.\n", superrank);
#if defined(MPI_COMPILE)
      MPI_Finalize();
#endif
      exit(1);
    }
  
  Snapshot[snapshot_number].Id--;   /* start with offset 1 */
  

  
	// Deallocate all of the above by calling free_snapshot_arrays(snapshot_number) etc., or collectively by free_snapshot(snapshot_number).

  if (feedback>2) printf("Done allocating memory for Snapshot %d.\n", snapshot_number);
}

void free_snapshot_arrays(int snapshot_number) // frees all large particle arrays in a snapshot
{
	Snapshot[snapshot_number].P++;
	free(Snapshot[snapshot_number].P);
	Snapshot[snapshot_number].P=NULL;
	Snapshot[snapshot_number].Id++;
	free(Snapshot[snapshot_number].Id);
	Snapshot[snapshot_number].Id=NULL;
}

void free_all_snapshots(void) // frees all snapshots (only information retained is in planes (struct Plane).
{
	int snapshot_number;
	printf("Superrank %d: Freeing allocated memory:\n Arrays of Snapshots:", superrank);
	fflush(stdout);
	for (snapshot_number=0; snapshot_number<parameters.snapshots; snapshot_number++)
	  {
	    printf("Superrank %d: Will now free snapshot_number %d.\n", superrank, snapshot_number);
	    fflush(stdout);
	if (Snapshot[snapshot_number].P!=NULL)
	{
		free_snapshot_arrays(snapshot_number); // frees large particle arrays only if they haven't already been freed (checking only first, assuming all collectively deleted or not, with "free_snapshot_arrays()".
		printf(" (sr %d)-->%d ", superrank, snapshot_number);
		fflush(stdout);
	}
	  }
	free(Snapshot); // frees snapshot_control, losing all information about snapshot boxes. All relevant information has to be reduced to 2D planes by now and stored in the Plane arrays.
	printf(", snapshot control array. Memory free (superrank %d).\n", superrank);
	fflush(stdout);
}




/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
void reordering(int snapshot_number) // Function not actually used in code.
{
  int i,j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i=1; i<=Snapshot[snapshot_number].NumPart; i++)
    {
      if(Snapshot[snapshot_number].Id[i] != i)
	{
	  psource= Snapshot[snapshot_number].P[i];
	  idsource=Snapshot[snapshot_number].Id[i];
	  dest=Snapshot[snapshot_number].Id[i];

	  do
	    {
	      psave= Snapshot[snapshot_number].P[dest];
	      idsave=Snapshot[snapshot_number].Id[dest];

	      Snapshot[snapshot_number].P[dest]= psource;
	      Snapshot[snapshot_number].Id[dest]= idsource;
	      
	      if(dest == i) 
		break;

	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Snapshot[snapshot_number].Id++;   // correcting for original offset Id-- after allocation;
  free(Snapshot[snapshot_number].Id);

  printf("space for particle ID freed\n");
}


  
  
/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.)
 */

// Uses header1 (variable) and Snapshot (pointer), both for reading (this function also allocates part of Snapshot).
void load_snapshot_multi(char *fname, int file_number, int snapshot_number, int plane_number)
{
  FILE *fd;
  char   buf[200];
  int    i,k,dummy,ntot_withmasses;
  int    n,pc,pc_new,pc_sph;

	// now local variables, to be assigned to snapshot_control structure:
  	int NumPart;
	int Ngas;
	int Nhalo;
	int files;
	files=1;  // Needs to be set to 1 (old rudiment, file number is counted by file_number now).
	int NumPartThisFile;
	int NgasThisFile;
	int NhaloThisFile;
	
	float vel_temp[3]; // velocity placeholder for read-in (overwritten by every particle read, because velocities not used)

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  
  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      //if(files>1)
	sprintf(buf,"%s.%d",fname,file_number);
      //else sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("ERROR: Process %d (superrank) can't open file `%s`\n", superrank, buf);
	  exit(0);
	}

      printf("Process %d (superrank) reading `%s' ...\n", superrank, buf); fflush(stdout);

	  endianness.byteswap=parameters.byteswap; // set byteswap flag to previously determined relative endianness of Gadget-2 snapshot.
	  // endianness.byteswap is the flag executed on reading of every file. parameters.endianness is the flag containing the relative endianness of the Gadget-2 snapshot (which could be different from other files read during the run).

      fread(&dummy, sizeof(dummy), 1, fd);
      header_fread(&header1, sizeof(header1), 1, fd);
      // Endianness replacement (old): fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

        {
            for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++) NumPart+= header1.npartTotal[k];
            Ngas= header1.npartTotal[0];
            Nhalo= header1.npartTotal[1];
	  
            for(k=0, NumPartThisFile=0, ntot_withmasses=0; k<5; k++) NumPartThisFile+= header1.npart[k];
            NgasThisFile= header1.npart[0];
            NhaloThisFile= header1.npart[1];
        }

        for(k=0, ntot_withmasses=0; k<5; k++)
        {
            if(header1.mass[k]==0)
                ntot_withmasses+= header1.npart[k];
        }

        printf("Process %d (superrank) about to allocate NumPartThisFile=%d for Snapshot Particle Array.\n", superrank, NumPartThisFile);
        fflush(stdout);

    
        allocate_snapshot_memory(NumPartThisFile, snapshot_number);
      
        printf("Process %d (superrank) has just allocated NumPartThisFile=%d for Snapshot Particle Array.\n", superrank, NumPartThisFile);
        fflush(stdout);
 
        parameters.snapshot_allocated=1;

        SKIP;
        for(k=0,pc_new=pc;k<6;k++)
        {
            for(n=0;n<header1.npart[k];n++)
            {
            
                float_fread(&Snapshot[snapshot_number].P[pc_new].Pos[0], sizeof(float), 3, fd);
                // Endianness replacement (old): fread(&Snapshot[snapshot_number].P[pc_new].Pos[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;
	
        printf("Process %d (superrank) finished reading the-positions.\n", superrank);
        fflush(stdout);

        SKIP;
        for(k=0,pc_new=pc;k<6;k++)
        {
            for(n=0;n<header1.npart[k];n++)
            {
                float_fread(&vel_temp[0], sizeof(float), 3, fd); // read into placeholder variable, won't need velocities in code.
                // Correct line for actual read in: float_fread(&Snapshot[snapshot_number].P[pc_new].Vel[0], sizeof(float), 3, fd);
                // Endianness replacement (old): fread(&Snapshot[snapshot_number].P[pc_new].Vel[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;
    
        printf("Process %d (superrank) finished reading the-velocities.\n", superrank);
        fflush(stdout);

        SKIP;
        for(k=0,pc_new=pc;k<6;k++)
        {
            for(n=0;n<header1.npart[k];n++)
            {
                int_fread(&Snapshot[snapshot_number].Id[pc_new], sizeof(int), 1, fd);
                // Endianness replacement (old): fread(&Snapshot[snapshot_number].Id[pc_new], sizeof(int), 1, fd);
                if (n==0)  printf("Process %d (superrank) read first ID: %d.\n", superrank, Snapshot[snapshot_number].Id[pc_new]);
                pc_new++;
            }
        }
        SKIP;

        printf("Process %d (superrank) finished reading the-IDs.\n", superrank);
      fflush(stdout);

        if(ntot_withmasses>0)
            SKIP;
        for(k=0, pc_new=pc; k<6; k++)
        {
            for(n=0;n<header1.npart[k];n++)
            {
                Snapshot[snapshot_number].P[pc_new].Type=k;

                if(header1.mass[k]==0)
                    float_fread(&Snapshot[snapshot_number].P[pc_new].Mass, sizeof(float), 1, fd);
                // Endianness replacement (old): fread(&Snapshot[snapshot_number].P[pc_new].Mass, sizeof(float), 1, fd);
                else
                    Snapshot[snapshot_number].P[pc_new].Mass= header1.mass[k];
                pc_new++;
            }
        }
        if(ntot_withmasses>0)
            SKIP;
      
        printf("Process %d (superrank) finished reading the-masses.\n", superrank);
        fflush(stdout);


      /* Skip this part as long as do not have any gas particles:
       *
       *

      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      float_fread(&Snapshot[snapshot_number].P[pc_sph].U, sizeof(float), 1, fd);
	      // Endianness replacement (old): fread(&Snapshot[snapshot_number].P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      float_fread(&Snapshot[snapshot_number].P[pc_sph].Rho, sizeof(float), 1, fd);
	      // Endianness replacement (old): fread(&Snapshot[snapshot_number].P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  float_fread(&Snapshot[snapshot_number].P[pc_sph].Ne, sizeof(float), 1, fd);
		  // Endianness replacement (old): fread(&Snapshot[snapshot_number].P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		Snapshot[snapshot_number].P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
	}

      printf("Process %d (superrank) finished reading the-gas.\n", superrank);
      fflush(stdout);

      *
      *
      * End of skipped gas part
      */

      fclose(fd);
    }

    printf("Process %d (superrank) closed the-file.\n", superrank);
    fflush(stdout);

	// Before the end of the function, pass on header values to permanent ones in snapshot_control array:

	Snapshot[snapshot_number].NumPart=NumPartThisFile;
	Snapshot[snapshot_number].Ngas=NgasThisFile;
	Snapshot[snapshot_number].Nhalo=NhaloThisFile;
	Snapshot[snapshot_number].Time= header1.time;
	Snapshot[snapshot_number].Redshift= header1.redshift;
	
	for (i=0;i<6;i++)
	{
		Snapshot[snapshot_number].NumPartTotal[i]=header1.npartTotal[i];
		Snapshot[snapshot_number].mass[i]=header1.mass[i];
	}
	Snapshot[snapshot_number].files=header1.num_files; // number of files from which N-body snapshot is composed.
	Snapshot[snapshot_number].Boxsize=header1.BoxSize; // in kpc/h comoving.
	Snapshot[snapshot_number].Omega0=header1.Omega0; // Total matter fraction in universe (CDM, baryons, etc.)
	Snapshot[snapshot_number].OmegaLambda=header1.OmegaLambda; // Dark energy fraction.
	Snapshot[snapshot_number].H_0=header1.HubbleParam*100.0; // Hubble constant today (H_0) of the N-body simulation run by Gadget-2. Conversion factor 100 is because Gadget-2 gives h instead of H_0.
	
	Snapshot[snapshot_number].w0=header1.w0;
	Snapshot[snapshot_number].wa=header1.wa;
	Snapshot[snapshot_number].ns=0.0; // These will be incorporated in the Gadget-2 snapshot header in the next generation of simulations.
	Snapshot[snapshot_number].sigma_8=0.0;
	Snapshot[snapshot_number].initial_condition=0.0; // This one needs to be implemented from simulation codename generation routine.
	Snapshot[snapshot_number].comoving_distance=header1.comoving_distance;

	Snapshot[snapshot_number].set=1; // Flag set, shows snapshot control array Snapshot has been properly initialized for this snapshot number.

    //////////////////////////////
    // Manual Adjustment Section:
    // These will be incorporated in the Gadget-2 snapshot header in the next generation of simulations.
    Snapshot[snapshot_number].ns=0.96;
	Snapshot[snapshot_number].sigma_8=0.8; // MANUAL: change here as appropriate.
    
    //////////////////////////////
    
    
    
	printf("Process %d (superrank) finished read-in function.\n", superrank);
	fflush(stdout);
}



void load_snapshot_multi_header(char *fname, int file_number, int snapshot_number)
{
  FILE *fd;
  char   buf[200];
  int    i,k,dummy,ntot_withmasses;
  int    n,pc,pc_new,pc_sph;

  // now local variables, to be assigned to snapshot_control structure:                   
  int NumPart;
  int Ngas;
  int Nhalo;
  int files;
  files=1;  // Needs to be set to 1 (old rudiment, file number is counted by file_number now).

  int NumPartThisFile;
  int NgasThisFile;
  int NhaloThisFile;

  struct io_header_1 header;


#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
                                                                             
  sprintf(buf,"%s.%d",fname,file_number);


  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      //if(files>1)                                                                             
      sprintf(buf,"%s.%d",fname,file_number);
      //else sprintf(buf,"%s",fname);                                                           

      if(!(fd=fopen(buf,"r")))
        {
          printf("can't open file `%s`\n",buf);
          exit(0);
        }


  printf("reading `%s' ...\n",buf); fflush(stdout);

  endianness.byteswap=parameters.byteswap; // set byteswap flag to previously determined relative endianness of Gadget-2 snapshot.                                                     
  // endianness.byteswap is the flag executed on reading of every file. parameters.endianness is the flag containing the relative endianness of the Gadget-2 snapshot (which could be different from other files read during the run).                                                

  fread(&dummy, sizeof(dummy), 1, fd);
  header_fread(&header, sizeof(header), 1, fd);
  // Endianness replacement (old): fread(&header1, sizeof(header1), 1, fd);                 
  fread(&dummy, sizeof(dummy), 1, fd);


  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
    NumPart+= header1.npartTotal[k];
  Ngas= header1.npartTotal[0];
  Nhalo= header1.npartTotal[1];

  for(k=0, NumPartThisFile=0, ntot_withmasses=0; k<5; k++)
    NumPartThisFile+= header1.npart[k];
  NgasThisFile= header1.npart[0];
  NhaloThisFile= header1.npart[1];

  for(k=0, ntot_withmasses=0; k<5; k++)
    {
      if(header1.mass[k]==0)
	ntot_withmasses+= header1.npart[k];
    }

  printf("Process %d (superrank) read in Snapshot header NumPartThisFile=%d\n", superrank, NumPartThisFile);
  fflush(stdout);

  fclose(fd);

    }

  // Before the end of the function, pass on header values to permanent ones in snapshot_control array:                                                                                  

  Snapshot[snapshot_number].NumPart=NumPartThisFile;
  Snapshot[snapshot_number].Ngas=NgasThisFile;
  Snapshot[snapshot_number].Nhalo=NhaloThisFile;
  Snapshot[snapshot_number].Time= header.time;
  Snapshot[snapshot_number].Redshift= header.redshift;

  for (i=0;i<6;i++)
    {
      Snapshot[snapshot_number].NumPartTotal[i]=header1.npartTotal[i];
      Snapshot[snapshot_number].mass[i]=header1.mass[i];
    }
  Snapshot[snapshot_number].files=header.num_files; // number of files from which N-body snapshot is composed.                                                                          
  Snapshot[snapshot_number].Boxsize=header.BoxSize; // in kpc/h comoving.                
  Snapshot[snapshot_number].Omega0=header.Omega0; // Total matter fraction in universe (CDM, baryons, etc.)                                                                             
  Snapshot[snapshot_number].OmegaLambda=header.OmegaLambda; // Dark energy fraction.     
  Snapshot[snapshot_number].H_0=header.HubbleParam*100.0; // Hubble constant today (H_0) of the N-body simulation run by Gadget-2. Conversion factor 100 is because Gadget-2 gives h instead of H_0.                                                                                   

  Snapshot[snapshot_number].w0=header.w0;
  Snapshot[snapshot_number].wa=header.wa;
  Snapshot[snapshot_number].ns=0.0; // These will be incorporated in the Gadget-2 snapshot header in the next generation of simulations.                                                 
  Snapshot[snapshot_number].sigma_8=0.0;
  Snapshot[snapshot_number].initial_condition=0.0; // This one needs to be implemented from simulation codename generation routine.                                                      
  Snapshot[snapshot_number].comoving_distance=header.comoving_distance;

  Snapshot[snapshot_number].set=1; // Flag set, showing that snapshot control array has been properly initialized for this snapshot number.

  return;

}

