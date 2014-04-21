/*                                                                                                                                                                                       
 *  mpi_support.c
 *  Inspector Gadget                                                                                                                                                                           
 *                                                                                                                                                                                       
 *  Created by Jan Michael Kratochvil at University of Miami on 11/10/09.                                                                                                                                                
 *  Copyright 2009. All rights reserved.                                                                                                                               
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
#include "allocation.h"
#include "2D-plane_multi.h"
#include "mpi_support.h"

#include "ini.h"

#if defined(MPI_COMPILE)
#include <mpi.h>
#endif

#define MAXNAME 200
#define MAXPATHNAME 1000


void prepare_MPI(int argc, char **argv) // Constructs parameters for MPI run (some of these parameters are otherwise constructed by the ig_condor_writer script for Condor and passed as arguments to the serial runs. argc and argv not really needed here (unlike in the Condor equivalent).
{

  FILE *input_file2;
  // char comoving_file[1000];
  char filename[1000], sampler_filename[1000];
  int number_of_planes, source_counter;
  int r;
  number_of_planes=0;


  #if defined(MPI_COMPILE)
  MPI_Comm_rank(MPI_COMM_WORLD, &superrank);
  MPI_Comm_size(MPI_COMM_WORLD, &supersize);

  if(superrank==0){

    //Prompt user for correct number of arguments
    if(argc<5){
      fprintf(stderr,"Mode 1 usage: %s <num_ics_to_process(N)> <num_snapshots_per_ic> <ini_parameter_file> <ic_identifier_1> ... <ic_identifier_N>\n",*argv);
      fprintf(stderr,"Mode 2 usage: %s <num_cosmologies(N)> <num_processors_per_cosmology> <ini_parameter_file> <cosmo_id_1> ... <cosmo_id_N>\n",*argv);
      MPI_Finalize();
      exit(1);
    }

    //Read in arguments:
    printf("Called with number of arguments (count: 1 means no additional argument): %d\n", argc);
    //printf("Argument content: %s\n", argv[0]);
    //printf("Argument 1 content: %s\n", argv[1]);
    fflush(stdout);

    printf("--------------------------------------------------\n");
    printf("Inspector Gadget, Version X-5.0\n");
    printf("--------------------------------------------------\n");
    fflush(stdout);

    printf("Read-in completed.\n");
    fflush(stdout);
    // Read-in completed.
  
  }

  parameters.number_of_cosmologies=atoi(argv[1]); // How many cosmological models in parallel.
  parameters.number_of_processes=atoi(argv[2]); // On how many processors does each simulation run.

  FILE *input_file1;
  // input_file1 = fopen ("IG_parameters_P.txt", "r");
  input_file1 = fopen (argv[3], "r"); // Third input argument of program is parameter filename. 
  
  //read_analysis_parameter_file(input_file1, 1);
  
  //<AP>
  //Parse INI parameter file using INI parsing library
  if(ini_parse_file(input_file1,handler,&parameters)<0){
    fprintf(stderr,"There was a problem reading the input parameter file\n");
    MPI_Finalize();
    exit(1);
  }
  //</AP>

  fclose(input_file1);
  
  //<AP>
  if(superrank==0){
    printf("Done reading main parameter file; here are the options that you gave me:\n");
    fflush(stdout);
    print_options(stdout,&parameters);
    fflush(stdout);
  }

  //</AP>

  int subprocess_in_cosmology;
  int number_of_subprocesses; // number of processes into which one plane or map generation is split (before Version X-4.0 there were no such subprocesses).
  int number_of_purposes;

  // These two values are hardcoded and optimized for the specific machine you're running on, in this case NYBlue Blue Gene at Brookhaven National Laboratory: 
  if (parameters.mode==1)
    {
      number_of_subprocesses=1; // This number is not used in this version. //64 // This number should best equal the number of files into which each Gadget-2 snapshot is split.
      number_of_purposes=2;
    }
  else if (parameters.mode==2)
    {
        number_of_subprocesses=1; // must be set to 1. Does not work anymore with other settings, and is essentially obsolete, because this functionality can be performed better with OpenMP during ray tracing, which has been incorporated.
    }
  else
    {
      printf("ERROR: Mode not supported.\n");
      exit(1);
    }

  if (parameters.mode==1)
    {
      
        parameters.cosmology_number = superrank / (parameters.number_of_processes);
        parameters.process_number = superrank % (parameters.number_of_processes);
        MPI_Comm_split(MPI_COMM_WORLD, parameters.cosmology_number, parameters.process_number, &sim_comm);
        parameters.ThisTask=0; // ThisTask enumerates subprocesses (in Mode 2), but no subprocesses exist in Mode 1 (i.e. there is only one such process).
        
        

    }
  else if (parameters.mode==2)
    {
      parameters.cosmology_number = superrank / (parameters.number_of_processes*number_of_subprocesses);
      subprocess_in_cosmology = superrank % (parameters.number_of_processes*number_of_subprocesses);

      // Split depending on cosmologies:
      MPI_Comm_split(MPI_COMM_WORLD, parameters.cosmology_number, subprocess_in_cosmology, &sim_comm);

      // parameters.process_number counts processes as one per map and within a cosmology (as if there were no subprocesses).
      // parameters.ThisTask shall count which task within a map, i.e. subprocesses within a map process (one map split by MPI).
      // parameters.NTasks=number_of_subprocesses.

      parameters.process_number = subprocess_in_cosmology / number_of_subprocesses;
      parameters.ThisTask = subprocess_in_cosmology % number_of_subprocesses;

      // Split multiple processes for each plane (mode=1) or map (mode=2):
      MPI_Comm_split(sim_comm, parameters.process_number, parameters.ThisTask, &map_comm);

    }

  if (parameters.cosmology_number>=parameters.number_of_cosmologies || parameters.process_number>=parameters.number_of_processes) // If superfluous process on a CPU, end the process right here.
    {
      printf("Process (superrank %d process_number %d map task %d) is superfluous. Being terminated without effect on result.\n", superrank, parameters.process_number, parameters.ThisTask);

      MPI_Finalize();
      exit(0);
    }


  strcpy(parameters.simulation_codename, argv[4+parameters.cosmology_number]);

  printf( "MPI subprocess %d of prosess %d in cosmology %d.\n", parameters.ThisTask, parameters.process_number, parameters.cosmology_number);
  

  #else

    printf("ERROR: Non-MPI running of Inspector Gadget not supported since Version X-4.0 (and possibly earlier ones). Aborting.\n");
    exit(1)
  #endif

  printf("Running Mode %d\n", parameters.mode);

    
    parameters.number_of_source_planes=0; // Having dedicated source planes does not work anymore since Inspector Gadget 5.0, in favor of being able to output at the end of any desired plane during the run (facilitates producing more redshifts, although it disables the ability to choose source redshifts exactly.


  //////////////////////////////////////
  // Derived parameters (automatic):
  //////////////////////////////////////

  // Paths:
  char temp[2000];
  sprintf(temp, "%s", parameters.snapshot_path);
  sprintf(parameters.snapshot_path, "%s/%s", temp, parameters.simulation_codename);
  sprintf(temp, "%s", parameters.plane_output_path);
  sprintf(parameters.plane_output_path, "%s/%s/%s", temp, parameters.simulation_codename, parameters.planes_folder);
  sprintf(parameters.plane_urpath, "%s", parameters.plane_path);
  sprintf(temp, "%s", parameters.plane_path);
  sprintf(parameters.plane_path, "%s/%s/%s", temp, parameters.simulation_codename, parameters.planes_folder);
  sprintf(temp, "%s", parameters.map_output_path);
  sprintf(parameters.map_output_path, "%s/%s/%s", temp, parameters.simulation_codename, parameters.maps_folder);

  if (parameters.mode==1)
    {
      parameters.generate_planes=1;
      parameters.generate_maps=0;

    }
  else if (parameters.mode==2)
    {
      parameters.generate_maps=1;
      parameters.generate_planes=0;
    }
  else
    {
      printf("ERROR: Selected Mode does not exist.\n");
#if defined(MPI_COMPILE)
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
      exit(1);
    }

  // Read in first_ and last_snapshot are stored directly in global_first_ and _last_snapshot variables (so only need to construct first_ and last_snapshot variables here).

  if (parameters.mode==1)
    {

        parameters.parallel=parameters.number_of_processes;
        number_of_planes=(parameters.global_last_snapshot-parameters.global_first_snapshot)/parameters.snapskip+1;
      
        int planes_per_process;
        // without source planes:
        // planes_per_process=number_of_planes/parameters.number_of_processes;
        // if (number_of_planes%parameters.number_of_processes!=0) planes_per_process++;
        // with source planes (one per process):
        planes_per_process=number_of_planes/(parameters.number_of_processes-parameters.number_of_source_planes);
        if (number_of_planes%(parameters.number_of_processes-parameters.number_of_source_planes)!=0) planes_per_process++;
        
        int number_of_regular_processes; // (as opposed to terminal plane processes)
        number_of_regular_processes=number_of_planes/planes_per_process;
        if (number_of_planes%(parameters.number_of_processes-parameters.number_of_source_planes)!=0) number_of_regular_processes++;
        
        
        // Worknote: FIX REGULAR PROCESSES FOR TERMINAL PLANES RUN TOGETHER WITH REGULAR.
        
        if (parameters.process_number*planes_per_process < number_of_planes) // here do regular planes:
        {
            parameters.source_redshift=-1;
            parameters.plane_before_source=-1;
            parameters.source_comoving_distance=-1;
            parameters.last_snapshot=plane_to_snapshot(parameters.process_number*planes_per_process, parameters.global_last_snapshot, parameters.snapskip); // note that planes are numbered in reverse order compared to snapshots, because the plane closest to the observer is Plane 0 which corresponds to the very last global snapshot.
            parameters.first_snapshot=plane_to_snapshot((parameters.process_number+1)*planes_per_process-1, parameters.global_last_snapshot, parameters.snapskip);
            if (parameters.first_snapshot<0) parameters.first_snapshot=0;
        }
        else // here do terminal (source) planes (one per parameters.process_number - this functionality has currently been lost):
        {
        
            printf("superrank %d is superfluous process (would be terminal plane, which currently doesn't work in Inspector Gadget 5.0.\n", superrank);
            MPI_Finalize();
            exit(0);
	}
        
      parameters.survey_angle=0.0; // survey angle is not needed for this mode (plane generation).

    }
  else if (parameters.mode==2)
    {
      int local_number_of_realizations;

      parameters.first_snapshot=parameters.global_first_snapshot;
      parameters.last_snapshot=parameters.global_last_snapshot;

      parameters.global_first_realization=parameters.first_realization;
      parameters.global_last_realization=parameters.last_realization;
      local_number_of_realizations=(int) ceil((double) (parameters.last_realization-parameters.first_realization+1)/((double) parameters.number_of_processes));
      parameters.parallel=(int) ceil((double) (parameters.last_realization-parameters.first_realization+1)/((double) local_number_of_realizations));
      if (superrank<4) printf("Special: parameters.parallel=%d\n", parameters.parallel);
      parameters.first_realization+=(parameters.process_number*local_number_of_realizations);
      parameters.last_realization=parameters.first_realization+local_number_of_realizations-1;
      if (parameters.last_realization>parameters.global_last_realization) parameters.last_realization=parameters.global_last_realization;

    }
    else
    {
        printf("ERROR: Selected mode does not exist. Only Modes 1 and 2 are valid.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }

    printf("MPI Param (Process %d): (global) first and last snapshot: (%d %d) %d %d, snapskip, number_of_planes: %d %d, parallel, number_of_source_planes: %d %d\n", parameters.process_number, parameters.global_first_snapshot, parameters.global_last_snapshot, parameters.first_snapshot, parameters.last_snapshot, parameters.snapskip, number_of_planes, parameters.parallel, parameters.number_of_source_planes);
    fflush(stdout);

    parameters.number_of_snapshots=parameters.last_snapshot-parameters.first_snapshot+1;
    parameters.snapshots=parameters.global_last_snapshot-parameters.global_first_snapshot+1;
    parameters.number_of_realizations=parameters.last_realization-parameters.first_realization+1;
    parameters.realization=parameters.first_realization;

    parameters.number_of_boxcenters=parameters.global_last_snapshot+3; // large enough number, more than last snapshot used (guarantees that enough, since there are no negatively numbered snapshots).

      // MPI-specific checks:

  if (parameters.parallel>parameters.number_of_processes)
    {
      printf("ERROR: Insufficient number of parallel MPI processes for requested run.\n Number of Processes: %d, Number of Processes needed: %d.\n", parameters.NTasks, parameters.parallel);
      fflush(stdout);
      exit(1);
    }


  if (parameters.process_number>=parameters.parallel)
    {
      printf("Process %d is superfluous. Being terminated without effect on result.\n", parameters.process_number);
      #if defined(MPI_COMPILE)
      MPI_Finalize();
      #endif
      exit(0);
      // ends superfluous processes that have nothing to do.
    }

  // Note:
  // parameters.number_of_cosmologies = how many cosmologies in parallel.
  // parameters.cosmology_number = the how many-th cosmology (count starting at zero)
  // parameters.number_of_processes = number of processes allocated for each cosmology (could be automatized later).
  // parameters.process_number = rank among processes of one cosmology.
  // parameters.parallel = number of processes pro cosmology actually needed (must be <= parameters.number_of_processes).


  printf("Done with mpi_support parameter construction.\n");
  fflush(stdout);


}



void cosmology_sampler_initialization(void)
{

  char sampler_filename[2000];
  int r, i, j;


  if (parameters.plane_shift==0)
    {
      // Load cosmology sampler file with ic and realization mixing numbers:                                                                                                               
      // This option uses a random number file which lists ic number and realization number consecutively for each lens plane (i.e. schematically: ic_plane0 realization_plane0 ic_plane1 realization_plane1 ic_plane2 realization_plane2, etc.).                                                                                                                                    
      cosmology_sampler=(int **)malloc(parameters.global_last_realization*sizeof(int *));
      assert(cosmology_sampler!=NULL);
      for (r=0; r<parameters.global_last_realization; r++)
	{
          cosmology_sampler[r]=(int *)malloc(2*(parameters.plane_before_source+1)*sizeof(int));
          assert(cosmology_sampler[r]!=NULL);
	}

      if (parameters.cosmology_number+1<=parameters.fiducial)
	{

          sprintf(sampler_filename, "%s", parameters.plane_randomizer_file_fiducial);
	}
      else
	{
          sprintf(sampler_filename, "%s", parameters.plane_randomizer_file_general);
	}
      printf("Reading in random number file (for plane selection):\n%s\n", sampler_filename);
      fflush(stdout);
      read_sampler(sampler_filename, cosmology_sampler, parameters.global_last_realization, 2*(parameters.plane_before_source+1), parameters.max_realizations, parameters.galaxy_subfield-1); // reads up to last realization lines.        

    
        for (i=0; i<parameters.global_last_realization; i++)
        {
            for (j=0; j<parameters.plane_before_source+1; j++)
            {
                if (cosmology_sampler[i][2*j]==0)
                {
                    printf("ERROR: Initial Condition (IC) should never be zero. realization, plane, ic_number, realization_number: %d %d %d %d\n", i, j, cosmology_sampler[i][j*2], cosmology_sampler[i][j*2+1]);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                    exit(1);
                }
            }
        }
    }
  //////////////////////////////////////////           
  else
    {
        // Load cosmology sampler file with ic and realization mixing numbers:
        // This option uses a different random number file, having 5 consecutive numbers per plane instead of just two (as above), the last three being (in order) one int to determine rotations and mirrorings of lens planes, and two values for shifts in horizontal and vertical directions upon ray tracing. Use this option if using a minimal number of lens planes (i.e. just those which mutually exclusively slice the snapshot cube instead of a separately sliced piece for each map realization).
        cosmology_sampler=(int **)malloc(parameters.global_last_realization*sizeof(int *));
        assert(cosmology_sampler!=NULL);
        for (r=0; r<parameters.global_last_realization; r++)
        {
            cosmology_sampler[r]=(int *)malloc(5*(parameters.plane_before_source+1)*sizeof(int));
            assert(cosmology_sampler[r]!=NULL);
        }

        if (parameters.cosmology_number+1<=parameters.fiducial)
        {
            sprintf(sampler_filename, "%s", parameters.plane_randomizer_file_fiducial);
        }
        else
        {
            sprintf(sampler_filename, "%s", parameters.plane_randomizer_file_general);
        }
        printf("Reading in random number file (for plane selection and shifting):\n%s\n", sampler_filename);
        fflush(stdout);
        read_sampler(sampler_filename, cosmology_sampler, parameters.global_last_realization, 5*(parameters.plane_before_source+1), parameters.max_realizations, parameters.galaxy_subfield-1); // reads up to last realization lines. subfield count starts at 1, but at zero in read_sampler routine, that's why called with subfield-1.
   
        for (i=0; i<parameters.global_last_realization; i++)
        {
            for (j=0; j<parameters.plane_before_source+1; j++)
            {

                if (cosmology_sampler[i][5*j]==0)
                {
                    printf("ERROR: Initial Condition (IC) should never be zero. realization, plane, ic_number, realization_number, mirrot, shiftx, shifty: %d %d %d %d %d %d %d \n", i, j, cosmology_sampler[i][j*5], cosmology_sampler[i][j*5+1], cosmology_sampler[i][j*5+2], cosmology_sampler[i][j*5+3], cosmology_sampler[i][j*5+4]);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                    exit(1);
                }
            }
        }

    } // end of plane_shift else clause.  


}





void diagnostic(void)
{
    FILE *diagnos_file;
    char dasfile[20];
    
    sprintf(dasfile, "./diagnos_%d.txt", parameters.process_number);
  diagnos_file=fopen(dasfile, "a");
  parameters.diagnostic++;
  fprintf(diagnos_file, "Diagnostic Writeout %d:\n %e %e %e %e %d %d %d %d\n\n",
	  parameters.diagnostic,
	  parameters.Omega_m,
	  parameters.Omega_Lambda,
	  parameters.source_redshift,
	  parameters.source_comoving_distance,
	  parameters.number_of_source_planes,
	  parameters.number_of_planes,
	  parameters.plane_before_source);
  fclose(diagnos_file);
  
}


