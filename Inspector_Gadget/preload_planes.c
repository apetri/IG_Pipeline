/*
 *  preload_planes.c
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at the University of KwaZulu-Natal on 03/26/2014.
 *  Copyright 2014. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <mpi.h>

#include "main.h"
// The above is needed only for the file read path global variables.
#include "2D-plane_multi.h"
// The above is needed by fits.h.
#include "fits.h"
// The above is needed to read potential planes FITS files.
#include "preload_planes.h"


/////////////////////////////////
// Note 1: process_number in this file is always within a cosmology, i.e. the MPI rank within sim_comm.
// Note 2: plane_number starts count at zero, while plane_realization and plane_sim_ic start counting at 1 in the argument of all functions in this file (but within the function this may be changed, until the output of the function starts at 1 again). This is the same convention as in the rest of Inspector Gadget.


// Do one-time read-in (preloading) of all potential planes, and put them into an RMA window on different MPI processes for one-sided MPI_Get calls later during code run:
double* read_in_all_planes(MPI_Comm sim_comm, MPI_Win *plane_storage_window, int number_of_planes, int number_of_plane_realizations, int number_of_sim_ics, int first_sim_ic, int nx, int ny, int convergence_direct)
{
    int process_number, number_of_processes, plane_instance, number_of_plane_instances, number_of_plane_instances_per_process, i;
    MPI_Aint plane_storage_length;
    int plane_number, plane_realization, plane_sim_ic;
    int plane_exists;
    MPI_Aint displacement_unit;

    MPI_Comm_rank(sim_comm, &process_number);
    MPI_Comm_size(sim_comm, &number_of_processes);

    // Determine numebr of plane instances needed for each process for the whole job to be able to hold all potential planes simultaneously, preloaded in RMA windows:
	number_of_plane_instances=number_of_planes*number_of_plane_realizations*number_of_sim_ics;
	number_of_plane_instances_per_process=number_of_plane_instances/number_of_processes;
	if (number_of_plane_instances%number_of_processes!=0) number_of_plane_instances_per_process++; // in this case, some processes will eventually contain one fewer than number_of_instances preloaded potential planes, while others will have the full complement of number_of_instances.

    // Allocate plane_storage array, which is the double array covered by the RMA window plane_storage_window and will be available for one-sided MPI communications during the run (in particular for MPI_Get):
	plane_storage_length=number_of_plane_instances_per_process*nx*ny; // number of doubles the plane_storage array has to span to be able to hold all preloaded potential planes.
	double *plane_storage=(double *)malloc(plane_storage_length*sizeof(double));
    assert(plane_storage!=NULL);
    for (i=0; i<plane_storage_length; i++) plane_storage[i]=0.0; // initialize plane_storage array with zeros.


	// Now load all required planes:

	for (plane_instance=0; plane_instance<number_of_plane_instances_per_process; plane_instance++)
	{
        plane_exists=plane_instance_to_plane_parameters(plane_instance, number_of_planes, number_of_plane_realizations, number_of_sim_ics, first_sim_ic, sim_comm, &plane_number, &plane_realization, &plane_sim_ic);
	if (plane_exists==1) preload_potential_plane(plane_storage, plane_instance, plane_number, plane_realization, plane_sim_ic, nx, ny, convergence_direct); // load potential plane only if it actually exists, otherwise skip this function call, so that this part of the plane_storage array stays initialized to zero.
	} 


        // Now create an RMA window on the plane_storage:
        displacement_unit=sizeof(double);

	printf("Superrank %d, process_number %d (before Window Create): Parameters to create window: number_of_planes: %d, number_of_plane_realizations: %d, number_of_sim_ics: %d, number_of_plane_instances: %d, number_of_plane_instances_per_process: %d, number_of_processes: %d.\n", superrank, process_number, number_of_planes, number_of_plane_realizations, number_of_sim_ics, number_of_plane_instances, number_of_plane_instances_per_process, number_of_processes);
	printf("Superrank %d, process_number %d (befroe Window Create): plane_storage_length: %ld (in units of maps: %d), displacement_unit: %ld, window size in bytes: %ld.\n", superrank, process_number, plane_storage_length, plane_storage_length/(nx*ny), displacement_unit, plane_storage_length*displacement_unit);
	fflush(stdout);
	

	MPI_Win_create(plane_storage, plane_storage_length*displacement_unit, displacement_unit, MPI_INFO_NULL, sim_comm, plane_storage_window);
    
    
	printf("Superrank %d, process_number %d (Created Window): Parameters to create window: number_of_planes: %d, number_of_plane_realizations: %d, number_of_sim_ics: %d, number_of_plane_instances: %d, number_of_plane_instances_per_process: %d, number_of_processes: %d.\n", superrank, process_number, number_of_planes, number_of_plane_realizations, number_of_sim_ics, number_of_plane_instances, number_of_plane_instances_per_process, number_of_processes);
    printf("Superrank %d, process_number %d (Created Window): plane_storage_length: %ld (in units of maps: %d), displacement_unit: %ld, window size in bytes: %ld.\n", superrank, process_number, plane_storage_length, plane_storage_length/(nx*ny), displacement_unit, plane_storage_length*displacement_unit);
    fflush(stdout);


    // MPI_Barrier(sim_comm);
    // MPI_Win_free(plane_storage_window);
    // MPI_Finalize();
    // exit(0);

    
    return plane_storage;

}


// This is the file read-in function during one-time preloading of potential planes:
void preload_potential_plane(double *plane_storage, int plane_instance,  int plane_number, int plane_realization, int plane_sim_ic, int nx, int ny, int convergence_direct)
{
    char filename[2000];
    MPI_Aint target_displacement;
    
	target_displacement=plane_instance*nx*ny; // plane_instance counts displacement from window beginning in numbers of potential planes, while target_displacement counts the same but in number of doubles (so need to multiply by number of pixels in potential plane).
    
    // If compute convergence maps directly from density planes as lens planes instead of via gravitational potential planes and transverse spatial derivatives (this is mostly for testing and code verification of the FFTW, etc., as this mode cannot produce shear maps, only convergence maps):
    if (convergence_direct!=0) sprintf(filename, "%s/%s_ic%d/%s/%s_%s_ic%d_%04dxy_%04dr_%03dp.%s", parameters.plane_urpath, parameters.simulation_codename, plane_sim_ic, parameters.planes_folder, parameters.density_basename, parameters.simulation_codename, plane_sim_ic, nx, plane_realization, plane_number, parameters.extension);
    // Regular case using gravitational potential planes as lens planes:
    else sprintf(filename, "%s/%s_ic%d/%s/%s_%s_ic%d_%04dxy_%04dr_%03dp.%s", parameters.plane_urpath, parameters.simulation_codename, plane_sim_ic, parameters.planes_folder, parameters.potential_basename, parameters.simulation_codename, plane_sim_ic, nx, plane_realization, plane_number, parameters.extension);
    
    // Read 2D image in FITS file and store it in the RMA plane preloading window with proper displacement from window beginning for one-sided MPI communications later:
	readFITSpotential_singleplane(filename, plane_storage+target_displacement);
	

}


// Based on plane_instance, this funtion returns the correct parameters to read in the proper plane (most parameters start count at 1, except for plane_number which starts at zero):
int plane_instance_to_plane_parameters(int plane_instance, int number_of_planes, int number_of_plane_realizations, int number_of_sim_ics, int first_sim_ic, MPI_Comm sim_comm, int *plane_number, int *plane_realization, int *plane_sim_ic)
{
    int process_number, number_of_processes, total_plane_nr, temp;
    
    MPI_Comm_rank(sim_comm, &process_number);
    MPI_Comm_size(sim_comm, &number_of_processes);
    
    total_plane_nr=number_of_processes*plane_instance+process_number;
    // total plane number is all planes computed consecutively in an array as follows (C-ordering, last index runs first, then the others): [plane_sim_ic][plane_number][plane_realization]. Planes are distributed one plane per process in sim_comm and then looping around after each process got one.
    *plane_sim_ic=total_plane_nr/(number_of_planes*number_of_plane_realizations);
    temp=total_plane_nr%(number_of_planes*number_of_plane_realizations);
    *plane_number=temp/number_of_plane_realizations;
    *plane_realization=temp%number_of_plane_realizations;
    
    // Make sure proper values are only returned for planes which actually exist (prevents nonsensical numbers to be read from memory, should such a request ever come):
    if (*plane_sim_ic<number_of_sim_ics && *plane_number<number_of_planes && *plane_realization<number_of_plane_realizations)
    {
        *plane_realization+=1;
        *plane_sim_ic+=first_sim_ic;
        return 1; // tells calling function that a valid plane was identified for preloading.
    }
    else // Set plane parameters to nonsensical negative values if no such plane exists in reality:
    {
        *plane_sim_ic=-1; *plane_number=-1; *plane_realization=-1; // nonsensical values.
        return 0; // tells calling function that there is no such plane and that no potential plane shall be preloaded (i.e. the subsequent call to the loading function is to be skipped).
    }
}


// This function gets the rank and displacement within window of the target process to perform a one-sided MPI communication (to be called before MPI_Get()):
void plane_parameters_to_displacement(int plane_number, int plane_realization, int plane_sim_ic, int number_of_planes, int number_of_plane_realizations, int first_sim_ic, int nx, int ny, MPI_Comm sim_comm, int *target_rank, MPI_Aint *target_displacement)
{
    int number_of_processes, plane_instance, total_plane_nr;
    MPI_Comm_size(sim_comm, &number_of_processes);
    
    // Adjust for local values in this function (which are the values in the RMA window):
    plane_realization-=1; // count starts at zero in RMA window, while plane_realization count always starts at 1 (unlike sim ICs which can start at any number for a given run with 1 being the lowest possible).
    plane_sim_ic-=first_sim_ic; // first used sim IC (does not have to be ic1) is in zero spot without target_displacement in plane_storage array.

    total_plane_nr=plane_sim_ic*number_of_planes*number_of_plane_realizations+plane_number*number_of_plane_realizations+plane_realization; // total plane number is all planes computed consecutively in an array as follows (C-ordering, last index runs first, then the others): [plane_sim_ic][plane_number][plane_realization]. Planes are distributed one plane per process in sim_comm and then looping around after each process got one.
    
	*target_rank=total_plane_nr%number_of_processes; // which process to do an RMA MPI_Get from.
	plane_instance=total_plane_nr/number_of_processes; // number of planes to skip in RMA window during the MPI_Get.
	*target_displacement=plane_instance*nx*ny; // plane_instance counts displacement from window beginning in numbers of potential planes, while target_displacement counts the same but in number of doubles (so need to multiply by number of pixels in potential plane).
}



