/*
 *  2D-plane_multi.c
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 3/19/08.
 *  Copyright 2008. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <time.h>

#include <complex.h>
#include <fftw3.h>

#include "main.h"
#include "mathematics.h"
#include "allocation.h"
#include "io_routines.h"
#include "2D-plane_multi.h"
#include "weak_lensing_multi.h"
#include "endianness.h"
#include "darkenergy.h"
#include "comoving_distance.h"

#ifdef HAVE_FITS
#include "fits.h"
#endif


// From 2D-plane.c:
struct plane_2D *Plane;
struct fitsheader *FITSheader;
double particles_written;
double deflection_angle[2];



void match_planes(struct plane_2D *Plane, int plane_number, int global_last_snapshot, int snapskip)
{

  int snapshot_number;
 
  // Plane data initialization:
  // Specifies which snapshot (in absolute counts, filename of snapshot) belongs to plane, and what its comoving distance from observer is (this is the only information provided in a separate file, and calculated separately from Gadget-2; all other info about the snapshot can be read from the snapshot header directly later):
  snapshot_number=plane_to_snapshot(plane_number, global_last_snapshot, snapskip);
  Plane[plane_number].snapshot=snapshot_number;

  if (feedback>2) printf("Plane_number, comoving distance: %d %e.\n", plane_number, Plane[plane_number].comoving_distance);

}


int plane_to_snapshot(int plane_number, int global_last_snapshot, int snapskip)
{
  int snapshot_number;
  snapshot_number=global_last_snapshot-(plane_number*snapskip);
  return snapshot_number;
}


double cloud_in_cell(double x_particle, double y_particle, double* density_array, int nx, int ny, double binwidth_x, double binwidth_y, int plane_assignment_mode)
     // assigns an N-body simulation particle to density bin on 2D plane, to create density grid for Poissons equation, which will be solved with FFTW
{

    double x_diff, y_diff;
    double x_bin, y_bin, x_frac, y_frac;
    double sector_A, sector_B, sector_C, sector_D, particle_fraction;
    int Xbin_A, Ybin_A, main_bin;
    int x_shift, y_shift, y_shift_sign, x_period, y_period;

    int nxny;
    
    nxny=nx*ny;


    // Sector distribution in bins:
    //   C  |D
    //--------
    //      |
    //   A  |B
    //(main)|

    // Determine main bin cell, into which particle will be inserted:
    Xbin_A= ((int) floor(x_particle/binwidth_x+0.5*nx+0.5)); // snapshot number refers here to the snapshot box, in which the plane is located, based on its redshift.
    Ybin_A= ((int) floor(y_particle/binwidth_y+0.5*ny+0.5));


    if (plane_assignment_mode==0)
    {
        // In this mode, compensate for boxshifting (crosshairs) by completing periodically:
        if (Xbin_A<0)
        {
            Xbin_A=Xbin_A+nx;
            x_particle=x_particle+nx*binwidth_x;
        }
        if (Xbin_A>=nx)
        {
            Xbin_A=Xbin_A-nx;
            x_particle=x_particle-nx*binwidth_x;
        }
        if (Ybin_A<0)
        {
            Ybin_A=Ybin_A+ny;
            y_particle=y_particle-ny*binwidth_y;
        }
        if (Ybin_A>=ny)
        {
            Ybin_A=Ybin_A-ny;
            y_particle=y_particle+ny*binwidth_y;
        }

        printf("Xbin_A, Ybin_A: %d %d\n", Xbin_A, Ybin_A);

        // If particle still out of range, something is wrong and need to abort to prevent array boundary violation:
        if (Xbin_A<0 || Xbin_A>=nx || Ybin_A<0 || Ybin_A>=ny) // checks that particle gets properly assigned and there is no violation of array boundaries.
        {
            printf("ERROR: Particle out of range of density array. Unable to cloud-in-cell particle.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(1); // Terminate program, as this is a fatal error that needs fixing, since in this plane_assignment_mode, all particles should be assignable to a plane.
        }

    }
    else // In this plane_assignment_mode, checking if particle is way out of grid (more than one bin distance away), so no fraction can touch:
    {
        if (Xbin_A<-1 || Xbin_A>nx || Ybin_A<-1 || Ybin_A>ny)
        {
            return 0;
        }
    }

    // Determine the other cells, in which sectors will be deposited:
    x_bin=(Xbin_A-0.5*nx)*binwidth_x;
    y_bin=-(Ybin_A-0.5*ny)*binwidth_y;

    x_diff=x_particle-x_bin; // Difference in physical coordinates between particle and bin center.
    y_diff=y_particle-y_bin;

    // Calculate sector distribution:
    x_frac=fabs(x_diff/binwidth_x);
    y_frac=fabs(y_diff/binwidth_y);
    sector_D=x_frac*y_frac;
    sector_B=x_frac-sector_D;
    sector_C=y_frac-sector_D;
    sector_A=1.0-(sector_B+sector_C+sector_D);

    if (sector_A<0 || sector_B<0 || sector_C<0 || sector_D<0) printf("x_particle %e y_particle %e, x_bin %e y_bin %e, x_frac %e y_frac %e, Sectors A, B, C, D: %e %e %e %e\n", x_particle, y_particle, x_bin, y_bin, x_frac, y_frac, sector_A, sector_B, sector_C, sector_D);

    if (plane_assignment_mode==0) // more complicated scheme doing periodic completion at boundaries:
    {
        x_period=0;
        y_period=0;

        if (Xbin_A==0 && x_diff<=0) x_period=nx;
        if (Xbin_A==(nx-1) && x_diff>0 ) x_period=-nx;
        if (Ybin_A==0 && y_diff>0) y_period=nxny;
        if (Ybin_A==(nx-1) && y_diff<=0) y_period=-nxny;
        if (x_diff>0) x_shift=1;
        else x_shift=-1;
        if (y_diff>0) y_shift=-nx;
        else y_shift=nx;

        // Fill distribution fractions (above sectors) in corresponding bins:
        main_bin=Ybin_A*nx+Xbin_A;
#ifdef HAVE_OpenMP
#pragma omp critical (CICupdate)
        {
#endif
        density_array[main_bin]=density_array[main_bin]+sector_A; // complex number addition, sector_A is real, but that's fine.
        density_array[main_bin+x_shift+x_period]=density_array[main_bin+x_shift+x_period]+sector_B;
        density_array[main_bin+y_shift+y_period]=density_array[main_bin+y_shift+y_period]+sector_C;
        density_array[main_bin+x_shift+y_shift+x_period+y_period]=density_array[main_bin+x_shift+y_shift+x_period+y_period]+sector_D;
#ifdef HAVE_OpenMP
        }
#endif
        
        return 1.0;
    }

    else // only particle fractions lying within the lattice actually assigned. particle_written taking into account fractionally written particles:
    {
        if (x_diff>0) x_shift=1;
        else x_shift=-1;
        if (y_diff>0)
        {
            y_shift=-nx;
            y_shift_sign=-1; // needed for quick checking below, because abs(y_shift)!=1, unlike x_shift above.
        }
        else
        {
            y_shift=nx;
            y_shift_sign=1;
        }

        // Fill distribution fractions (above sectors) in corresponding bins and return proper collective fraction to counter:
        main_bin=Ybin_A*nx+Xbin_A;
        
            
        if (Xbin_A>0 && Xbin_A<(nx-1) && Ybin_A>0 && Ybin_A<(ny-1))
        {
            
#ifdef HAVE_OpenMP
#pragma omp critical (CICupdate)
            {
#endif
	    particle_fraction = 0.0;
            // If in center of array with at least one cell row/column from edge, assign all four components without further questioning:
            density_array[main_bin]=density_array[main_bin]+sector_A; // complex number addition, sector_A is real, but that's fine.
            particle_fraction=particle_fraction+sector_B;
            density_array[main_bin+y_shift]=density_array[main_bin+y_shift]+sector_C;
            density_array[main_bin+x_shift+y_shift]=density_array[main_bin+x_shift+y_shift]+sector_D;
#ifdef HAVE_OpenMP
            }
#endif

            return 1; // whole particle assigned.
        }
        else
        {
            particle_fraction=0;

#ifdef HAVE_OpenMP
#pragma omp critical (CICupdate)
            {
#endif
            if (Xbin_A>=0 && Xbin_A<nx && Ybin_A>=0 && Ybin_A<ny)
            {
                density_array[main_bin]=density_array[main_bin]+sector_A; // complex number addition, sector_A is real, but that's fine.
                particle_fraction=particle_fraction+sector_A;
            }
            if (Xbin_A+x_shift>=0 && Xbin_A+x_shift<nx && Ybin_A>0 && Ybin_A<ny)
            {
                density_array[main_bin+x_shift]=density_array[main_bin+x_shift]+sector_B;
                particle_fraction=particle_fraction+sector_B;
            }
            if (Xbin_A>=0 && Xbin_A<nx && Ybin_A+y_shift_sign>0 && Ybin_A+y_shift_sign<ny)
            {
                density_array[main_bin+y_shift]=density_array[main_bin+y_shift]+sector_C;
                particle_fraction=particle_fraction+sector_C;
            }
            if ((Xbin_A+x_shift)>=0 && (Xbin_A+x_shift)<nx && (Ybin_A+y_shift_sign)>0 && (Ybin_A+y_shift_sign)<ny)
            {
                density_array[main_bin+x_shift+y_shift]=density_array[main_bin+x_shift+y_shift]+sector_D;
                particle_fraction=particle_fraction+sector_D;
            }
#ifdef HAVE_OpenMP
            }
#endif
            
     
            return particle_fraction; // possibly only a fraction of the particle actually assigned to plane grid.
        }
            


    }

}


////////////////////////////////////////////////////////////////
// TSC:
////////////////////


// Plane[plane_number].binwidth_x
double TSC_assign(double x_particle, double y_particle, double* density_array, int nx, int ny, double binwidth_x, double binwidth_y)
{
    int i, j;
    double gridpoint_x, gridpoint_y;
    int Xbin, Ybin, Xbin_A, Ybin_A, points_check, main_bin, adjust;

    double particle_fraction, total_particle_fraction;

    double weightcorrection=1000.0; // corrects to megaparsec units, since that's what Premadi's W_TSC takes (inferred from paper, although not explicitly stated there).
    double hh;
    
    int nxny;
    
    nxny=nx*ny;
    
    hh=binwidth_x*binwidth_y/(weightcorrection*weightcorrection); // divides TSC weight function out properly, so that inserted particle fraction is 1.
    

    // Determine main bin cell, into which particle will be inserted:
    Xbin_A= ((int) floor(x_particle/binwidth_x+0.5*nx+0.5)); // snapshot number refers here to the snapshot box, in which the plane is located, based on its redshift.
    Ybin_A= ((int) floor(y_particle/binwidth_y+0.5*ny+0.5));

    points_check=2;
    main_bin=Ybin_A*nx+Xbin_A;

    double test_weightwidth;
    test_weightwidth=1.5*weightcorrection*(1.0/binwidth_x);


    particle_fraction=0.0;
    total_particle_fraction=0.0;

    for (i=-points_check;i<=points_check;i++)
    {
        for (j=-points_check;j<=points_check;j++)
        {
            Xbin=Xbin_A+j;
            Ybin=Ybin_A+i;

            // Determine comoving coordinate of gridpoint, without periodic completion! This serves just as a distance measure for the weight function:
            gridpoint_x=((Xbin)-0.5*nx-0.5)*binwidth_x;
            gridpoint_y=((Ybin)-0.5*ny-0.5)*binwidth_y;

            adjust=i*nx+j;

            // Adjust for periodicity of boundary conditions for actual indices in array:
            if (Xbin<0)
            {
                adjust+=nx;
            }
            else if (Xbin>=nx)
            {
                adjust-=nx;
            }
            if (Ybin<0)
            {
                adjust+=nxny;
            }
            else if (Ybin>=ny)
            {
                adjust-=nxny;
            }

            // Insert fraction of particle into that gridpoint:
            particle_fraction=W_TSC(fabs((x_particle-gridpoint_x))/binwidth_x)*W_TSC(fabs((y_particle-gridpoint_y)/binwidth_y));
#ifdef HAVE_OpenMP
#pragma omp critical (TSCupdate)
{
#endif

            density_array[main_bin+adjust]=density_array[main_bin+adjust]+particle_fraction;
#ifdef HAVE_OpenMP
}
#endif
            total_particle_fraction+=particle_fraction;

        }
    }

    if (total_particle_fraction<0.99999 || total_particle_fraction>1.00001)
    {
        printf("ERROR: TSC particle fraction inserted: %e (instead of =1). Aborting.\n", total_particle_fraction);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }

    return total_particle_fraction; // Return total particle fraction inserted for that particle (should always be 1 for fully periodic grid).
}


double W_TSC(double s)
{
    double w_tsc;

    if (s<=0.5) w_tsc=0.75-s*s;
    else if (s>1.5) w_tsc=0;
    else w_tsc=0.5*(1.5-s)*(1.5-s);

    return w_tsc;
}


/////////////////////////////////////////////
/////////////////////////////////////////////
/////////////////////////////////////////////




void compute_plane (struct plane_2D *Plane, int plane_number, int nx, int ny, int seed_block, int realization, int last_realization, int number_of_planes, double source_redshift, double plane_before_source, int scramble_mode, int species, int cell_embedding, int plane_assignment_mode)
{
		
    int i, j, k, r, upper_r_limit; // r is the local-block realization counter, running from 0 to upper_r_limit-1, where upper_r_limit<=parameters.seed_block. (Local-block is the clustering during snapshot reading to process several realizations (but not all assigned to this process) at the same time with one snapshot file read.)

    char plane_file[200];
    int redshift_tag;
    
    int nxny;
    
    nxny=nx*ny;
    

    // Initialize local variables containing grids for Fourier transforms:
	double **density_array;
	fftw_complex *Phi_ft_array, *Phi_array; // density contrast, gravitational potential, first and second derivatives of gravitational potential, Fourier transforms of all quantities.
	fftw_plan plan_density, plan_Phi_ft; // plans for FFTW Fourier Transforms.


	density_array = Matrix0(seed_block, nxny);
	
	if (realization+seed_block<=last_realization) upper_r_limit=seed_block;
	else upper_r_limit=last_realization-realization+1; // caps length of seed block when less than that realizations left.

	
	// Initialize density:
	printf("Process %d (superrank) %d (local process number): Initializing density contrast planes...\n", superrank, parameters.process_number);
	initialize_density_plane_multi(Plane, plane_number, density_array, nx, ny, seed_block, realization, last_realization, number_of_planes, source_redshift, plane_before_source, scramble_mode, species, cell_embedding, plane_assignment_mode);
	printf("Process %d (superrank) %d (local process number): Done initializing density constrast planes.\n", superrank, parameters.process_number);
	timing();
	fflush(stdout);
	


	printf("Process %d (superrank) %d (local process number): Size of allocated FFTW plan plan_density: %d.\n", superrank, parameters.process_number, ((int) seed_block* sizeof(fftw_plan)));
	
	
	Phi_ft_array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);
	Phi_array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny);
	assert(Phi_ft_array!=NULL);
	assert(Phi_array!=NULL);
	
	// Temporary:
	plan_density = fftw_plan_dft_2d(ny, nx, Phi_ft_array, Phi_array, FFTW_FORWARD, FFTW_ESTIMATE); // will actually do density into density_ft.
    plan_Phi_ft = fftw_plan_dft_2d(ny, nx, Phi_ft_array, Phi_array, FFTW_BACKWARD, FFTW_ESTIMATE);
	

  
    for (r=0; r<upper_r_limit; r++)
    {


        // First write density plane, then do Fourier transform and Poisson equation to obtain gravitational potential Phi, and write that.

        // Write density array:

        //////////////////////////////////////////////////////////////////////////
        // Writing density planes (only for a few realizations, since not needed for weak lensing map generation):
        printf("Writing density plane into FITS file...\n");

        redshift_tag=((int) floor(100.0*source_redshift));
        if (plane_before_source<0) sprintf(plane_file, "%s/%s_%s_%04dxy_%04dr_%03dp.%s", parameters.plane_output_path, parameters.density_basename, parameters.simulation_codename, nx, realization+r, plane_number, parameters.extension); // regular plane.
        else sprintf(plane_file, "%s/%s_%s_%04dxy_%04dr_%03dp_%04dz.%s", parameters.plane_output_path, parameters.density_basename, parameters.simulation_codename, nx, realization+r, plane_number, redshift_tag, parameters.extension); // last plane before source plane.

        long mapdims[2]={nx, ny};
#ifdef HAVE_FITS
        // NOTE: Uncomment the lines below if want to output density planes (used only for convergence_direct!=0 during ray tracing, otherwise potential planes are used, which are output further below.
        // writePlaneFITSimage_f(plane_file, 2, mapdims, density_array[r], plane_number);
        // if(parameters.realization+r<=10) writePlaneFITSimage_f(plane_file, 2, mapdims, imagearray, plane_number); // writes density planes for the first 10 realizations (for cluster mass estimation in line-of-sight investigation later).
#else
        exit(1); // alternatives currently not implemented, must have FITS (cfitsio library) installed.
#endif

        // Now calculate Potential with Fourier transforms and Poisson equation:

        // Zero out receiving arrays:
		for (i=0; i<nxny; i++)
		{
            Phi_ft_array[i]=density_array[r][i]; // transfering real density_array into fftw_complex array Phi_ft_array.
            Phi_array[i]=0.0;
		}


        // Caluculate Fourier Transform of density:
        if (feedback>2) printf("Executing FFTW-forward of density...\n");
        fflush(stdout);
        fftw_execute(plan_density);
        
        if (feedback>2) printf("Done executing FFTW-forward of density.\n");
        fflush(stdout);

        double H_0=0.1;   // H_0=0.1 here, since H_0 needs to be in units of km/s/kpc*h in the Poisson equation below, because Nabla operator in front of potential has spatial derivatives with length measures in h^(-1)*kpc (Plane[plane_number].binwidth_x is in kpc) !!!!!
        double Omega_m=parameters.Omega_m;
        double light_speed=299792.458; // speed of light in km/s. Keeping it in km/s here cancels nicely the H_0 in km/s/Mpc above. The residual /Mpc is intentional, because physical distances in this code are calculated in units of Mpc.

        double Three_OmH_0c2=3.0*Omega_m*(H_0/light_speed)*(H_0/light_speed); // 3*Omega_m*(H_0/c)^2, enters Poisson equation. Make sure density_array and density_ft_array are density contrast, i.e. mean density subtracted and divided.
	
        // Calculating Phi from FT-ed Poisson Equation, as well as second derivative of Phi in Fourier space:
        for (i=0; i<ny; i++)
        {
            for (j=0; j<nx; j++)
            {
                k=i*nx+j;
			
                if (i==0 && j==0) Phi_ft_array[k]=0; // setting zero mode to zero; not specified by Poisson Equation.
                else Phi_ft_array[k]=(Three_OmH_0c2*Phi_array[k])/((-4.0*sin(pi*j/nx)*sin(pi*j/nx)/(Plane[plane_number].binwidth_x*Plane[plane_number].binwidth_x))+(-4.0*sin(pi*i/parameters.ny)*sin(pi*i/parameters.ny)/(Plane[plane_number].binwidth_y*Plane[plane_number].binwidth_y)));  // continuous: /(k1sq*isq+k2sq*jsq); // Fourier Transformed Poisson Equation.
                // Phi_array on RHS is actually density_ft_array, we're just reusing the variable.
            }
        }
	
	
        // Transforming first and second derivative of Phi back to physical space:
        if (feedback>2) printf("Executing FFTW-backward of Phi, DPhi, and DDPhi...\n");
        fflush(stdout);
        // Temporary: // why is this line titled temporary?? Seems vital.
        fftw_execute(plan_Phi_ft);
		
        // Normalizing whole FT cycle, since FFTW is unnormalized:
        for (i=0; i<nxny; i++)
        {
            Phi_array[i] = Phi_array[i] / ((double) nxny);
        }
	
	

        printf("Process %d (superrank) %d (local process number): Done with Fourier.\n", superrank, parameters.process_number);
        timing();
        fflush(stdout);

        ////////////////////////////////////////////////////////////
        // Write potential plane into FITS image file:
        printf("Process %d (superrank) %d (local process number): Writing potential plane into FITS image file...\n", superrank, parameters.process_number);
	       
        redshift_tag=((int) floor(100.0*source_redshift));
        // printf("Redshift tag, source redshift: %d %e\n", redshift_tag, parameters.source_redshift);
        if (plane_before_source<0) sprintf(plane_file, "%s/%s_%s_%04dxy_%04dr_%03dp.%s", parameters.plane_output_path, parameters.potential_basename, parameters.simulation_codename, nx, realization+r, plane_number, parameters.extension); // regular plane.
        else sprintf(plane_file, "%s/%s_%s_%04dxy_%04dr_%03dp_%04dz.%s", parameters.plane_output_path, parameters.potential_basename, parameters.simulation_codename, nx, realization+r, plane_number, redshift_tag, parameters.extension); // last plane before source plane.
        // printf("Realization: %d, i: %d\n", parameters.realization, r);
	
        // Transfer complex array into real array before writing (note: density array will now actually contain the gravitational potential):
        for (i=0; i<ny; i++)
        {
            for (j=0; j<nx; j++)
            {
                k=2*(i*nx+j);
                density_array[r][i*nx+j]=((double *)Phi_array)[k];
            }
        }
	
        Plane[plane_number].particles_written[seed_block]=Plane[plane_number].particles_written[r]; // write parameter for this realization into last unused spot in array which is the write-out position in subsequents FITS file writing (this is not a very nice way of doing it).
        for (i=0; i<3; i++) // ditto for Rotation and Shift:
        {
            Plane[plane_number].rRot[i][seed_block]=Plane[plane_number].rRot[i][r];
            Plane[plane_number].rShift[i][seed_block]=Plane[plane_number].rShift[i][r];
        }
	
        
#ifdef HAVE_FITS
        writePlaneFITSimage_f(plane_file, 2, mapdims, density_array[r], plane_number); // writes a FITS image file with signed float pixel values.
#else
        exit(1); // Alternative file formats currently not implemented, must have FITS (cfitsio library) installed.
#endif
	
	
    } // ends r-writing loop over seeds

    
	// Free allocated FFTW plans and complex arrays:
	fftw_destroy_plan(plan_Phi_ft);
	fftw_destroy_plan(plan_density);

	fftw_free(Phi_ft_array);
	fftw_free(Phi_array);
	fftw_cleanup();
	

	//<AP> Failsafe to prevent stack corruption from crashing the execution
	seed_block = parameters.seed_block;
	//</AP>

	free_Matrix0(density_array, seed_block);

	printf("Process %d (superrank) %d (local process number): Done writing images.\n", superrank, parameters.process_number);
        timing();
	fflush(stdout);
    
}


void initialize_density_plane_multi(struct plane_2D *Plane, int plane_number, double** density_array, int nx, int ny, int seed_block, int realization, int last_realization, int number_of_planes, double source_redshift, int plane_before_source, int scramble_mode, int species, int cell_embedding, int plane_assignment_mode)
// WARNING: This function uses a global file pointer at the beginning, do not parallelize calls to this function (but the later part of this function can be parallelized).
{
	
    int snapshot_number, snapshot_file_number;
    int number_of_files_per_snapshot;
	int i, j, k, species_offset;
	double plane_zmax, plane_zmin;

	char input_fname[1000], input_fname2[1000];
	int file_number; // variable running over the snapshot box parts (box is devided into several files for memory reasons).
	double Position[3]; // used to be float.
	double halfboxsize;
	plane_zmax=0; plane_zmin=0; halfboxsize=0;
	Position[0]=0.0; Position[1]=0.0; Position[2]=0.0;

	double x1, x2, x3;
	int ii;
	int r, upper_r_limit; // index summing over realizations.
	double temp;

	double *particles_written = malloc(seed_block * sizeof(double));
    double particle_fraction;
					
	// Allocate random variable arrays:
	struct scramble *scrambled = malloc(seed_block * sizeof(struct scramble));
	double *random_number1 = malloc(seed_block * sizeof(double));
	double *random_number2 = malloc(seed_block * sizeof(double));
	double *random_number3 = malloc(seed_block * sizeof(double));
		
	double source_comoving_distance;
    int nxny;
    
    nxny=nx*ny;
    
	if (realization+seed_block<=last_realization) upper_r_limit=seed_block;
	else upper_r_limit=last_realization-realization+1; // caps length of seed block when less than that realizations left.

		
	snapshot_number=Plane[plane_number].snapshot;
	
	if (feedback>3) printf("Investigating Plane number, Snapshot Number %d %d \n", plane_number, snapshot_number); 
		

    // Turns snapshot boxes and shifts them, in order to avoid artifacts from reusing the same box:
    box_rotation_random_number_assignment(scrambled, random_number1, random_number2, random_number3, scramble_mode, upper_r_limit);
    // WARNING: The above function uses a global pointer to a file and is not thread safe! Do not parallelize using shared memory (OpenMP); distributed memory parallelization (MPI), where each process has its own independent memory content with all duplications, is ok.
    
		
    number_of_files_per_snapshot=1; // Default set, at least one snapshot file; gets updated in the following loop as soon as first file is read.
		
		for (file_number=0; file_number<number_of_files_per_snapshot; file_number++)
		{
		  
			snapshot_file_number=snapshot_number;
			sprintf(input_fname, "%s/%s_%03d", parameters.snapshot_path, parameters.snapshot_name, snapshot_file_number); // only place where snapshot_file_number needs to be used, all rest of code done in internal standard snapshot_number count, starting at 0 with the snapshot farthest away.
			
			printf("Process %d (superrank) %d (local process number): Begin loading snapshot...\n%s\n", superrank, parameters.process_number, input_fname);
			timing();
			fflush(stdout);
			if (parameters.endianness_set==0)
			  {
			    // Determine endianness of Gadget-2 snapshot (do this only once per run):
			    parameters.byteswap=snapshotEndianness(input_fname, file_number);
			    parameters.endianness_set=1;
				// parameters.byteswap explicitly stores the relative endianness of the Gadget-2 snapshot.
			    printf("Done setting endianness.\n");
			  }
			fflush(stdout);
			load_snapshot_multi(input_fname, file_number, snapshot_number, plane_number);  // snapshot_number is the internal count used by IG here, not snapshot_file_number. (After rewrite update for snapshot skipping, these became the same.)
            // WARNING GLOBAL: The above function uses global structures *Snapshot (contains particle positions, among other things), and header1.
            // Because of the above function (and its allocation of parts of Snapshot within), Snapshot must be a global variable in initialize_density_plane_multi.
			printf("Process %d (superrank) %d (local process number): Done loading snapshot.\n", superrank, parameters.process_number);
			timing();
			fflush(stdout);

			if (file_number==0) // First file of a snapshot read in:
			{
				number_of_files_per_snapshot=Snapshot[snapshot_number].files; // read number of files from snapshot header and transfer to loop parameter.
				
				// Quantities variable between snapshots (but const. among files from same snapshot):
				parameters.boxsize=Snapshot[snapshot_number].Boxsize; // NOTE: parameters.boxsize should be phased out to make sure code can handle simulations with variable boxsizes.
				halfboxsize=0.5*parameters.boxsize;
				// Quantities constant for all snapshots of one simulation:
				parameters.Omega_m=Snapshot[snapshot_number].Omega0;
				parameters.Omega_Lambda=Snapshot[snapshot_number].OmegaLambda;
				parameters.H_0=Snapshot[snapshot_number].H_0;
				parameters.h=parameters.H_0/100.0; // parameters.h is needed for proper results of deflection angle and convergence (although they are independent, the finite sum needs it). 
				parameters.w0=Snapshot[snapshot_number].w0;
				parameters.wa=Snapshot[snapshot_number].wa;
				parameters.ns=Snapshot[snapshot_number].ns;
				parameters.sigma_8=Snapshot[snapshot_number].sigma_8;

				initialize_darkenergy(); // initializes w(z) profile for calculation of comoving distances (those which are not taken from Gadget snapshots directly, like for the sources).
				setup_plane(Plane, Snapshot, plane_number); // initializes values in Plane[] structure for this plane.
				// Also need data for neighboring planes, to get slicing right:
				if (plane_number>0 && Plane[plane_number-1].set!=1)
				  {
				    sprintf(input_fname2, "%s/%s_%03d", parameters.snapshot_path, parameters.snapshot_name, Plane[plane_number-1].snapshot);
				    // printf("Will read closer header:\n%s\n", input_fname2);
				    load_snapshot_multi_header(input_fname2, file_number, Plane[plane_number-1].snapshot);
				    setup_plane(Plane, Snapshot, plane_number-1);
				  }
				if (plane_number<number_of_planes-1 && Plane[plane_number+1].set!=1)
				  {
				    sprintf(input_fname2, "%s/%s_%03d", parameters.snapshot_path, parameters.snapshot_name, Plane[plane_number+1].snapshot);
				    // printf("Will read farther header:\n%s\n", input_fname2);
				    load_snapshot_multi_header(input_fname2, file_number, Plane[plane_number+1].snapshot);
                                    setup_plane(Plane, Snapshot, plane_number+1);
				  }


				if (plane_before_source!=-1)
				  {
				    source_comoving_distance=calculate_comoving_distance((1.0/(source_redshift+1.0))); // if end plane, need to calculate comoving distance of sources, now that we know the cosmological parameters.
				  }

				if (feedback>3) printf("Cosmological Parameter Transfer from Snapshot: Boxsize=%e, Omega_m=%e, Omega_Lambda=%e, H_0=%e.\n", parameters.boxsize, parameters.Omega_m, parameters.Omega_Lambda, parameters.H_0);
				
				
				if (plane_number==0)
				{
					plane_zmax=(Plane[plane_number+1].comoving_distance-Plane[plane_number].comoving_distance)/2.0+halfboxsize;
					plane_zmin=halfboxsize-Plane[plane_number].comoving_distance;
				}
				else if (plane_number==number_of_planes-1)
				{
					plane_zmin=-(Plane[plane_number].comoving_distance-Plane[plane_number-1].comoving_distance)/2.0+halfboxsize;
					plane_zmax=(Plane[plane_number].comoving_distance-Plane[plane_number-1].comoving_distance)/2.0+halfboxsize;
				}
				else
				{	
					plane_zmax=(Plane[plane_number+1].comoving_distance-Plane[plane_number].comoving_distance)/2.0+halfboxsize;
					plane_zmin=-(Plane[plane_number].comoving_distance-Plane[plane_number-1].comoving_distance)/2.0+halfboxsize;
				}


				// Adjustment for last plane before source plane (overwrites above plane_zmax with comoving distance of source plane):
				if (plane_before_source>=0) plane_zmax=(source_comoving_distance-Plane[plane_number].comoving_distance)+halfboxsize;
				// The source plane always lies at greater comoving distance than the center of the last plane before source plane. The source plane is not actually a file in this code, it's the farthest limit from which particles are included into the generation of the last plane before source plane.

				Plane[plane_number].catchment_far=plane_zmax-halfboxsize;
				Plane[plane_number].catchment_close=plane_zmin-halfboxsize;
				if (plane_before_source>=0)
				  {
				    printf("Superrank: %d, Plane catchment far, Plane comoving distance, source_comoving_distance: %e %e %e \n", superrank, Plane[plane_number].catchment_far, Plane[plane_number].comoving_distance, source_comoving_distance);
				    fflush(stdout);
				  }
			}
			else if (number_of_files_per_snapshot!=Snapshot[snapshot_number].files)
			{
					printf("ERROR: Inconsistent number of files indicated in headers of files composing one snapshot. Aborting.");
					exit(1);
			}

			if ((plane_zmin<0.0 || plane_zmax>parameters.boxsize) && file_number==0) printf("WARNING: Gap between snapshot boxes at Box %d (Plane %d) and the one behind it.\n", snapshot_number, plane_number);

			// From here on want to keep for every process:
		
			if (file_number==0) 
			{
				for (r=0; r<seed_block; r++) particles_written[r]=0.0;
                particle_fraction=0.0; // this is new variable for OpenMP version.
				if (feedback>2) printf("Comoving distance of plane: %d %d %e %e %e \n", number_of_planes, plane_number, Plane[plane_number].comoving_distance, Plane[plane_number-1].comoving_distance, Plane[plane_number+1].comoving_distance);
				if (feedback>2) printf("Plane close and far Einzugs comoving distance boundary: %e %e %e \n", plane_zmin, plane_zmax, halfboxsize);
			}
			
			species_offset=0;
			
			printf("Inserting realizations...\n");

			for (r=0; r<upper_r_limit; r++) // loop over realizations:
			{
			
			  if (file_number==0) // do scaling once, when processing first file.
			  {
			    random_number1[r] *= Plane[plane_number].boxsize; 
			    random_number2[r] *= Plane[plane_number].boxsize; 
			    random_number3[r] *= Plane[plane_number].boxsize; 
			    printf("Random Numbers for box shift: %e %e %e \n", random_number1[r], random_number2[r], random_number3[r]);


			    // Record rotation and shift parameters (x axis =1, y axis =2, z axis =3, mirror an axis: sign (minus)):
			    // also incorporates switch from right-handed (z, y, x) coordinate system to (x, y, z), where x is horizontally to the right, y is vertically up, and z is away from observer.
			    Plane[plane_number].rRot[0][r]=scrambled[r].signx1*(scrambled[r].x1+1);
			    Plane[plane_number].rRot[1][r]=scrambled[r].signx2*(scrambled[r].x2+1);
			    Plane[plane_number].rRot[2][r]=scrambled[r].signx3*(scrambled[r].x3+1);
			    Plane[plane_number].rShift[0][r]=random_number1[r]/1000.0; // Store shift in Mpc/h, so convert from kpc/h used by Gadget-2.
			    Plane[plane_number].rShift[1][r]=random_number2[r]/1000.0;
                Plane[plane_number].rShift[2][r]=random_number3[r]/1000.0;

			  }

				if (feedback>1) printf("Inserting realization %d...\n", r);
				fflush(stdout);
			
			for (i=species;i<(species+1);i++)
			{
				if (species!=0)
				{
					for (k=0; k<species;k++)
					{
						species_offset=species_offset+header1.npart[k];
					}
				}

#ifdef HAVE_OpenMP
# pragma omp parallel private(j, x1, x2, x3, particle_fraction) shared(density_array, particles_written)
// WARNING: Must still resolve parallel writing to density_array[r]. The current implementation may be too blocking and may run much slower than the non-OpenMP version.
// WARNING2: particles_written is not properly computed in the OpenMP version. This doesn't affect runs other than that the wrong number for particles written into a density/potential plane will appear in the header of the plane file (it is currently not further used by any codes and just information for the user.
                {
# pragma omp for schedule (static)
				for (j=(species_offset+1);j<=(species_offset+header1.npart[i]);j=j+1)
				{

					if (scramble_mode!=0) // scramble coordinates of particles by axes interchange (rotation), axes sign (mirroring), and shifts with periodic completion:
					{

						x1=scrambled[r].signx1*Snapshot[snapshot_number].P[j].Pos[scrambled[r].x1]+random_number1[r];
						x2=scrambled[r].signx2*Snapshot[snapshot_number].P[j].Pos[scrambled[r].x2]+random_number2[r];
						x3=scrambled[r].signx3*Snapshot[snapshot_number].P[j].Pos[scrambled[r].x3]+random_number3[r];


						// Readjust particle coordinates so that they are within box again:
                        while(x1<0) x1 += parameters.boxsize;
                        while(x1>=parameters.boxsize) x1 -= parameters.boxsize;
                        while(x2<0) x2 += parameters.boxsize;
                        while(x2>=parameters.boxsize) x2 -= parameters.boxsize;
                        while(x3<0) x3 += parameters.boxsize;
                        while(x3>=parameters.boxsize) x3 -= parameters.boxsize;
                    
					}
                    else // if no snapshot box rotations and shifts:
                    {
                        x1=Snapshot[snapshot_number].P[j].Pos[0];
                        x2=Snapshot[snapshot_number].P[j].Pos[1];
                        x3=Snapshot[snapshot_number].P[j].Pos[2];
                        
                    }
                    
					// Write particle into density array (via cloud-in-cell method):
					if (x1>=plane_zmin && x1<plane_zmax)
                    {
					    x1-=halfboxsize;
					    x2-=halfboxsize;
					    x3-=halfboxsize;
                          

					    if (cell_embedding==1) particle_fraction+=cloud_in_cell(x3, x2, density_array[r], nx, ny, Plane[plane_number].binwidth_x, Plane[plane_number].binwidth_y, plane_assignment_mode);
					    else if (cell_embedding==2) particle_fraction+=TSC_assign(x3, x2, density_array[r], nx, ny, Plane[plane_number].binwidth_x, Plane[plane_number].binwidth_y);
#pragma omp critical (particlefraction)
                        {
                            particles_written[r]+=particle_fraction;
                        }
                    }
                } // end j loop.
                
                } // end of OpenMP parallelization.
#else
                // This is the way it was before OpenMP. It will become obsolete once the OpenMP version is verified to have better performance (which will likely still need some changes above):
                for (j=(species_offset+1);j<=(species_offset+header1.npart[i]);j=j+1)
                {
                    
					if (scramble_mode!=0) // scramble coordinates of particles by axes interchange (rotation), axes sign (mirroring), and shifts with periodic completion:
					{
                        
						x1=scrambled[r].signx1*Snapshot[snapshot_number].P[j].Pos[scrambled[r].x1];
						x2=scrambled[r].signx2*Snapshot[snapshot_number].P[j].Pos[scrambled[r].x2];
						x3=scrambled[r].signx3*Snapshot[snapshot_number].P[j].Pos[scrambled[r].x3];
                        
						Position[0]=x1+random_number1[r];
						Position[1]=x2+random_number2[r];
						Position[2]=x3+random_number3[r];
                        
						// Readjust particle coordinates so that they are within box again:
						for (ii=0; ii<3; ii++)
						{
                            while(Position[ii]<0) Position[ii] += parameters.boxsize;
                            while(Position[ii]>=parameters.boxsize) Position[ii] -= parameters.boxsize;
						}
					}
                    else // if no snapshot box rotations and shifts:
                    {
                        Position[0]=Snapshot[snapshot_number].P[j].Pos[0];
                        Position[1]=Snapshot[snapshot_number].P[j].Pos[1];
                        Position[2]=Snapshot[snapshot_number].P[j].Pos[2];
                        
                    }
                    
					// Write particle into density array (via cloud-in-cell method):
					if (Position[0]>=plane_zmin && Position[0]<plane_zmax)
                    {
					    Position[0]-=halfboxsize;
					    Position[1]-=halfboxsize;
					    Position[2]-=halfboxsize;
                        
                        
					    if (cell_embedding==1) particles_written[r]=particles_written[r]+cloud_in_cell(Position[2], Position[1], density_array[r], nx, ny, Plane[plane_number].binwidth_x, Plane[plane_number].binwidth_y, plane_assignment_mode);
					    else if (cell_embedding==2) particles_written[r]=particles_written[r]+TSC_assign(Position[2], Position[1], density_array[r], nx, ny, Plane[plane_number].binwidth_x, Plane[plane_number].binwidth_y);
                        
                    }

    
				}
#endif
				
			}
			}
			
			free_snapshot_arrays(snapshot_number); // Free snapshot file memory before reading in next snapshot file.
		}

		// Swap header coordinates x and z (because used differently in post-processing):
		for (r=0; r<upper_r_limit; r++)
		  {
		    temp=Plane[plane_number].rRot[0][r];
		    Plane[plane_number].rRot[0][r]=Plane[plane_number].rRot[2][r];
		    Plane[plane_number].rRot[2][r]=temp;
		    temp=Plane[plane_number].rShift[0][r];
		    Plane[plane_number].rShift[0][r]=Plane[plane_number].rShift[2][r];
		    Plane[plane_number].rShift[2][r]=temp;
		  }

		if (feedback>2) for (r=0; r<upper_r_limit; r++) printf("Particles written %e\n", particles_written[r]);		
		fflush(stdout);


			for (r=0; r<upper_r_limit; r++)
		  {
		    Plane[plane_number].particles_written[r]=particles_written[r];

		    for (i=0;i<nxny;i++) density_array[r][i]=density_array[r][i]*nxny*(Plane[plane_number].boxsize/1000.0)/Snapshot[snapshot_number].NumPartTotal[species]-(fabs(Plane[plane_number].catchment_far/1000.0)+fabs(Plane[plane_number].catchment_close/1000.0)); // Divisions by 1000.0 because want chi_thickness units in Mpc/h, and values are in kpc/h natively on first run of Inspector Gadget.

		  }

		// Deallocate arrays (have mostly to do with random numbers and box rotations and shifts):
		free(particles_written);
		free(scrambled);
		free(random_number1);
		free(random_number2);
		free(random_number3);


		fflush(stdout);
}



// Turns snapshot boxes and shifts them, in order to avoid artifacts from reusing the same box:
void box_rotation_random_number_assignment(struct scramble *scrambled, double *random_number1, double *random_number2, double *random_number3, int scramble_mode, int upper_r_limit)
// WARNING GLOBAL: This function uses global file pointer input_file_ran. It is not thread safe!
{
    char line[2000];
    int line_length;
	float re1, re2, re3, re4, re5, re6, re7, re8, re9;
	double ra1, ra2, ra3, ra4, ra5, ra6, ra7, ra8, ra9;
    
    int r;
    
    for (r=0; r<upper_r_limit; r++)
	{
        
        if ((line_length=fgetline(input_file_ran, line, sizeof(line))) > 0) // WARNING: This line uses global file pointer input_file_ran.
        {
            
            if (sscanf(line, "%e %e %e %e %e %e %e %e %e", &re1, &re2, &re3, &re4, &re5, &re6, &re7, &re8, &re9) == 9) // scan line for proper format of random number file (9 floating point random numbers).
            {
                ra1=0; ra2=0; ra3=0; ra4=0; ra5=0; ra6=0; ra7=0; ra8=0; ra9=0;
                ra1=((double) re1); ra2=((double) re2); ra3=((double) re3); ra4=((double) re4); ra5=((double) re5); ra6=((double) re6); ra7=((double) re7); ra8=((double) re8); ra9=((double) re9); // convert to doubles.
                printf("Random Number Read In: %e %e %e %e %e %e %e %e %e\n", ra1, ra2, ra3, ra4, ra5, ra6, ra7, ra8, ra9);
            }
            else // if line does not have proper format:
            {
                printf("ERROR: Ran out of random numbers in random number file! Aborting.\n");
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD, 1);
                exit(1);
            }
    
            scrambled[r]=Scrambler(scramble_mode, ra1, ra2, ra3, ra4, ra5, ra6);
	
            if (scramble_mode!=0)
            {
                random_number1[r]=ra7; random_number2[r]=ra8; random_number3[r]=ra9;
                if (feedback>3) printf("Random Numbers: %e %e %e \n", random_number1[r], random_number2[r], random_number3[r]);
            }
            else
            {
                random_number1[r]=0.0; random_number2[r]=0.0; random_number3[r]=0.0;
            }
            
        }
        else // if line_length not > 0:
        {
            printf("ERROR: Ran out of random numbers in random number file! Aborting.\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(1);
        }
        
	} // end of r for-loop.
    
}




void setup_plane(struct plane_2D *Plane, struct snapshot_control *Snapshot, int plane_number) // requires match_planes() to have been called earlier, such that Plane[plane_number].snapshot is already set properly (the rest of the Plane structure is set by this function).
{
  int snapshot_number, i;
  
  snapshot_number=Plane[plane_number].snapshot; // this has already been initialized earlier for all planes via match_planes() in main.c.

  if (Snapshot[snapshot_number].set!=1)
    {
      printf("ERROR: Snapshot[snapshot_number] has to be initialized before setup_plane(plane_number) can be called.\nError occurred for plane=%d and snapshot=%d.\n", plane_number, snapshot_number);
      fflush(stdout);
      exit(1);
    }

  for (i=0;i<6;i++)
    {
      Plane[plane_number].NumPartTotal[i]=Snapshot[snapshot_number].NumPartTotal[i];
      Plane[plane_number].mass[i]=Snapshot[snapshot_number].mass[i];
    }
  Plane[plane_number].scale_factor=Snapshot[snapshot_number].Time;
  Plane[plane_number].Omega0=Snapshot[snapshot_number].Omega0;
  Plane[plane_number].OmegaLambda=Snapshot[snapshot_number].OmegaLambda;
  Plane[plane_number].H_0=Snapshot[snapshot_number].H_0;
  Plane[plane_number].boxsize=Snapshot[snapshot_number].Boxsize;
  Plane[plane_number].binwidth_x=Plane[plane_number].boxsize/((double) parameters.nx);
  Plane[plane_number].binwidth_y=Plane[plane_number].boxsize/((double) parameters.ny);

  Plane[plane_number].w0=Snapshot[snapshot_number].w0;
  Plane[plane_number].wa=Snapshot[snapshot_number].wa;
  Plane[plane_number].ns=Snapshot[snapshot_number].ns;
  Plane[plane_number].sigma_8=Snapshot[snapshot_number].sigma_8;
  Plane[plane_number].initial_condition=Snapshot[snapshot_number].initial_condition;
  Plane[plane_number].comoving_distance=Snapshot[snapshot_number].comoving_distance;

  Plane[plane_number].set=1; // this sets flag that parameters for this plane have been initialized properly.

}
