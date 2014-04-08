/*
 *  weak_lensing_multi.c
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 4/18/08.
 *  Copyright 2008. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <time.h>

#include <mpi.h>

#include <complex.h>
#include <fftw3.h>

#include "main.h"
#include "mathematics.h"
#include "allocation.h"
#include "2D-plane_multi.h"
#include "preload_planes.h"
#include "weak_lensing_multi.h"
#include "mpi_support.h"
#include "darkenergy.h"
#include "comoving_distance.h"
#include "galaxy.h"
#include "initialization.h"


#ifdef HAVE_OpenMP
#include <omp.h>
#endif


// #ifdef HAVE_FITS
#include "fits.h"
// #endif

// Number of OpenMP threads (not needed, can do with environment variable OMP_NUM_THREADS):
// #define  NR_OF_OMP_THREADS 64


#define MAX_RANDOM_NUMBER 10000000.0
// make sure the above random number is the same as the one used in the random number file generation application "Random Numbers for Inspector Gadget 5".


#define CONSIST (0) 
// consistency shift with previous top down y axis indexing. set to zero for production runs.

// WARNING: In this file, ray tracing arrays (like A11, F11, etc.) with Xbin and Ybin as index counters have FORTRAN ordering of indices by accident. Xbin still stands for x axis (horizontal to right) and Ybin for y axis (vertically up). Calculation is not impacted, except for writeout into writeout_array this needs to be taken into account by reversing the indices, since the 1D-(pseudo-2D) writeout_array feeding the FITS file writing routine needs the numbers in standard C ordering (Y*imagedim+X, i.e. [Y][X] as the C arrays are realized). The indexing of the 1's in A11, theta1, etc. is still correct, 1 denoting first component or x component as appropriate).


// double* potential_array;


// will need:
// int nx, int ny, double binwidth_x, double binwidth_y, int ray_tracing, int plane_assignment_mode, int raypoint_averaging, struct Plane
void weak_lensing2(int NbinsX_global, int NbinsY_global, int nx, int ny, int ray_tracing, int plane_assignment_mode, int raypoint_averaging, struct plane_2D *Plane, int convergence_direct, int galaxy_catalogue_type, int realization, double survey_angle_in_rad, double source_comoving_distance, int number_of_planes, int plane_before_source, double source_redshift, int plane_shift, int preload_planes, int number_of_plane_realizations, int first_sim_ic, MPI_Comm sim_comm, MPI_Win *plane_storage_window)
{

 int pixelcounter, Xbin, Ybin, plane_number, planes;
 double survey_xres, survey_yres;
 double x, y;
 double delPhi1, delPhi2;
 double chi, a, temp;
 
 double **theta1_ini, **theta2_ini, **theta1, **theta2, **f1, **g1, **f2, **g2; // NxN arrays
 //struct matrix2x2 **A, **F, **G; // NxNx2x2 arrays
 double **A11, **A12, **A21, **A22, **F11, **F12, **F21, **F22, **G11, **G12, **G21, **G22;
 //struct matrix2x2 U, IdentityM;
 double U11, U12, U21, U22, IdentityM11, IdentityM12, IdentityM21, IdentityM22;
 double *writeout_array;
 // double **PrevA11, **PrevA12, **PrevA21, **PrevA22, **SingleA11, **SingleA12, **SingleA21, **SingleA22, **CumulA11, **CumulA12, **CumulA21, **CumulA22;
 // double *ray_pos;
    double **SC; // array containing source comoving distances of galaxies (or rays if do maps to a certain redshift plane).

    
    FILE *galaxy_input_file, *galaxy_output_file; // files to read in and write galaxy catalogues (used if do not want to do simple maps).
    char galaxy_input_filename[2000], galaxy_output_filename[2000];
    struct galaxy_struct galaxy;
    int number_of_galaxies, galaxy_number, number_of_redshifts_per_galaxy;
    int last_plane_number;
    int galaxy_redshift_counter;
    int dummy_int; // placeholder variable not needed for anything, but required for some function calls.
    
    int nxny;
    nxny=nx*ny;
    
    double *potential_array;
    
    int NbinsX, NbinsY;
    

    ///////////////////////////////
    // Needed only for convergence direct:
    double convergence_prefactor;
    double H_0=100.0; // H_0 in km/s/Mpc, but still in units of Mpc/h, so the factor of h cancels.
    double light_speed=299792.458; // speed of light in km/s.
    convergence_prefactor=3.0/2.0*parameters.Omega_m*(H_0/light_speed)*(H_0/light_speed);
    //////////////////////////////

    

 if (convergence_direct!=0) ray_tracing=0; // if want to compute convergence directly from matter density (rather via gravitational potential) then potential plane is actually density plane, and deflection angle cannot be estimated.
    

 // Allocate arrays (matrices and array-matrices):
 printf("Allocating weak lensing map creation arrays (includes matrix arrays)...\n");
 fflush(stdout);
    
    

 
 int ThisTask, number_of_subprocesses;
 MPI_Comm_rank(map_comm, &ThisTask);
 MPI_Comm_size(map_comm, &number_of_subprocesses);


    
 if (galaxy_catalogue_type==0)
 {
     NbinsX=NbinsX_global;
     NbinsY=NbinsY_global/number_of_subprocesses;
 }
 else // WARNING: subprocesses split not properly implemented for this case yet (would have to tinker with skipping lines in catalogue file during galaxy read in, otherwise it works).
 {

     // Input and output filenames and paths for galaxy catalogue: 
     sprintf(galaxy_input_filename, "%s/%s%d.txt", parameters.galaxy_catalogue_path, parameters.galaxy_catalogue_basename, parameters.galaxy_subfield);
     sprintf(galaxy_output_filename, "%s/%s/%s%d_WL-only_%s_%04dxy_%04dr.%s", parameters.galaxy_catalogue_output_path, parameters.simulation_codename, parameters.galaxy_catalogue_basename, parameters.galaxy_subfield, parameters.simulation_codename, nx, realization, parameters.extension);
     
     number_of_galaxies=get_number_of_galaxies(galaxy_input_filename, &number_of_redshifts_per_galaxy);
     // The above line determines automatically how many galaxies are in the catalogue file and how many redshifts are stored for them (make sure there is the same number of redshifts for each galaxy, as the function only returns the maximum number of redshifts detected).
     
     NbinsX=floor(sqrt(number_of_galaxies*number_of_redshifts_per_galaxy))+1;
     NbinsY=NbinsX/number_of_subprocesses;
 }
    
 /*
 printf("NbinsX %d NbinsY %d parameters.NbinsX %d parameters.NbinsY %d parameters.nx %d parameters.ny %d, parameters.nxny %d\n", NbinsX, NbinsY, parameters.NbinsX, parameters.NbinsY, parameters.nx, parameters.ny, parameters.nxny);


 if (parameters.NbinsX!=2048 || parameters.NbinsY!=2048 || parameters.nx!=4096 || parameters.ny!=4096)
   {
     printf("ERROR: NbinsX %d NbinsY %d parameters.NbinsX %d parameters.NbinsY %d parameters.nx %d parameters.ny %d, parameters.nxny %d\n", NbinsX, NbinsY, parameters.NbinsX, parameters.NbinsY, parameters.nx, parameters.ny, parameters.nxny);
     exit(1001);
   }
 */

 if (NbinsY_global%number_of_subprocesses!=0)
   {
     printf("ERROR: Number of Pixels of WL map in vertical dimension (NbinsY) must be a multiple of number_of_subprocesses. Aborting.\n");
     MPI_Abort(MPI_COMM_WORLD, 1);
     exit(1);
   }

 // The line below should be included only by ThisTask==0 and then communicated to the others.
 potential_array=Vector0(nxny);


 theta1_ini=Matrix0(NbinsY, NbinsX);
 theta2_ini=Matrix0(NbinsY, NbinsX);
 theta1=Matrix0(NbinsY, NbinsX);
 theta2=Matrix0(NbinsY, NbinsX);
 f1=Matrix0(NbinsY, NbinsX);
 f2=Matrix0(NbinsY, NbinsX);
 g1=Matrix0(NbinsY, NbinsX);
 g2=Matrix0(NbinsY, NbinsX);
 

 A11=Matrix0(NbinsY, NbinsX);
 A12=Matrix0(NbinsY, NbinsX);
 A21=Matrix0(NbinsY, NbinsX);
 A22=Matrix0(NbinsY, NbinsX);

 F11=Matrix0(NbinsY, NbinsX);
 F12=Matrix0(NbinsY, NbinsX);
 F21=Matrix0(NbinsY, NbinsX);
 F22=Matrix0(NbinsY, NbinsX);

 G11=Matrix0(NbinsY, NbinsX);
 G12=Matrix0(NbinsY, NbinsX);
 G21=Matrix0(NbinsY, NbinsX);
 G22=Matrix0(NbinsY, NbinsX);
    
 // Source Galaxy Comoving Distances:
    SC=Matrix0(NbinsY, NbinsX);
    

 // For ray positions on planes:

 // RAY-PLANES:  ray_pos=Vector0(2*parameters.NbinsX*parameters.NbinsY);


 // Use for line-of-sight investigation readouts (these work and can be activated as needed, including all occurences below and corresponding writeout funtion calls):

 /*
 PrevA11=Matrix0(parameters.NbinsX, parameters.NbinsY);
 PrevA12=Matrix0(parameters.NbinsX, parameters.NbinsY);
 PrevA21=Matrix0(parameters.NbinsX, parameters.NbinsY);
 PrevA22=Matrix0(parameters.NbinsX, parameters.NbinsY);
 

 SingleA11=Matrix0(parameters.NbinsX, parameters.NbinsY);
 SingleA12=Matrix0(parameters.NbinsX, parameters.NbinsY);
 SingleA21=Matrix0(parameters.NbinsX, parameters.NbinsY);
 SingleA22=Matrix0(parameters.NbinsX, parameters.NbinsY);

 
 CumulA11=Matrix0(parameters.NbinsX, parameters.NbinsY);
 CumulA12=Matrix0(parameters.NbinsX, parameters.NbinsY);
 CumulA21=Matrix0(parameters.NbinsX, parameters.NbinsY);
 CumulA22=Matrix0(parameters.NbinsX, parameters.NbinsY);
 */

 
 // if (ThisTask==0) writeout_array=Vector0(parameters.NbinsX*parameters.NbinsY);

 printf("Done allocating weak lensing map creation arrays.\n");
 fflush(stdout);
 timing();

 ////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
 
 survey_xres=survey_angle_in_rad/((double) NbinsX_global);
 survey_yres=survey_angle_in_rad/((double) NbinsY_global);
 // IdentityM.m11=1.0; IdentityM.m12=0.0; IdentityM.m21=0.0; IdentityM.m22=1.0;
 IdentityM11=1.0; IdentityM12=0.0; IdentityM21=0.0; IdentityM22=1.0;

 
        if (galaxy_catalogue_type!=0)
        {
            galaxy_input_file=fopen(galaxy_input_filename, "r");
        }
    
    galaxy_number=0;
    galaxy_redshift_counter=number_of_redshifts_per_galaxy;
    
    printf("Redshifts per galaxy: %d\n", number_of_redshifts_per_galaxy);

  
    
	// Values on Plane 0:
	for (Ybin=0; Ybin<NbinsY; Ybin++) {
	for (Xbin=0; Xbin<NbinsX; Xbin++) {	

        
            if (galaxy_catalogue_type!=0 && galaxy_number<number_of_galaxies && galaxy_redshift_counter>=number_of_redshifts_per_galaxy)
            {
                get_galaxy_parameters(galaxy_input_file, &galaxy, &dummy_int); // read in galaxies one by one from galaxy catalogue file.
                galaxy_redshift_counter=0;
            }
        
        if (galaxy_catalogue_type==0)
        {
        
            theta1_ini[Ybin][Xbin]=((double) (Xbin-0.5*NbinsX_global)*survey_xres); // initial angle of lightray as it hits first plane traveling in reverse (and the angle at which it is eventually seen by observer)
            theta2_ini[Ybin][Xbin]=((double) ((Ybin+ThisTask*NbinsY)-0.5*NbinsY_global)*survey_yres);
            
            SC[Ybin][Xbin]=source_comoving_distance; // all light rays have same comoving distance if creating regular lensing maps.
        }
        else
        {
            if (galaxy_number<number_of_galaxies)
            {
                theta1_ini[Ybin][Xbin]=galaxy.theta[0]; // apparent angle in flat sky (with respect to center of map (field of view) of galaxy from galaxy catalogue).
                theta2_ini[Ybin][Xbin]=galaxy.theta[1];
                
                if (galaxy.redshift[galaxy_redshift_counter]>source_redshift) galaxy.redshift[galaxy_redshift_counter]=source_redshift; // safeguards against too distant galaxies by placing them at the maximally allowed simulated redshift (can skip those in the analysis later if one wants).
                SC[Ybin][Xbin]=get_chi_in_Mpc(galaxy.redshift[galaxy_redshift_counter]);
                // Note: must call the spline-based function above to get the comoving distance for each galaxy based on its redshift. Using the line below instead would work too slowly (computes the comoving distance integral from scratch).
                // SC[Ybin][Xbin]=calculate_comoving_distance((1.0/(galaxy.redshift[galaxy_redshift_counter]+1.0)))/1000.0; // comoving distance of galaxy, computed from galaxy redshift from entry in galaxy catalogue. Comoving distance caclulation function computes in kpc/h, thus divide by 1000.
                // Do not use the above commented out line (it's correct, but too slow for this purpose).
                galaxy_redshift_counter++; // increment now, since have just placed a redshift.
            }
            else
            {
                theta1_ini[Ybin][Xbin]=0.0;
                theta2_ini[Ybin][Xbin]=0.0;
                SC[Ybin][Xbin]=0.0;
                
            }
            
        }
        
            
		theta1[Ybin][Xbin]=theta1_ini[Ybin][Xbin];
		theta2[Ybin][Xbin]=theta2_ini[Ybin][Xbin];
		f1[Ybin][Xbin]=0.0;
		g1[Ybin][Xbin]=0.0;
		f2[Ybin][Xbin]=0.0;
		g2[Ybin][Xbin]=0.0;

		A11[Ybin][Xbin]=1.0; A12[Ybin][Xbin]=0.0; A21[Ybin][Xbin]=0.0; A22[Ybin][Xbin]=1.0;
		F11[Ybin][Xbin]=0.0; F12[Ybin][Xbin]=0.0; F21[Ybin][Xbin]=0.0; F22[Ybin][Xbin]=0.0;
		G11[Ybin][Xbin]=0.0; G12[Ybin][Xbin]=0.0; G21[Ybin][Xbin]=0.0; G22[Ybin][Xbin]=0.0;
        
        if (galaxy_redshift_counter==number_of_redshifts_per_galaxy) galaxy_number++;
        
	}
	}

    
    if (galaxy_catalogue_type!=0) fclose(galaxy_input_file); // Close galaxy input file. Will be reopened again in writeout routine.
    
    // Obsolete lines from Inspector Gadget version which still had terminal planes for arbitrary source redshifts.
    /*
	planes=parameters.number_of_planes+1; // This includes all intermediate planes, as well as the source plane (which is not a physical FITS plane per se). If change this line, have to change line further down as well, where parameters.number_of_planes is restored from planes (the line with the bad programming note after save_WL_results() is called).

	if (planes<parameters.plane_before_source+1)
	  {
	    printf("ERROR: Last plane before source is disconnected from contiguous plane stack starting at observer. Weak lensing maps will be wrong. Aborting.\n"); 
	    exit(1);
	  }
     */
    planes=number_of_planes; // in new Inspector Gadget version without terminal planes use this value instead. Make sure parameters.number_of_planes is set to emcompass at least one plane beyond parameters.source_comoving_distance when doing full maps and beyond the furthest galaxy in the galaxy catalogue, or a segmentation fault will occur when Planes[plane_number+1] is called in the loops below (this outmost plane does not necessarily have to satisfy the field of view opening, because no interpolation will be done on it, only the comoving distance of it will be read. Save unnecessary read-in of unused planes by choosing the number of planes close to the number minimally needed (typically ~49 planes for maximum redshift of z=2 with 80Mpc/h plane spacings). The desired redshift can be read in the headers of FITS files of the planes, so determining the optimal farthest plane is easy.

 for (plane_number=1; plane_number<planes; plane_number++)
 {
   timing();
   printf("Calculating weak lensing maps including planes up to Plane %d...\n", plane_number-1);
     
     
   if (ThisTask==0) load_potential_plane(potential_array, plane_number-1, nx, ny, realization, plane_before_source, source_redshift, plane_shift, convergence_direct, preload_planes, number_of_planes, number_of_plane_realizations, first_sim_ic, sim_comm, plane_storage_window);  // Load potential values from plane FITS file into potential_array containing these values for calculation.
	// The above function is written such, that is will automatically load a terminal plane (last plane before source plane) ending at source redshift, when that plane number is called.
	// If want to create WL maps for multiple redshifts in one calculation, will need to modify the above.


   // Testing: Replace potential plane with something known:
   // (Comment out for actual result run!)
   /*
   int ii, jj, kk, iii, jjj;
   if (ThisTask==0)
     {
       for (kk=0; kk<parameters.nxny; kk++)
	 {
	   if (plane_number==3 || plane_number==15)
	     {
	       ii=kk/parameters.nx;
	       jj=kk%parameters.nx;
	       iii=ii-parameters.ny/2;
	       jjj=jj-parameters.nx/2;
	       // if (plane_number==3) potential_array[kk]=-0.00001*cos(2*3.141*sqrt(iii*iii+jjj*jjj)/8.0);

	       if (plane_number==3) potential_array[kk]=-0.00001*96.0*cos(2*3.141*sqrt(iii*iii+jjj*jjj)/96.0);
	       // iii+=parameters.ny/16;
	       // jjj+=parameters.nx/16;
	       if (plane_number==15) potential_array[kk]=-0.00001*50.0/6.0*cos(2*3.141*sqrt(iii*iii+jjj*jjj)/50.0);
	       // potential_array[kk]=-0.00001*((ii-parameters.ny/2)*(ii-parameters.ny/2)+(jj-parameters.nx/2)*(jj-parameters.nx/2));
	     }
	   else potential_array[kk]=0.0;
	 }
     }
   */
   // End of Testing.
   ///////////////////////////////////////////////////////

	
   // Now communicate the read in plane to the other subprocesses for this map:
   // Note: each process needs the whole plane, because the light rays could be anywhere on it.
   MPI_Bcast(potential_array, nxny, MPI_DOUBLE, 0, map_comm);

     
     
     // Computation of values for Plane plane_number includes planes up to Plane plane_number-1.

     
#ifdef HAVE_OpenMP

// #pragma omp parallel num_threads(NR_OF_OMP_THREADS) private(pixelcounter, Ybin, Xbin, delPhi1, delPhi2, x, y, U11, U12, U21, U22, temp, chi)
#pragma omp parallel private(pixelcounter, Ybin, Xbin, delPhi1, delPhi2, x, y, U11, U12, U21, U22, temp, chi)
     // Split will be based on Ybin, order is going to be as fast as possibly (one at a time)
     
     { // Begin OpenMP parallelization:
         
     int tid, nthreads;
     tid = omp_get_thread_num();
     if (superrank == 0 && tid==0)
     {
         nthreads = omp_get_num_threads();
         printf("Number of OpenMP threads = %d\n", nthreads);
     }
     
#pragma omp for schedule(dynamic,1)
#endif
	// for (Ybin=0; Ybin<NbinsY; Ybin++) {
	// for (Xbin=0; Xbin<NbinsX; Xbin++) {
    for (pixelcounter=0; pixelcounter<NbinsX*NbinsY; pixelcounter++)
    {
            // Determine pixel coordinates in 2D arrays (instead of two separate for-loops over Ybin and Xbin):
            Ybin=pixelcounter/NbinsX;
            Xbin=pixelcounter%NbinsX;
            
	  
        if (SC[Ybin][Xbin]+0.001>Plane[plane_number-1].comoving_distance) // If ray is already beyond the source galaxy (i.e. plane comoving distance is larger than source galaxy comoving distance, don't propagate deflection angles and lensing quantities further. This if statement also eleminates immediately empty "ray grid pixels" when working with a galaxy catalogue, because these have SC[Ybin][Xbin] (source comoving distance) set to zero and fail this condition on Plane 1. The small incremental value is intended to protect against roundoff errors if galaxies in catalogue have been placed precisely at plane distance (which they never should be intentionally, because the catchment area of the plane extends beyond the plane comoving distance, half way to the next.
        {
            
        
	// Advance construction functions of theta function components:	
	x=theta1[Ybin][Xbin]*Plane[plane_number-1].comoving_distance; // transversal comoving position on Plane "plane_number-1". Needed to compute first end second potential derivatives there.
	y=theta2[Ybin][Xbin]*Plane[plane_number-1].comoving_distance;

	// Record transverse positions of light rays on Plane plane_number-1:
	// RAY-PLANES:	ray_pos[Ybin*parameters.NbinsX+Xbin]=x;
	// RAY-PLANES:     ray_pos[parameters.NbinsX*parameters.NbinsY+Ybin*parameters.NbinsX+Xbin]=y;

            
	if (convergence_direct!=0)
	  {
	    delPhi1=0.0; delPhi2=0.0;
	  }
	else // regular case:
	  {
	    delPhi1=DPhi_2D_multi(potential_array, x, y, plane_number-1, 0, nx, ny, Plane[plane_number-1].binwidth_x, Plane[plane_number-1].binwidth_y, ray_tracing, plane_assignment_mode, raypoint_averaging); // arguments: comoving transversal postitions x and y, derivative direction. Evaluated at Plane "plane_number-1".
	    delPhi2=DPhi_2D_multi(potential_array, x, y, plane_number-1, 1, nx, ny, Plane[plane_number-1].binwidth_x, Plane[plane_number-1].binwidth_y, ray_tracing, plane_assignment_mode, raypoint_averaging);
	  }

	f1[Ybin][Xbin]+=delPhi1/Plane[plane_number-1].scale_factor; // *parameters.h; // parameters.h compensates for fact that derivative in delPhi is in units of Mpc/h.
	g1[Ybin][Xbin]+=Plane[plane_number-1].comoving_distance*delPhi1/Plane[plane_number-1].scale_factor; // *parameters.h;
	f2[Ybin][Xbin]+=delPhi2/Plane[plane_number-1].scale_factor; // *parameters.h;
	g2[Ybin][Xbin]+=Plane[plane_number-1].comoving_distance*delPhi2/Plane[plane_number-1].scale_factor; // *parameters.h;

	
	// Advance construction matrices of distortion tensor components:

	if (convergence_direct!=0)
	  {

	    // Put density contrast from density plane on diagonal:
	    U11=convergence_prefactor*Phi_2D_multi(potential_array, x, y, plane_number-1, nx, ny, Plane[plane_number-1].binwidth_x, Plane[plane_number-1].binwidth_y, ray_tracing, plane_assignment_mode); // evaluated at Plane "plane_number-1".
        U12=0.0;
        U21=U12;
        U22=convergence_prefactor*Phi_2D_multi(potential_array, x, y, plane_number-1, nx, ny, Plane[plane_number-1].binwidth_x, Plane[plane_number-1].binwidth_y, ray_tracing, plane_assignment_mode);

	    // Make old distortion tensor array identity matrix (don't include "feedback" term in equation):
	    A11[Ybin][Xbin]=1.0; A12[Ybin][Xbin]=0.0; A21[Ybin][Xbin]=0.0; A22[Ybin][Xbin]=1.0;
	  }
	else // regular case:
	  {
	    // Calculate Derivative Matrix U at ray position:
	    U11=DDPhi_2D_multi(potential_array, x, y, plane_number-1, 0, nx, ny, Plane[plane_number-1].binwidth_x, Plane[plane_number-1].binwidth_y, ray_tracing, plane_assignment_mode, raypoint_averaging); // evaluated at Plane "plane_number-1".
	    U12=DDPhi_2D_multi(potential_array, x, y, plane_number-1, 1, nx, ny, Plane[plane_number-1].binwidth_x, Plane[plane_number-1].binwidth_y, ray_tracing, plane_assignment_mode, raypoint_averaging);
	    U21=U12;
	    U22=DDPhi_2D_multi(potential_array, x, y, plane_number-1, 2, nx, ny, Plane[plane_number-1].binwidth_x, Plane[plane_number-1].binwidth_y, ray_tracing, plane_assignment_mode, raypoint_averaging);
	  }

	// Advance F and G (calculating F_{n+1} and G_{n+1} from F_n and G_n, where n+1=plane_number):
	chi=Plane[plane_number-1].comoving_distance;
	a=Plane[plane_number-1].scale_factor;
	temp=chi/a*(U11*A11[Ybin][Xbin]+U12*A21[Ybin][Xbin]);
	// temp*=parameters.h; // Puts in h back (compensation for comoving distances and spatial derivatives being in units of Mpc/h) - neutralises h-factors, convergence is h-independent.
	F11[Ybin][Xbin]+=temp;
	G11[Ybin][Xbin]+=chi*temp;
	temp=chi/a*(U11*A12[Ybin][Xbin]+U12*A22[Ybin][Xbin]);
	// temp*=parameters.h; // See above.
	F12[Ybin][Xbin]+=temp;
	G12[Ybin][Xbin]+=chi*temp;
	temp=chi/a*(U21*A11[Ybin][Xbin]+U22*A21[Ybin][Xbin]);
	// temp*=parameters.h;
	F21[Ybin][Xbin]+=temp;
	G21[Ybin][Xbin]+=chi*temp;
	temp=chi/a*(U21*A12[Ybin][Xbin]+U22*A22[Ybin][Xbin]);
	// temp*=parameters.h;
	F22[Ybin][Xbin]+=temp;
	G22[Ybin][Xbin]+=chi*temp;

	// Now advance the actual angles theta1, theta2 and distortion matrix A (these are the final quantities of the step):

        /*
	if (plane_number==parameters.plane_before_source+1) // if plane is the source plane, and calculation is ending, use comoving distance of source plane for chi:
	  {
	    chi=parameters.source_comoving_distance;
	  }
	else // normal case, regular plane in stack, take comoving distance of target to be next such plane:
	  {
	    chi=Plane[plane_number].comoving_distance; // this is chi_{n+1} for A_{n+1}
	  }
         */

     

        
        if (SC[Ybin][Xbin]+0.001<Plane[plane_number].comoving_distance)
        {
            
            if (galaxy_catalogue_type==0) chi=SC[Ybin][Xbin];
            else chi=(Plane[plane_number-1].comoving_distance+Plane[plane_number].comoving_distance)/2.0; // for real galaxy catalogue, snap galaxy positions to plane distance grid instead of using actual source galaxy redshifts here (but using chi=SC[Ybin][Xbin] would also be an option here).
            
            // The first of the two if conditions implies this is the last plane closer to the observer than the galaxy, and therefore the check higher above whether to continue the calculation for this ray will fail, making this the last plane. Thus, use the actual source galaxy redshift for propagation if doing the last plane for this light ray. This is very important, or the algorithm would not compute the proper convergence/shear equation. Do not use the comoving distance of the last plane here, even if the source galaxy comoving distance has been "snapped to plane spacing", because the catchment area of a plane extends beyond the plane comoving distance half way to the next, so the comoving distance of the far catchemt boundary of the last lens plane (or the source galaxy comoving distance itself, stored in the SC-array) is the most consisten value for the calculation.

            last_plane_number=plane_number-1; // number of last plane, used to put into header of WL maps FITS file. For maps this corresponds to the last lens plane still used (there are no terminal planes anymore), for galaxy catalogues this corresponds to the last plane used for the farthest galaxy (but this is currently not recorded in the output.
        }
        else chi=Plane[plane_number].comoving_distance; // use the plane comoving distance for all other light rays.
        
        
	if (ray_tracing!=0) // change angle only if ray-tracing activated.
	  {
	    // Advance theta1 and theta2:
	    theta1[Ybin][Xbin]=theta1_ini[Ybin][Xbin]-f1[Ybin][Xbin]+g1[Ybin][Xbin]/chi; // new angle for Plane "plane_number"
	    theta2[Ybin][Xbin]=theta2_ini[Ybin][Xbin]-f2[Ybin][Xbin]+g2[Ybin][Xbin]/chi;
	  }

	// Advance A: //////////////////////

	
	// Save old values before advancing:
	//PrevA11[Xbin][Ybin]=A11[Xbin][Ybin];
	//PrevA12[Xbin][Ybin]=A12[Xbin][Ybin];
	//PrevA21[Xbin][Ybin]=A21[Xbin][Ybin];
	//PrevA22[Xbin][Ybin]=A22[Xbin][Ybin];
	

	// Advance:
	A11[Ybin][Xbin]=IdentityM11-F11[Ybin][Xbin]+G11[Ybin][Xbin]/chi;
	A12[Ybin][Xbin]=IdentityM12-F12[Ybin][Xbin]+G12[Ybin][Xbin]/chi;
	A21[Ybin][Xbin]=IdentityM21-F21[Ybin][Xbin]+G21[Ybin][Xbin]/chi;
	A22[Ybin][Xbin]=IdentityM22-F22[Ybin][Xbin]+G22[Ybin][Xbin]/chi;

	///////////////////////////////////
        
    } // end decision to propagate if.

	// Calculate line-of-sight investigation values:
	// (Work only for distortion tensor, not deflection angle, so do not read out angles for these in save_WL_maps().)

	
	// Moving sources further and further out:
        //PrevA11[Xbin][Ybin]=A11[Xbin][Ybin]-PrevA11[Xbin][Ybin];
        //PrevA12[Xbin][Ybin]=A12[Xbin][Ybin]-PrevA12[Xbin][Ybin];
        //PrevA21[Xbin][Ybin]=A21[Xbin][Ybin]-PrevA21[Xbin][Ybin];
        //PrevA22[Xbin][Ybin]=A22[Xbin][Ybin]-PrevA22[Xbin][Ybin];
	

	// Recording contribution of only current plane (with plane number  plane_number):
	//SingleA11[Xbin][Ybin]=IdentityM.m11-(Plane[plane_number-1].comoving_distance*(parameters.source_comoving_distance-Plane[plane_number-1].comoving_distance)/(Plane[plane_number-1].scale_factor*parameters.source_comoving_distance))*U.m11;
	//SingleA12[Xbin][Ybin]=IdentityM.m12-(Plane[plane_number-1].comoving_distance*(parameters.source_comoving_distance-Plane[plane_number-1].comoving_distance)/(Plane[plane_number-1].scale_factor*parameters.source_comoving_distance))*U.m12;
	//SingleA21[Xbin][Ybin]=IdentityM.m21-(Plane[plane_number-1].comoving_distance*(parameters.source_comoving_distance-Plane[plane_number-1].comoving_distance)/(Plane[plane_number-1].scale_factor*parameters.source_comoving_distance))*U.m21;
        //SingleA22[Xbin][Ybin]=IdentityM.m22-(Plane[plane_number-1].comoving_distance*(parameters.source_comoving_distance-Plane[plane_number-1].comoving_distance)/(Plane[plane_number-1].scale_factor*parameters.source_comoving_distance))*U.m22;

	
	// Recording all contributions of planes closer and equal to plane_number-1:
        //CumulA11[Xbin][Ybin]=IdentityM.m11-F11[Xbin][Ybin]+G11[Xbin][Ybin]/parameters.source_comoving_distance;
        //CumulA12[Xbin][Ybin]=IdentityM.m12-F12[Xbin][Ybin]+G12[Xbin][Ybin]/parameters.source_comoving_distance;
        //CumulA21[Xbin][Ybin]=IdentityM.m21-F21[Xbin][Ybin]+G21[Xbin][Ybin]/parameters.source_comoving_distance;
        //CumulA22[Xbin][Ybin]=IdentityM.m22-F22[Xbin][Ybin]+G22[Xbin][Ybin]/parameters.source_comoving_distance;
	// This is similar to regular advance, except for current plane (denoted by plane_number) is placed as sources.
	

	//////////////////
	
    } // end pixelcounter loop.
            
	// } // end Xbin loop.

	//	printf("Ybin: %d\n", Ybin); // run progress readout for testing.
	//      fflush(stdout);

	// } // end Ybin loop.
         
#ifdef HAVE_OpenMP
    } // end of OpenMP parallelization.
#endif
     
     
	// Write out difference to previous step for along-the-line-of-sight study:
	if (realization<=0)
	  {
	    // save_WL_results(theta1, theta2, theta1_ini, theta2_ini, PrevA11, PrevA12, PrevA21, PrevA22, writeout_array,plane_number-1, 1);
	    // save_WL_results(theta1, theta2, theta1_ini, theta2_ini, SingleA11, SingleA12, SingleA21, SingleA22, writeout_array,plane_number-1, 2);
	    // save_WL_results(theta1, theta2, theta1_ini, theta2_ini, CumulA11, CumulA12, CumulA21, CumulA22, writeout_array,plane_number-1, 3);
	    // RAY-PLANES:	    save_ray_pos(ray_pos, plane_number-1);
	  }

	// Now we have theta1, theta2, A (distortion tensor matrix) calculated for Plane plane_number.
	// Can either calculate derived quantities like convergence and save results to file, or proceed directly to next plane without saving interim result.
     
     /* // These lines are moved out below, because in new version of Inspector Gadget there are no terminal planes anymore for arbitrary source redshifts; instead allowed source redshifts are restricted to planes located precisely midway between two potential planes.
	// if (plane_number==parameters.plane_before_source+1) // if plane_number points to source plane, we're done, need to save WL map and exit the planes loop:
	  {
	    free_Vector0(potential_array); // with current version of code not needed anymore, so free before allocating writeout arrays.
	    save_WL_results(theta1, theta2, theta1_ini, theta2_ini, A11, A12, A21, A22, writeout_array, plane_number-1, 0); // call with plane_number-1, because that's last FITS potential plane, plane_number is currently pointing at the source plane, which is not a potential plane and file. 
	    plane_number=planes; // With the current version of the code, want to calculate only to the (only) source plane. By setting the plane_number larger than planes, the for loop of the iteration over planes terminates.
	  }
      */

   
} // end plane number loop.
 
    
    free_Matrix0(SC, NbinsY);
    
    free_Vector0(potential_array); // with current version of code not needed anymore, so free before allocating writeout arrays.
    if (galaxy_catalogue_type==0) save_WL_results(theta1, theta2, theta1_ini, theta2_ini, A11, A12, A21, A22, writeout_array, last_plane_number, 0, NbinsX_global, NbinsY_global, nx, realization); // save WL maps.
    // else save_WL_galaxy_catalogue(galaxy_input_filename, galaxy_output_filename, theta1, theta2, theta1_ini, theta2_ini, A11, A12, A21, A22, writeout_array, NbinsX, NbinsY);
    else save_WL_galaxy_catalogue_output_only(galaxy_output_filename, theta1, theta2, theta1_ini, theta2_ini, A11, A12, A21, A22, writeout_array,last_plane_number, NbinsX, NbinsY, number_of_galaxies, number_of_redshifts_per_galaxy); // save WL values at positions of galaxies in catalogue (and don't write galaxy positions into the file to save disk space).
    // NOTE: The galaxy catalogue does not work with the MPI subprocess split (but it does work with OpenMP).



 	 // Free arrays (matrices and array-matrices):
    
    free_Matrix0(theta1_ini, NbinsY);
    free_Matrix0(theta2_ini, NbinsY);
	free_Matrix0(theta1, NbinsY);
	free_Matrix0(theta2, NbinsY);
	free_Matrix0(f1, NbinsY);
	free_Matrix0(f2, NbinsY);
	free_Matrix0(g1, NbinsY);
	free_Matrix0(g2, NbinsY);

	free_Matrix0(A11, NbinsY);
	free_Matrix0(A12, NbinsY);
	free_Matrix0(A21, NbinsY);
	free_Matrix0(A22, NbinsY);
	free_Matrix0(F11, NbinsY);
	free_Matrix0(F12, NbinsY);
	free_Matrix0(F21, NbinsY);
	free_Matrix0(F22, NbinsY);
	free_Matrix0(G11, NbinsY);
	free_Matrix0(G12, NbinsY);
	free_Matrix0(G21, NbinsY);
	free_Matrix0(G22, NbinsY);

	//	free_Vector0(ray_pos);

	/*
        free_Matrix0(PrevA11, parameters.NbinsX);
        free_Matrix0(PrevA12, parameters.NbinsX);
        free_Matrix0(PrevA21, parameters.NbinsX);
        free_Matrix0(PrevA22, parameters.NbinsX);
	
        free_Matrix0(SingleA11, parameters.NbinsX);
        free_Matrix0(SingleA12, parameters.NbinsX);
        free_Matrix0(SingleA21, parameters.NbinsX);
        free_Matrix0(SingleA22, parameters.NbinsX);
	
        free_Matrix0(CumulA11, parameters.NbinsX);
        free_Matrix0(CumulA12, parameters.NbinsX);
        free_Matrix0(CumulA21, parameters.NbinsX);
        free_Matrix0(CumulA22, parameters.NbinsX);
	*/

	//	if (ThisTask==0) free_Vector0(writeout_array);

	printf("Weak lensing arrays freed. Returning from weak lensing routine.\n");
	fflush(stdout);

}

void save_ray_pos(double *ray_pos, int plane_number)
{
  char filename[200];
  long mapdims[3]={parameters.NbinsX, parameters.NbinsY, 2};
  int redshift_tag;

  parameters.first_plane=0;
  parameters.last_plane=plane_number;
  parameters.source_scale_factor=1.0/(parameters.source_redshift+1.0);


  if (plane_number==parameters.plane_before_source) // if last plane before source:
    {
      redshift_tag=((int) floor(100.0*parameters.source_redshift));

      sprintf(filename, "%s/ray_%s_%04dxy_%04dr_%04dp_%04dz.%s", parameters.map_output_path, parameters.simulation_codename, parameters.NbinsX, parameters.realization, plane_number, redshift_tag, parameters.extension);
    }
  else
    {
      sprintf(filename, "%s/ray_%s_%04dxy_%04dr_%04dp.%s", parameters.map_output_path, parameters.simulation_codename, parameters.NbinsX, parameters.realization, plane_number, parameters.extension);
    }

#ifdef HAVE_FITS
  writeRayPlaneFITSimage_f(filename, 3, mapdims, ray_pos, plane_number);
#else
    exit(1); // currently not implemented.
#endif
  printf("Done writing ray plane for Plane %d.\n", plane_number);
  fflush(stdout);

}

void save_WL_results(double **theta1, double **theta2, double **theta1_ini, double **theta2_ini, double **A11, double **A12, double **A21, double **A22, double *writeout_array, int plane_number, int write_mode, int NbinsX_global, int NbinsY_global, int nx, int realization)
{

	printf("Commencing weak lensing subroutine (distortion tensor calculation)...\n");

	int fits_filetype;
	int Xbin, Ybin;
	
	char basename[200]; // local variable temporarily containing basename for type of FITS file that is written (during loop over fits file types).
	char map_file[200]; // local variable containing the full name and path of FITS file that is written.
	long mapdims[2]={NbinsX_global, NbinsY_global}; // variable needed for FITS writing routine.

	int NbinsX, NbinsY;
	int ThisTask, number_of_subprocesses;
	// Adjust for split between subprocesses:
	MPI_Comm_rank(map_comm, &ThisTask);
	MPI_Comm_size(map_comm, &number_of_subprocesses);
	NbinsX=NbinsX_global;
	NbinsY=NbinsY_global/number_of_subprocesses;
	if (NbinsY_global%number_of_subprocesses!=0)
	  {
	    printf("ERROR: Number of Pixels of WL map in vertical dimension (NbinsY) must be a multiple of number_of_subprocesses. Aborting.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
	    exit(1);
	  }

	double *send_array;
	if (ThisTask==0) writeout_array=Vector0(NbinsX_global*NbinsY_global);
	send_array=Vector0(NbinsX*NbinsY);



	
	//struct matrix2x2 Lensing_Matrix;
	
	double gamma1, gamma2, deflection1, deflection2;


	//printf("Done allocating distortion tensor analyzers.\n");
	//timing();
	

	// ************************************************************ //
	// Get proper writeout into header of fits file (proper listing of used planes):
	parameters.first_plane=0;
	parameters.last_plane=plane_number;
	parameters.source_scale_factor=1.0/(parameters.source_redshift+1.0);
	// source_redshift and source_comoving_distance are supplied directly by the Condor script now. 
	// ************************************************************ //


	printf("Writing weak lensing maps into FITS files...\n");
	timing();

	// Exist up to 10 file types below, but to save output space, restrict to first three: convergence and the two shear components:
	for (fits_filetype=1; fits_filetype<=10; fits_filetype++) // caluclate all file types separately and write them into fits file. This way need only one readout array.
	{

		printf("Writing fits_filetype %d ... \n", fits_filetype);
		//timing();

		
		for (Ybin=0; Ybin<NbinsY; Ybin++) {
		for (Xbin=0; Xbin<NbinsX; Xbin++) {	

			/////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////
			// For implementation using standard 2D ray-tracing algorithm with 2D lens planes: //
						
			// kappa:
			if (fits_filetype==1) send_array[Ybin*NbinsX+Xbin]=-(A11[Ybin][Xbin]+A22[Ybin][Xbin])/2.0+1.0;
			// gamma1:
			else if (fits_filetype==2) send_array[Ybin*NbinsX+Xbin]=-(A11[Ybin][Xbin]-A22[Ybin][Xbin])/2.0;
			// gamma2:
			else if (fits_filetype==3) send_array[Ybin*NbinsX+Xbin]=-(A12[Ybin][Xbin]+A21[Ybin][Xbin])/2.0;
			// omega:
			else if (fits_filetype==4) send_array[Ybin*NbinsX+Xbin]=-(A12[Ybin][Xbin]-A21[Ybin][Xbin])/2.0;
			// gamma_modulus:
			else if (fits_filetype==5)
			{
				gamma1=-(A11[Ybin][Xbin]-A22[Ybin][Xbin])/2.0;
				gamma2=-(A12[Ybin][Xbin]+A21[Ybin][Xbin])/2.0;
				send_array[Ybin*NbinsX+Xbin]=sqrt(gamma1*gamma1+gamma2*gamma2);
			}
			// gamma_angle:
			else if (fits_filetype==6)
			{
				gamma1=-(A11[Ybin][Xbin]-A22[Ybin][Xbin])/2.0;
				gamma2=-(A12[Ybin][Xbin]+A21[Ybin][Xbin])/2.0;
				send_array[Ybin*NbinsX+Xbin]=atan(gamma2/gamma1);
			}
			// deflection1:
			else if (fits_filetype==7) send_array[Ybin*NbinsX+Xbin]=theta1[Ybin][Xbin]-theta1_ini[Ybin][Xbin];
			// deflection2:
			else if (fits_filetype==8) send_array[Ybin*NbinsX+Xbin]=theta2[Ybin][Xbin]-theta2_ini[Ybin][Xbin];
			// deflection_total:
			else if (fits_filetype==9) 
			  {
			    deflection1=theta1[Ybin][Xbin]-theta1_ini[Ybin][Xbin];
			    deflection2=theta2[Ybin][Xbin]-theta2_ini[Ybin][Xbin];
			    send_array[Ybin*NbinsX+Xbin]=sqrt(deflection1*deflection1+deflection2*deflection2);
			  }
			// deflection_winkel:
			else if (fits_filetype==10)
			  {
                            deflection1=theta1[Ybin][Xbin]-theta1_ini[Ybin][Xbin];
                            deflection2=theta2[Ybin][Xbin]-theta2_ini[Ybin][Xbin];
			    send_array[Ybin*NbinsX+Xbin]=atan(deflection2/deflection1);
			  }
											
		}}
		

		printf("Superrank: %d right before Gather.\n", superrank);
		fflush(stdout);

		MPI_Gather(send_array, NbinsX*NbinsY, MPI_DOUBLE, writeout_array, NbinsX*NbinsY, MPI_DOUBLE, 0, map_comm);

		printf("Superrank: %d right after Gather.\n", superrank);
		fflush(stdout);

		
		if (ThisTask==0) // If lead subprocess, write actual file:
		  {
	
	
		// Set basename properly for writing FITS file:
	
		// kappa:
		if (fits_filetype==1) strncpy(basename, parameters.convergence_basename, MAXNAME); 
		// gamma1:
		else if (fits_filetype==2) strncpy(basename, parameters.shear1_basename, MAXNAME);
		// gamma2:
		else if (fits_filetype==3) strncpy(basename, parameters.shear2_basename, MAXNAME);
		// omega:
		else if (fits_filetype==4) strncpy(basename, parameters.omega_basename, MAXNAME);
		// gamma_modulus:
		else if (fits_filetype==5) strncpy(basename, parameters.shear_modulus_basename, MAXNAME);
		// gamma_angle:
		else if (fits_filetype==6) strncpy(basename, parameters.shear_angle_basename, MAXNAME);
		// deflection1:
		else if (fits_filetype==7) strncpy(basename, parameters.deflection1_basename, MAXNAME);
		// deflection2:
		else if (fits_filetype==8) strncpy(basename, parameters.deflection2_basename, MAXNAME);
		// deflection_total:
		else if (fits_filetype==9) strncpy(basename, parameters.deflection_total_basename, MAXNAME);
		// deflection_winkel:
		else if (fits_filetype==10) strncpy(basename, parameters.deflection_winkel_basename, MAXNAME);

		// Write weak lensing map as FITS image file for file type:
		int redshift_tag;
		redshift_tag=((int) floor(100.0*parameters.source_redshift));
		if (write_mode==0) // Regular weak lensing maps to be analyzed with Anacondor.
		  {
		    // Includes both resolution labels below:
		    // sprintf(map_file, "%s/%s%s_%s_%04dxy_%04nx_%04dr_%04dp_%04dz_og.%s", parameters.map_output_path, parameters.map_basename, basename, parameters.simulation_codename, parameters.NbinsX, parameters.nx, parameters.realization, plane_number, redshift_tag, parameters.extension);
		    // To have xy in file name represent resolution of lattice on lens planes (plane dimension): (This is the preferred option, as plane resolution impacts the power spectrum at high l (see Sato et al. (2009)), so it needs to be tested with different values. The number of light-rays can be kept at 2048x2048 (increasing it would just add higher l's to the power spectrum, not improve its quality).)
                    sprintf(map_file, "%s/%s%s_%s_%04dxy_%04dr_%04dp_%04dz_og.%s", parameters.map_output_path, parameters.map_basename, basename, parameters.simulation_codename, nx, realization, plane_number, redshift_tag, parameters.extension);
		    // To have xy in file name represent resolution of light rays in map (map dinemsion) instead:
		    // sprintf(map_file, "%s/%s%s_%s_%04dxy_%04dr_%04dp_%04dz_og.%s", parameters.map_output_path, parameters.map_basename, basename, parameters.simulation_codename, parameters.NbinsX, parameters.realization, plane_number, redshift_tag, parameters.extension);
		    if (fits_filetype<=3 || fits_filetype==7 || fits_filetype==8)
			{
#ifdef HAVE_FITS				
				writeWLmapFITSimage_f(map_file, 2, mapdims, writeout_array); // writes a FITS image file with signed float pixel values.
                // WARNING: Note that the above function still makes use of a lot of global variables, in particular of the parameters structure to write into the header.
#else
				exit(1); // currently not implemented, must have FITS.
#endif
			}
		  }
		else // Line-of-sight special writeout series:
		  {
		    if (write_mode==1) sprintf(map_file, "%s/LoS-Prev_%s%s_%s_%04dxy_%04dr_%04dp.%s", parameters.map_output_path, parameters.map_basename, basename, parameters.simulation_codename, NbinsX_global, realization, plane_number, parameters.extension);
		    else if (write_mode==2) sprintf(map_file, "%s/LoS-Single_%s%s_%s_%04dxy_%04dr_%04dp.%s", parameters.map_output_path, parameters.map_basename, basename, parameters.simulation_codename, NbinsX_global, realization, plane_number, parameters.extension);
		    else if (write_mode==3) sprintf(map_file, "%s/LoS-Cumul_%s%s_%s_%04dxy_%04dr_%04dp.%s", parameters.map_output_path, parameters.map_basename, basename, parameters.simulation_codename, NbinsX_global, realization, plane_number, parameters.extension);
		    else
		      {
			printf("Write_mode %d not implemented. Aborting\n", write_mode);
            MPI_Abort(MPI_COMM_WORLD, 1001);
			exit(1001);
		      }
		    if (fits_filetype==1)
			{
#ifdef HAVE_FITS
				writeWLmapFITSimage_f(map_file, 2, mapdims, writeout_array); // writes a FITS image file with signed float pixel values.
#else
                exit(1); // currently not implemented, must have FITS.
#endif			 
			}
		  }

		  } // end if clause for ThisTask==0.
		
	} // ending loop over FITS file types that were written.
	
	if (ThisTask==0) free_Vector0(writeout_array);
	free_Vector0(send_array);
	
	printf("Finished writing weak lensing FITS files.\n");
	timing();
	/////////////////////////////////////////////////////////////
	
	
	
    printf("Weak lensing quantities calculated successfully.\n");
    fflush(stdout);
	
			
    return;
}

//Plane[plane_number].binwidth_x
double Phi_2D_multi(double *potential_array, double x, double y, int plane_number, int nx, int ny, double binwidth_x, double binwidth_y, int ray_tracing, int plane_assignment_mode)
     // Name is historical, this actually reads out plane data without taking derivatives, and is used for density planes and convergence, when calculated via approximation.
{
  double Phi;
  int Xbin_A, Ybin_A;

  Xbin_A= ((int) floor(x/binwidth_x+0.5*nx)); // snapshot number refers here to the snapshot box, in which the plane is located, based on its redshift.
  Ybin_A= ((int) floor(y/binwidth_y+CONSIST+0.5*ny));
  // FLIP-CEIL:  Ybin_A= ((int) ceil(-y/Plane[plane_number].binwidth_y+0.5*parameters.ny)); // no recentering here, want to know between which four points it lies.

  if (ray_tracing!=0 && plane_assignment_mode==0)
    {
      // Compensate for light-rays bent out of range of array by completing periodically once (since array periodic anyway):
      if (Xbin_A<0)
	{
	  Xbin_A=Xbin_A+nx;
	  x=x+nx*binwidth_x;
	}
      if (Xbin_A>=nx)
	{
	  Xbin_A=Xbin_A-nx;
	  x=x-nx*binwidth_x;
	}
      if (Ybin_A<0)
	{
	  Ybin_A=Ybin_A+ny;
	  y=y-ny*binwidth_y;
	}
      if (Ybin_A>=ny)
	{
	  Ybin_A=Ybin_A-ny;
	  y=y+ny*binwidth_y;
	}
    }


  if (Xbin_A<0 || Xbin_A>=nx || Ybin_A<0 || Ybin_A>=ny) // checks that light-ray terminated properly in array after one periodic correction to avoid violation of array boundaries (second corrction not allowed, because it would mean that the light ray wrapped around the whole survey width once.
    {
      printf("WARNING: Light ray out of range of DPhi_array at Plane %d. Requested bin values are (%d, %d). Unable to calculate first derivative of potential for this ray. Setting to zero.\n", plane_number, Xbin_A, Ybin_A);
      return 0;
    }

  // TSC raypoint averaging done in any case:
  Phi=raypoint_average_multi(potential_array, x, y, Xbin_A, Ybin_A, plane_number, 0 , 0, nx, ny, binwidth_x, binwidth_y); // 0 in last argument extracts zeroth derivative (function value without derivative) when calling raypoint_average. Sixth argument (derivative_type), is zero here because not needed anyway.

  return Phi;

}

// Plane[plane_number].binwidth_x
double DPhi_2D_multi(double *potential_array, double x, double y, int plane_number, int derivative_type, int nx, int ny, double binwidth_x, double binwidth_y, int ray_tracing, int plane_assignment_mode, int raypoint_averaging)
// returns second derivative of potential, i.e. FFTW-backward of grid get_DDPhi-ft, and then initialize 2D spline over that 2D-DDPhi field, to be evaluated by add_line_of_sight readout.
{
	double DPhi;
	
	// Using no interpolation or smoothening, whatever Fourier grid box line-of-sight hits:
	// (More sophisticated would be to average over all Fourier boxes that fall into larger bin box of weak lensing analysis (Xbin Ybin in weak_lensing.c, NbinX instead of nx, etc.))
	int Xbin_A, Ybin_A;
	
    // on the fly derivative calculation instead of via pre-calculated Plane[plane_number].DPhi_array:
		
		double DPhi_1, DPhi_2, DPhi_3, DPhi_4;
		int Xbin_Aplus, Ybin_Aplus;
		double Xbin_A_phys, Ybin_A_phys; // used for bilinear and bicubic interpolation.
		double Xbin_A_phys_plus, Ybin_A_phys_plus; // used only for bicubic interpolation.
		double t, u;


		
		Xbin_A= ((int) floor(x/binwidth_x+0.5*nx)); // snapshot number refers here to the snapshot box, in which the plane is located, based on its redshift.
		Ybin_A= ((int) floor(y/binwidth_y+CONSIST+0.5*ny));
		// FLIP-CEIL: Ybin_A= ((int) ceil(-y/Plane[plane_number].binwidth_y+0.5*parameters.ny)); // no recentering here, want to know between which four points it lies.
		
		if (ray_tracing!=0 && plane_assignment_mode==0)
		{
			// Compensate for light-rays bent out of range of array by completing periodically once (since array periodic anyway):
			if (Xbin_A<0)
			{
				Xbin_A=Xbin_A+nx;
				x=x+nx*binwidth_x;
			}
			if (Xbin_A>=nx)
			{
				Xbin_A=Xbin_A-nx;
				x=x-nx*binwidth_x;
			}
			if (Ybin_A<0)
			{
				Ybin_A=Ybin_A+ny;
				y=y+ny*binwidth_y;
			}	
			if (Ybin_A>=ny)
			{
				Ybin_A=Ybin_A-ny;
				y=y-ny*binwidth_y;
			}
		}
	
	
		if (Xbin_A<0 || Xbin_A>=nx || Ybin_A<0 || Ybin_A>=ny) // checks that light-ray terminated properly in array after one periodic correction to avoid violation of array boundaries (second corrction not allowed, because it would mean that the light ray wrapped around the whole survey width once.
		{
			printf("WARNING: Light ray out of range of DPhi_array at Plane %d. Requested bin values are (%d, %d). Unable to calculate first derivative of potential for this ray. Setting to zero.\n", plane_number, Xbin_A, Ybin_A);
			return 0;
		}
		
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Now various schemes of averaging of derivatives at gridpoints for value at where ray pierces the potential plane: 
		if (raypoint_averaging==1) // linear ray averaging:
		{
			if (Xbin_A>=nx-1) Xbin_Aplus=0;
			else Xbin_Aplus=Xbin_A+1;
			// if (Ybin_A<=0) Ybin_Aminus=parameters.ny-1;
			// else Ybin_Aminus=Ybin_A-1;
			if (Ybin_A>=ny-1) Ybin_Aplus=0;
			else Ybin_Aplus=Ybin_A+1;
		
			DPhi_1=FD_DPhi_multi(potential_array, Xbin_A, Ybin_A, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
			DPhi_2=FD_DPhi_multi(potential_array, Xbin_Aplus, Ybin_A, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
			DPhi_3=FD_DPhi_multi(potential_array, Xbin_Aplus, Ybin_Aplus, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
			DPhi_4=FD_DPhi_multi(potential_array, Xbin_A, Ybin_Aplus, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
		
			Xbin_A_phys=(Xbin_A-0.5*nx)*binwidth_x;
			Ybin_A_phys=(Ybin_A-0.5*ny)*binwidth_y;
			// FLIP-CEIL: Ybin_A_phys=-(Ybin_A-0.5*parameters.ny)*Plane[plane_number].binwidth_y;
		
			t=(x-Xbin_A_phys)/binwidth_x; // equation for linear interpolation from Numerical Recipes in C, Section 3.6.
			u=(y-Ybin_A_phys)/binwidth_y; // equation for linear interpolation from Numerical Recipes in C, Section 3.6.
		
			if (t<0 || u<0)
			{
				printf("ERROR: t or u is smaller than 0. t, u: %e %e \n.", t, u); 
			}
		
			DPhi=(1-t)*(1-u)*DPhi_1+t*(1-u)*DPhi_2+t*u*DPhi_3+(1-t)*u*DPhi_4; // main equation for linear interpolation from Numerical Recipes in C, Section 3.6.
			// Turn off interpolation (do not use for actual result runs!):
			// DPhi=DPhi_1; // test line without bilinear interpolation, corresponds to nearest lower-left grid point (rounded down), and corresponds to test line in bicubic interpolation case.
		}
		else if (raypoint_averaging==2) // bicubic averaging of ray position on grid points:
		  {


            if (Xbin_A>=nx-1) Xbin_Aplus=0;
            else Xbin_Aplus=Xbin_A+1;
		    if (Ybin_A>=ny-1) Ybin_Aplus=0;
		    else Ybin_Aplus=Ybin_A+1;

		    Xbin_A_phys=(Xbin_A-0.5*nx)*binwidth_x;
		    Ybin_A_phys=(Ybin_A-0.5*ny)*binwidth_y;

		    Xbin_A_phys_plus=(Xbin_A+1-0.5*nx)*binwidth_x;
		    Ybin_A_phys_plus=(Ybin_A+1-0.5*ny)*binwidth_y;

		    int ii, a, b;
		    double y0[5], y1[5], y2[5], y12[5];
		    // CAUTION: Usage of above arrays y0, y1, y2, y12 runs from 1,...,4
		    // (Numerical Recipes in C convention to start array at index 1 instead of at 0 as is standard in C).
		    double ansy, ansy1, ansy2;
		    ansy=0.0;
		    ansy1=0.0;
		    ansy2=0.0;


                    // Field values to be interpolated:                            
                    // Read values (second derivatives of potential), first derivatives and cross derivatives of the values at the four corners of grid point box around actual point to be interpolated:       
		    for (ii=1; ii<=4; ii++)
		      {
			if (ii==1) // lower left grid point                        
			  {
			    a=Xbin_A;
			    b=Ybin_A;
			  }
			else if (ii==2) // lower right grid point                  
			  {
			    a=Xbin_A+1;
			    b=Ybin_A;
			  }
			else if (ii==3) // upper right grid point                  
			  {
			    a=Xbin_A+1;
			    b=Ybin_A+1;
			  }
			else if (ii==4) // upper left grid point                   
			  {
			    a=Xbin_A;
			    b=Ybin_A+1;
			  }
			else // non-existent choice (error for safety)             
			  {
			    printf("ERROR: In bicubic interpolation, this case does not exist.\n");
			    exit(1);
			  }

                        // Values of second derivative of potential at corners of squares:
                        y0[ii]=FD_DPhi_multi(potential_array, a, b, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
                        // Values of first derivatives (used for bicubic interpolation only) of first derivative of potential:
                        y1[ii]=(FD_DPhi_multi(potential_array, a+1, b, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)-FD_DPhi_multi(potential_array, a-1, b, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y))/(2.0*binwidth_x);
                        y2[ii]=(FD_DPhi_multi(potential_array, a, b+1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)-FD_DPhi_multi(potential_array, a, b-1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y))/(2.0*binwidth_y);
                        // Values of cross-derivatives (used for bicubic interpolation only) of first derivatives of potential:
			y12[ii]=(FD_DPhi_multi(potential_array, a+1, b+1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)-FD_DPhi_multi(potential_array, a+1, b-1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)-FD_DPhi_multi(potential_array, a-1, b+1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)+FD_DPhi_multi(potential_array, a-1, b-1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y))/(4.0*binwidth_x*binwidth_y);

			// Some of the above can probably be coded torun faster by calculating the second (and third) derivative of the potential directly in some cases, but the above makes for a simpler code, since the first derivative used for the bicubic interpolation may be in a different direction than the other first derivative specified by derivative_type.

                      }

                    //              printf("About to do bicubic interpolation of light ray positions...\n\");
                    // fflush(stdout);                                             

		    // Execute bicubic interpolation routine from Numerical Recipes in C:
		    bcuint(y0, y1, y2, y12, Xbin_A_phys, Xbin_A_phys_plus, Ybin_A_phys, Ybin_A_phys_plus, x, y, &ansy, &ansy1, &ansy2);

		    // printf("Done bicubic interpolation.\n");                    
		    // fflush(stdout);                                             

                    DPhi=ansy;
		    // Turn off interpolation (do not use for actual result runs!): 
                    // DPhi=y0[1]; // test line without bicubic interpolation, corresponds to nearest lower-left grid point (rounded down), and corresponds to test line in bilinear interpolation case.

		  }
		else if (raypoint_averaging==3) // TSC averaging of ray position on grid points:
		{
			DPhi=raypoint_average_multi(potential_array, x, y, Xbin_A, Ybin_A, plane_number, derivative_type, 1, nx, ny, binwidth_x, binwidth_y); // 1 in last argument extracts first derivative when calling raypoint_average.
		}
		else
		  {
		    printf("ERROR: Selection of raypoint_averaging is not valid. Please choose one of 1: bilinear, 2: bicubic, 3: TSC (the latter contains slight smoothing over nearest neighbors). Bicubic is preferred and the most accurate selection.\n");
		    exit(1);
		  }
	
	
	return DPhi;
}

// Plane[plane_number].binwidth_x
double DDPhi_2D_multi(double *potential_array, double x, double y, int plane_number, int derivative_type, int nx, int ny, double binwidth_x, double binwidth_y, int ray_tracing, int plane_assignment_mode, int raypoint_averaging)
// returns second derivative of potential, i.e. FFTW-backward of grid get_DDPhi-ft, and then initialize 2D spline over that 2D-DDPhi field, to be evaluated by add_line_of_sight readout.
{
	//int i, j;
	double DDPhi;

	// Using no interpolation or smoothening, whatever Fourier grid box line-of-sight hits:
	// (More sophisticated would be to average over all Fourier boxes that fall into larger bin box of weak lensing analysis (Xbin Ybin in weak_lensing.c, NbinX instead of nx, etc.))
	int Xbin_A, Ybin_A;
	//printf("x_particle, y_particle, binwidth_x: %e %e %e.\n", x_particle, y_particle, binwidth_x);
	
	
	
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// on the fly derivative calculation instead of via pre-calculated Plane[plane_number].DPhi_array:
		
		double DDPhi_1, DDPhi_2, DDPhi_3, DDPhi_4;
		int Xbin_Aplus, Ybin_Aplus;
		double Xbin_A_phys, Ybin_A_phys; // used for bilinear and bicubic interpolation. 
		double Xbin_A_phys_plus, Ybin_A_phys_plus; // used only for bicubic interpolation.
		double t, u;
		
		Xbin_A= ((int) floor(x/binwidth_x+0.5*nx)); // snapshot number refers here to the snapshot box, in which the plane is located, based on its redshift.
		Ybin_A= ((int) floor(y/binwidth_y+CONSIST+0.5*ny));
		// FLIP-CEIL: Ybin_A= ((int) ceil(-y/Plane[plane_number].binwidth_y+0.5*parameters.ny)); // no recentering here, want to know between which four points it lies.
		
		if (ray_tracing!=0)
		{
			// Compensate for light-rays bent out of range of array by completing periodically once (since array periodic anyway):
			if (Xbin_A<0)
			{
				Xbin_A=Xbin_A+nx;
				x=x+nx*binwidth_x;
			}
			if (Xbin_A>=nx)
			{
				Xbin_A=Xbin_A-nx;
				x=x-nx*binwidth_x;
			}
			if (Ybin_A<0)
			{
				Ybin_A=Ybin_A+ny;
				y=y+ny*binwidth_y;
			}	
			if (Ybin_A>=ny)
			{
				Ybin_A=Ybin_A-ny;
				y=y-ny*binwidth_y;
			}
		}
	
		if (Xbin_A<0 || Xbin_A>=nx || Ybin_A<0 || Ybin_A>=ny) // checks that light-ray terminated properly in array after one periodic correction to avoid violation of array boundaries (second corrction not allowed, because it would mean that the light ray wrapped around the whole survey width once.
		{
			printf("WARNING: Light ray out of range of DDPhi_array at Plane %d. Bin requested is (%d, %d). Unable to calculate second derivative of potential for this ray. Setting to zero.\n", plane_number, Xbin_A, Ybin_A);
			printf("Plane, associated Snapshot, Plane.physical_distance, boxcenter.physical_distance %d %d %e \n", plane_number, Plane[plane_number].snapshot, Plane[plane_number].physical_distance);
			return 0;
		}

		// Now various schemes of averaging of derivatives at gridpoints for value at where ray pierces the potential plane: 
		if (raypoint_averaging==1) // linear ray averaging:
		{
			if (Xbin_A>=nx-1) Xbin_Aplus=0;
			else Xbin_Aplus=Xbin_A+1;
			//			if (Ybin_A<=0) Ybin_Aminus=parameters.ny-1;
			// else Ybin_Aminus=Ybin_A-1;

                        if (Ybin_A>=ny-1) Ybin_Aplus=0;
                        else Ybin_Aplus=Ybin_A+1;

		
			DDPhi_1=FD_DDPhi_multi(potential_array, Xbin_A,Ybin_A, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
			DDPhi_2=FD_DDPhi_multi(potential_array, Xbin_Aplus,Ybin_A, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
			DDPhi_3=FD_DDPhi_multi(potential_array, Xbin_Aplus,Ybin_Aplus, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
			DDPhi_4=FD_DDPhi_multi(potential_array, Xbin_A,Ybin_Aplus, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
		
			Xbin_A_phys=(Xbin_A-0.5*nx)*binwidth_x;
			Ybin_A_phys=(Ybin_A-0.5*ny)*binwidth_y;
			// FLIP-CEIL: Ybin_A_phys=-(Ybin_A-0.5*parameters.ny)*Plane[plane_number].binwidth_y;
		
			t=(x-Xbin_A_phys)/binwidth_x; // equation for linear interpolation from Numerical Recipes in C, Section 3.6.
			u=(y-Ybin_A_phys)/binwidth_y; // equation for linear interpolation from Numerical Recipes in C, Section 3.6.
		
			if (t<0 || u<0)
			{
				printf("ERROR: t or u is smaller than 0. t, u: %e %e \n.", t, u); 
			}
		
			DDPhi=(1-t)*(1-u)*DDPhi_1+t*(1-u)*DDPhi_2+t*u*DDPhi_3+(1-t)*u*DDPhi_4; // main equation for linear interpolation from Numerical Recipes in C, Section 3.6.
			// Turn off interpolation (do not use for actual result runs!):
			// DDPhi=DDPhi_1; // test line without bilinear interpolation, corresponds to nearest lower-left grid point (rounded down), and corresponds to test line in bicubic interpolation case. 
		}
		else if (raypoint_averaging==2) // bicubic interpolation.
		  {

		    if (Xbin_A>=nx-1) Xbin_Aplus=0;
		    else Xbin_Aplus=Xbin_A+1;
		    //                      if (Ybin_A<=0) Ybin_Aminus=parameters.ny-1;     
		    // else Ybin_Aminus=Ybin_A-1;                                           

		    if (Ybin_A>=ny-1) Ybin_Aplus=0;
		    else Ybin_Aplus=Ybin_A+1;

		    Xbin_A_phys=(Xbin_A-0.5*nx)*binwidth_x;
		    Ybin_A_phys=(Ybin_A-0.5*ny)*binwidth_y;

                    Xbin_A_phys_plus=(Xbin_A+1-0.5*nx)*binwidth_x;
                    Ybin_A_phys_plus=(Ybin_A+1-0.5*ny)*binwidth_y;

		    int ii, a, b;
		    double y0[5], y1[5], y2[5], y12[5];
                    // CAUTION: Usage of above arrays y0, y1, y2, y12 runs from 1,...,4
                    // (Numerical Recipes in C convention to start array at index 1 instead of at 0 as is standard in C).  
		    double ansy, ansy1, ansy2;
		    ansy=0.0;
		    ansy1=0.0;
		    ansy2=0.0;

		    // Field values to be interpolated:
		    // Read values (second derivatives of potential), first derivatives and cross derivatives of the values at the four corners of grid point box around actual point to be interpolated:
		    for (ii=1; ii<=4; ii++)
		      {
			if (ii==1) // lower left grid point
			  {
			    a=Xbin_A;
			    b=Ybin_A;
			  }
			else if (ii==2) // lower right grid point
			  {
			    a=Xbin_A+1;
			    b=Ybin_A;
			  }
			else if (ii==3) // upper right grid point
			  {
			    a=Xbin_A+1;
			    b=Ybin_A+1;
			  }
			else if (ii==4) // upper left grid point
			  {
			    a=Xbin_A;
			    b=Ybin_A+1;
			  }
			else // non-existent choice (error for safety)
			  {
			    printf("ERROR: In bicubic interpolation, this case does not exist.\n");
			    exit(1);
			  }

			// Values of second derivative of potential at corners of squares:
			y0[ii]=FD_DDPhi_multi(potential_array, a, b, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
			// Values of first derivatives (used for bicubic interpolation only) of second derivative of potential: 
			y1[ii]=(FD_DDPhi_multi(potential_array, a+1, b, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)-FD_DDPhi_multi(potential_array, a-1, b, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y))/(2.0*binwidth_x);
			y2[ii]=(FD_DDPhi_multi(potential_array, a, b+1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)-FD_DDPhi_multi(potential_array, a, b-1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y))/(2.0*binwidth_y);
			// Values of cross-derivatives (used for bicubic interpolation only) of second derivatives of potential:
			y12[ii]=(FD_DDPhi_multi(potential_array, a+1, b+1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)-FD_DDPhi_multi(potential_array, a+1, b-1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)-FD_DDPhi_multi(potential_array, a-1, b+1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y)+FD_DDPhi_multi(potential_array, a-1, b-1, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y))/(4.0*binwidth_x*binwidth_y);

			// Some of the above can probably be coded torun faster by calculating the third (and fourth) derivative of the potential directly in some cases, but this makes for a simpler code, since the first derivative used for bicubic interpolation may be in a different direction than the second derivative specified by derivative_type.

		      }

		    //		    printf("About to do bicubic interpolation for second derivative of potential (weak lensing contribution)...\n");
		    // fflush(stdout);

		    		    bcuint(y0, y1, y2, y12, Xbin_A_phys, Xbin_A_phys_plus, Ybin_A_phys, Ybin_A_phys_plus, x, y, &ansy, &ansy1, &ansy2);

		    // printf("Done bicubic interpolation.\n");
		    // fflush(stdout);

		    DDPhi=ansy;
		    // Turn off interpolation (do not use for actual result runs!): 
		    // DDPhi=y0[1]; // test line without bicubic interpolation, corresponds to nearest lower-left grid point (rounded down), and corresponds to test line in bilinear interpolation case.  

		  }
		else if (raypoint_averaging==3) // TSC averaging of ray position on grid points (warning: this contains a smoothing effect, bicubic interpolation preferred):
		{
			DDPhi=raypoint_average_multi(potential_array, x, y, Xbin_A, Ybin_A, plane_number, derivative_type, 2, nx, ny, binwidth_x, binwidth_y); // 2 in last argument extracts second derivative when calling raypoint_average.
		}
		else
                  {
                    printf("ERROR: Selection of raypoint_averaging is not valid. Please chooseone of 1: bilinear, 2: bicubic, 3: TSC (the latter contains slight smoothing over nearest neighbors). Bicubic is preferred and the most accurate selection.\n");
                    exit(1);
                  }

	
	
	return DDPhi;
}


// Plane[plane_number].binwidth_x
double raypoint_average_multi(double *potential_array, double x, double y, int Xbin_A, int Ybin_A, int plane_number, int derivative_type, int derivative_order, int nx, int ny, double binwidth_x, double binwidth_y)
{
	int i,j;
	int Xbin, Ybin;
	double gridpoint_x, gridpoint_y;
	double derivative, derivative_fraction, total_derivative_fraction, derivative_value;
	derivative=0.0;
	derivative_fraction=0.0;
	total_derivative_fraction=0.0;
	
	int points_check;
	points_check=2;

	for (i=-points_check;i<=points_check;i++)
	{
		for (j=-points_check;j<=points_check;j++)
		{
			Xbin=Xbin_A+j;
			Ybin=Ybin_A+i;
			
			// Determine comoving coordinate of gridpoint, without periodic completion! This serves just as a distance measure for the weight function:
			gridpoint_x=((Xbin)-0.5*nx-0.5)*binwidth_x;
			gridpoint_y=((Ybin-CONSIST)-0.5*ny-0.5)*binwidth_y;
			// FLIP-CEIL: gridpoint_y=(-(Ybin)+0.5*parameters.ny-0.5)*Plane[plane_number].binwidth_y;

			
			// Adjust for periodicity of boundary conditions for actual indices in array:
			if (Xbin<0)
			{
				Xbin+=nx;
			}
			else if (Xbin>=nx)
			{
				Xbin-=nx;
			}
			if (Ybin<0)
			{
				Ybin+=ny;
			}
			else if (Ybin>=ny)
			{
				Ybin-=ny;
			}
			
			// Insert fraction of particle into that gridpoint:
			derivative_fraction=W_TSC(fabs((x-gridpoint_x))/binwidth_x)*W_TSC(fabs((y-gridpoint_y)/binwidth_y));
			if (derivative_order==0) derivative_value=potential_array[Ybin*nx+Xbin]; // direct readout (used for density planes and not necessary for ray-tracing operation mode. Note that for potential_array, indices are swapped.
			else if (derivative_order==1) derivative_value=FD_DPhi_multi(potential_array, Xbin,Ybin, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
			else if (derivative_order==2) derivative_value=FD_DDPhi_multi(potential_array, Xbin,Ybin, plane_number, derivative_type, nx, ny, binwidth_x, binwidth_y);
			else
			{
				printf("This derivative order is not implemented. You chose Order %d.\n", derivative_order);
				exit(1);
			}
			derivative=derivative+derivative_fraction*derivative_value;
			total_derivative_fraction+=derivative_fraction;
			
		}
	}
	
	if (total_derivative_fraction<0.99999 || total_derivative_fraction>1.00001)
	{
		printf("ERROR: TSC particle fraction inserted: %e (instead of =1). Aborting.\n", total_derivative_fraction);
        MPI_Abort(MPI_COMM_WORLD, 1);
		exit(1);
	}
	
	return derivative;

}


// Plane[plane_number].binwidth_x
double FD_DPhi_multi(double *potential_array, int X, int Y, int plane_number, int derivative_type, int nx, int ny, double binwidth_x, double binwidth_y)
{
	double DPhi_FD;
	int i, j, iplus, iminus, jplus, jminus;
	
	// i=Y+1; // Adjust for Numerical Recipes counting convention, in which Plane[plane_number].Phi_array is given.
	// j=X+1; // Adjust for Numerical Recipes counting convention, in which Plane[plane_number].Phi_array is given.
	i=Y;
	j=X;

	if (i==ny-1) iplus=0;
	else iplus = i+1;
		
	if (i==0) iminus=ny-1;
	else iminus = i-1;
	
	if (j==nx-1) jplus=0;
	else jplus = j+1;
			
	if (j==0) jminus=nx-1;
	else jminus = j-1;

	
	if (derivative_type==0)
	// x-derivative:
	DPhi_FD=(potential_array[i*nx+jplus]-potential_array[i*nx+jminus])/(2.0*binwidth_x);
	else if (derivative_type==1)
	// y-derivative:
	DPhi_FD=(potential_array[iplus*nx+j]-potential_array[iminus*nx+j])/(2.0*binwidth_y);
	else
	{
		printf("ERROR: Wrong derivative type supplied.\n");
		exit(1);
	}

	return DPhi_FD;
}

double FD_DDPhi_multi(double *potential_array, int X, int Y, int plane_number, int derivative_type, int nx, int ny, double binwidth_x, double binwidth_y)
{
	double DDPhi_FD;
	int i, j, iplus, iminus, jplus, jminus;
	
	// i=Y+1; // Adjust for Numerical Recipes counting convention, in which Plane[plane_number].Phi_array is given.
	// j=X+1; // Adjust for Numerical Recipes counting convention, in which Plane[plane_number].Phi_array is given.
	i=Y;
	j=X;

	if (i==ny-1) iplus=0;
	else iplus = i+1;
		
	if (i==0) iminus=ny-1;
	else iminus = i-1;
	
	if (j==nx-1) jplus=0;
	else jplus = j+1;
			
	if (j==0) jminus=nx-1;
	else jminus = j-1;

	if (derivative_type==0)
	// xx-derivative:
	DDPhi_FD=(potential_array[i*nx+jplus]+potential_array[i*nx+jminus]-2.0*potential_array[i*nx+j])/(binwidth_x*binwidth_x);
	else if (derivative_type==1)
	// xy-derivative:
	DDPhi_FD=(potential_array[iminus*nx+jminus]+potential_array[iplus*nx+jplus]-potential_array[iplus*nx+jminus]-potential_array[iminus*nx+jplus])/(4.0*binwidth_x*binwidth_y);
	else if (derivative_type==2)
	// yy-derivative:
	DDPhi_FD=(potential_array[iplus*nx+j]+potential_array[iminus*nx+j]-2.0*potential_array[i*nx+j])/(binwidth_y*binwidth_y);
	else
	{
		printf("ERROR in FD_DDPhi: Wrong derivative type supplied. Supplied %d\n", derivative_type);
		exit(1);
	}

	return DDPhi_FD;
}


void readpotentialplanes_header(struct fitsheader *FITSheader, struct plane_2D *Plane, int number_of_planes, int plane_shift, int convergence_direct, int nx, int realization, double source_redshift, double *source_comoving_distance)
{
	int i;
	int plane_number;
	char filename[1000];
	// double max_survey_angle;

	int ic_number;
	int realization_number; // realization number of plane (not of WL map being generated, the latter is the variable realization given as an argument to this function).
    int mirrot, shiftx, shifty;
    

		FILE *filelist;
		//		filelist=fopen("Plane_List.txt", "r");
	
		plane_number=number_of_planes;
		printf("Number of Planes to read in externally: %d\n", number_of_planes);

		while (plane_number>0)
		{
			//if (sscanf(line, "%s", filename) == 1)
			plane_number--;

            if (plane_shift==0)
            {
                // Select proper plane (mix ic's and realizations):
                ic_number=cosmology_sampler[realization-1][2*plane_number];
                realization_number=cosmology_sampler[realization-1][2*plane_number+1];
            }
            else // do this with minimal number of planes (few mutually exclusive boxslices) and plane shifting (use specific random number file for this purpose with five integer entries per plane):
            {
                ic_number=cosmology_sampler[realization-1][5*plane_number];
                realization_number=cosmology_sampler[realization-1][5*plane_number+1];
                mirrot=cosmology_sampler[realization-1][5*plane_number+2];
                shiftx=cosmology_sampler[realization-1][5*plane_number+3];
                shifty=cosmology_sampler[realization-1][5*plane_number+4];
            }

            
//<AP>:this line needs to be removed once we run more simulations
    ic_number = 1;
//</AP>
            
            //////////////////////////
            // TEMP:
            /*
            // while (ic_number>=9) { ic_number=ic_number/9; }
            while (realization_number>=9) { realization_number=realization_number/9; }
            // ic_number++;
            realization_number++;
            ic_number=1;
             */
            //////////////////////////
            

			printf("plane_number, parameters.realization, ic_number, realization_number: %d %d %d %d\n", plane_number, realization, ic_number, realization_number);

			if (convergence_direct!=0) sprintf(filename, "%s/%s_ic%d/%s/%s_%s_ic%d_%04dxy_%04dr_%03dp.%s", parameters.plane_urpath, parameters.simulation_codename, ic_number, parameters.planes_folder, parameters.density_basename, parameters.simulation_codename, ic_number, nx, realization_number, plane_number, parameters.extension); // if doing convergence directly from density planes instead of potential planes (with derivatives).
			else sprintf(filename, "%s/%s_ic%d/%s/%s_%s_ic%d_%04dxy_%04dr_%03dp.%s", parameters.plane_urpath, parameters.simulation_codename, ic_number, parameters.planes_folder, parameters.potential_basename, parameters.simulation_codename, ic_number, nx, realization_number, plane_number, parameters.extension); // regular case.
			/*
			if (parameters.convergence_direct!=0) sprintf(filename, "%s/%s_%s_%04dxy_%04dr_%03dp.%s", parameters.plane_path, parameters.density_basename, parameters.simulation_codename, parameters.nx, parameters.realization, plane_number, parameters.extension); // if doing convergence directly from density planes instead of potential planes (with derivatives).
			else sprintf(filename, "%s/%s_%s_%04dxy_%04dr_%03dp.%s", parameters.plane_path, parameters.potential_basename, parameters.simulation_codename, parameters.nx, parameters.realization, plane_number, parameters.extension); // regular case.
			*/
			{
				printf("Reading Header of Plane File %s \n...", filename);
				timing();

#ifdef HAVE_FITS
				readFITSheader(filename, plane_number, FITSheader, Plane); // NOTE: This function already assigns some values to the Plane structure, not just to the FITSheader structure.
#else
                exit(1); // currently not implemented, must have FITS.
#endif
                
                FITSheader_to_parameters(FITSheader, Plane, plane_number, number_of_planes); // NOTE: This function transfers values from FITSheader to parameters and Plane structures.
                
                
 
                
                
				if (plane_number==number_of_planes-1) // First plane is the farthest plane, determines the max possible opening angle of the simulated survey.
                    // Re-adjust quantities from initial parameters file automatically accordingly:
				{
                
                
                    ////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////
                    // Initialize dark energy and comoving distance caluculation:
                    // (makes splines, which can then be evaluated with quick calls)
                    
					initialize_darkenergy(); // initializes w(z) dark energy profile for calculation of comoving distances which do not come from Gadget-2 snapshots (like the one for the sources).
					*source_comoving_distance=calculate_comoving_distance(1.0/(1.0+source_redshift))/1000.0; // comoving distance is originally calculated in kpc/h.
					// Include others here as well of the above overwritten quantities, as needed.
                    
                    initialize_chi_in_Mpc(); // initializes spline lookup of comoving distance with get_chi() instead of computing it from scratch each time with calculate_comoving_distance(a). This is for later use (necessary for computation with galaxy catalogue, where every galaxy has its own comoving distance to be calculated.

                    // Note: both initialization functions above use flags to prevent reinitialization if called again. initialize_darkenergy() uses parameters.darkenergy_initialized while initialize_chi_in_Mpc() sets parameters.chi_initialized. They are initially set to zero in main() and then changed to 1 upon first call of the above functions. When they are set to 1, the above functions return immediately when called, without executing anything.
                    
                    // Note 2: Note that initialize_darkenergy() needs to know all cosmological parameters, which are only known to the code after the first potential plane is read in. So this is the first place in the code when the above functions can be called.
                    
                    
                    // Small incomplete check that comoving distance spline was initialized properly (comparison to the direct computation value obtained above):
                    double comoving_check;
                    comoving_check=get_chi_in_Mpc(source_redshift);
                    printf("Comoving distance of source: parameters.source_comoving_distance (direct computation): %e, via spline: %e (at source redshift %e).\n", *source_comoving_distance, comoving_check, source_redshift);
                    if (fabs(*source_comoving_distance-comoving_check)>0.0001)
                    {
                        printf("Error: The comoving distance computed from the spline (%e) via get_chi_in_Mpc(z) is not within requested precision compared to the from-scratch calculation (5e) with calculate_comoving_distance(a) at source redshift %e, for which this has been checked. Something is likely wrong with spline initialization or evaluation. Aborting.\n");
                        MPI_Abort(MPI_COMM_WORLD, 1);
                        exit(1);
                        // Note: this error can also occur is parameters.source_comoving distance is chosen very small, ~0.001, because the spline there is not very precise with the default spline parameters chosen. In that case decrease the precision requirement above of adjust the spline parameters to cover that z-region accurately.
                    }
 

				}
				else
				{
					// Check for consistency of lens plane set (all need to have same cosmological parameters):
					if (parameters.H_0!=FITSheader[plane_number].H_0 || parameters.Omega_m!=FITSheader[plane_number].Omega_m || parameters.Omega_Lambda!=FITSheader[plane_number].Omega_Lambda || parameters.w0!=FITSheader[plane_number].w0 || parameters.wa!=FITSheader[plane_number].wa || parameters.ns!=FITSheader[plane_number].ns || parameters.sigma_8!=FITSheader[plane_number].sigma_8)
					{
						printf("ERROR: Lens planes used from inconsistent cosmologies:\nCosmology (H_0=%e, Omega_m=%e, Omega_Lambda=%e, w0=%e, wa=%e, ns=%e, sigma_8=%e),\nPlane %d: (H_0=%e, Omega_m=%e, Omega_Lambda=%e, w0=%e, wa=%e, ns=%e, sigma_8=%e).\nAborting.\n", parameters.H_0, parameters.Omega_m, parameters.Omega_Lambda, parameters.w0, parameters.wa, parameters.ns, parameters.sigma_8, plane_number, FITSheader[plane_number].H_0, FITSheader[plane_number].Omega_m, FITSheader[plane_number].Omega_Lambda, FITSheader[plane_number].w0, FITSheader[plane_number].wa, FITSheader[plane_number].ns, FITSheader[plane_number].sigma_8);
                        MPI_Abort(MPI_COMM_WORLD, 1);
						exit(1);
					}
				}
				
				
			}
		}
		if (plane_number!=0)
		{
			printf("WARNING: parameters.number_of_planes does not equal number of planes in list of plane files in file Plane_List.txt.\n");
			// exit(1);
		}
		parameters.pass=1;
}


void FITSheader_to_parameters(struct fitsheader *FITSheader, struct plane_2D *Plane, int plane_number, int number_of_planes)
{
    int i;
    double max_survey_angle;
    
    
    parameters.nx=FITSheader[plane_number].NAXIS1;
    parameters.ny=FITSheader[plane_number].NAXIS2;
    parameters.nxny=parameters.nx*parameters.ny;
    
    for (i=0;i<6;i++)
    {
        parameters.NumPartTotal[i]=FITSheader[plane_number].NumPartTotal[i];
        parameters.mass[i]=FITSheader[plane_number].mass[i];
    }
    
    
    if (feedback>3) printf("Done allocating planes arrays.\n");
	
    Plane[plane_number].boxsize=FITSheader[plane_number].BOXSIZE;
    Plane[plane_number].binwidth_x=((double) FITSheader[plane_number].BOXSIZE)/((double) FITSheader[plane_number].NAXIS1);
    Plane[plane_number].binwidth_y=((double) FITSheader[plane_number].BOXSIZE)/((double) FITSheader[plane_number].NAXIS2);
    Plane[plane_number].physical_distance=FITSheader[plane_number].CHI; // physical distance = comoving distance, since a=1 today.
    Plane[plane_number].comoving_distance=FITSheader[plane_number].CHI;
    Plane[plane_number].scale_factor=FITSheader[plane_number].A; // scale factor at plane.
    
    Plane[plane_number].H_0=FITSheader[plane_number].H_0;
    Plane[plane_number].Omega0=FITSheader[plane_number].Omega_m;
    Plane[plane_number].OmegaLambda=FITSheader[plane_number].Omega_Lambda;
    Plane[plane_number].w0=FITSheader[plane_number].w0;
    Plane[plane_number].wa=FITSheader[plane_number].wa;
    Plane[plane_number].ns=FITSheader[plane_number].ns;
    Plane[plane_number].sigma_8=FITSheader[plane_number].sigma_8;
    Plane[plane_number].initial_condition=FITSheader[plane_number].initial_condition;
    // catchment_close, catchment_far, Rot[], and Shift[] are initialized directly by functions readFITSheader().
    Plane[plane_number].set=1; // Set flag that plane parameters initialized properly.
    
    
    if (plane_number==number_of_planes-1) // First plane is the farthest plane, determines the max possible opening angle of the simulated survey.
        // Re-adjust quantities from initial parameters file automatically accordingly:
    {
        
        
        parameters.boxsize=FITSheader[plane_number].BOXSIZE;
        
        if (parameters.survey_angle==0.0) // Automatic maximization of survey angle only if set to zero initially
            // To analyze on simulation, automatic maximization to fill farthest plane is desirable, not so though to compare different cosmological models,
            // where one wants to keep the same survey angle for all of them.
        {
            parameters.survey_angle_in_rad=(parameters.nx*Plane[plane_number].binwidth_x)/Plane[plane_number].physical_distance/(1+parameters.plane_padding); // survey andgle of survey set such that spans widest possible angle in simulation up to given redshift of sources (last lens plane).
            
            //Optional:
            parameters.survey_angle=parameters.survey_angle_in_rad/(2.0*pi)*360.0;
            printf("Automatically readjusted quantities from pre-fab plane read-in:\n");
            printf("Full Survey Angle (in degrees): %e\n", parameters.survey_angle);
        }
        else
        {
            // If survey angle chosen manually, check if it's not bigger than farthest snapshot/plane of simulation allows:
            max_survey_angle=(parameters.nx*Plane[plane_number].binwidth_x)/Plane[plane_number].physical_distance/(2.0*pi)*360.0;
            if (parameters.survey_angle>max_survey_angle && plane_number==parameters.plane_before_source)
            {
                printf("ERROR: Manually selected survey angle exceeds maximally allowed survey angle (given by size of fartherst plane/snapshot box): Survey Angle > Max Angle %e %e.  Aborting.\n", parameters.survey_angle, max_survey_angle);
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD, 1);
                exit(1);
            }
        }
        
        // Further values to be written into header of WL lensing map FITS files:
        parameters.ASPP_resolution=((double) parameters.survey_angle*3600.0)/((double) parameters.NbinsX); // Arcseconds per pixel in WL map (resolution).
        parameters.PPAM_resolution=((double) parameters.NbinsX)/((double) parameters.survey_angle*60.0); // Pixels per arcminute in WL map (resolution).
        
        // Read cosmological parameters from farthest used lens plane header (overwriting values from parameter file if necessary):
        if (parameters.H_0!=FITSheader[plane_number].H_0 || parameters.Omega_m!=FITSheader[plane_number].Omega_m || parameters.Omega_Lambda!=FITSheader[plane_number].Omega_Lambda || parameters.w0!=FITSheader[plane_number].w0 || parameters.wa!=FITSheader[plane_number].wa || parameters.ns!=FITSheader[plane_number].ns || parameters.sigma_8!=FITSheader[plane_number].sigma_8) printf("Overwriting cosmology parameter values in parameter file by values in farthest lens plane: H_0=%e, Omega_m=%e, Omega_Lambda=%e, w0=%e, wa=%e, ns=%e, sigma_8=%e.\n", FITSheader[plane_number].H_0, FITSheader[plane_number].Omega_m, FITSheader[plane_number].Omega_Lambda, FITSheader[plane_number].w0, FITSheader[plane_number].wa, FITSheader[plane_number].ns, FITSheader[plane_number].sigma_8);
        parameters.H_0=FITSheader[plane_number].H_0;
        parameters.Omega_m=FITSheader[plane_number].Omega_m;
        parameters.Omega_Lambda=FITSheader[plane_number].Omega_Lambda;
        parameters.w0=FITSheader[plane_number].w0;
        parameters.wa=FITSheader[plane_number].wa;
        parameters.ns=FITSheader[plane_number].ns;
        parameters.sigma_8=FITSheader[plane_number].sigma_8;
        parameters.initial_condition=FITSheader[plane_number].initial_condition;
        
    } // end of plane_number if clause.

    return;
}





void load_potential_plane(double *potential_array, int plane_number, int nx, int ny, int realization, int plane_before_source, double source_redshift, int plane_shift, int convergence_direct, int preload_planes, int number_of_planes, int number_of_plane_realizations, int first_sim_ic, MPI_Comm sim_comm, MPI_Win *plane_storage_window)
{
	// Construct filename for readin and read in plane file:
	char filename[1000];
	int redshift_tag;

	// Select proper plane (mix ic's and realizations):
	int ic_number;
	int realization_number;
    int mirrot, shiftx, shifty;
    
    int nxny;
    nxny=nx*ny;
    
    if (plane_shift==0) // use different random number file in this case.
    {
        ic_number=cosmology_sampler[realization-1][2*plane_number];
        realization_number=cosmology_sampler[realization-1][2*plane_number+1];
    }
    else // do this with minimal number of planes (few mutually exclusive boxslices) and plane shifting (use specific random number file for this purpose with five integer entries per plane):
    {
        ic_number=cosmology_sampler[realization-1][5*plane_number];
        realization_number=cosmology_sampler[realization-1][5*plane_number+1];
        mirrot=cosmology_sampler[realization-1][5*plane_number+2];
        shiftx=cosmology_sampler[realization-1][5*plane_number+3];
        shifty=cosmology_sampler[realization-1][5*plane_number+4];
    }
    
//<AP>:this line needs to be removed once we run more simulations
    ic_number = 1;
//</AP>
    
    
    //////////////////////////
    // TEMP:
    /*
    // while (ic_number>=9) { ic_number=ic_number/9; }
    while (realization_number>=9) { realization_number=realization_number/9; }
    // ic_number++;
    realization_number++;
    ic_number=1;
     */
    //////////////////////////
    

	printf("load plane_number, parameters.realization, ic_number, realization_number: %d %d %d %d\n", plane_number, realization, ic_number, realization_number);

    
    if (preload_planes==0)
    {
	if (plane_number==plane_before_source) // if last plane before source:
	  {
	    redshift_tag=((int) floor(100.0*source_redshift));

            if (parameters.convergence_direct!=0) sprintf(filename, "%s/%s_ic%d/%s/%s_%s_ic%d_%04dxy_%04dr_%03dp_%04dz.%s", parameters.plane_urpath, parameters.simulation_codename, ic_number, parameters.planes_folder, parameters.density_basename, parameters.simulation_codename, ic_number, nx, realization_number, plane_number, redshift_tag, parameters.extension); // when calculating convergence direct from density planes (approximation without using Poisson equation to convert to gravitational potential and derivatives).
	    else sprintf(filename, "%s/%s_ic%d/%s/%s_%s_ic%d_%04dxy_%04dr_%03dp_%04dz.%s", parameters.plane_urpath, parameters.simulation_codename, ic_number, parameters.planes_folder, parameters.potential_basename, parameters.simulation_codename, ic_number, nx, realization_number, plane_number, redshift_tag, parameters.extension); // regular case for weak lensing.

	    /*
	    if (parameters.convergence_direct!=0) sprintf(filename, "%s/%s_%s_%04dxy_%04dr_%03dp_%04dz.%s", parameters.plane_path, parameters.density_basename, parameters.simulation_codename, parameters.nx, parameters.realization, plane_number, redshift_tag, parameters.extension); // when calculating convergence direct from density planes (approximation without using Poisson equation to convert to gravitational potential and derivatives).
	    else sprintf(filename, "%s/%s_%s_%04dxy_%04dr_%03dp_%04dz.%s", parameters.plane_path, parameters.potential_basename, parameters.simulation_codename, parameters.nx, parameters.realization, plane_number, redshift_tag, parameters.extension); // regular case for weak lensing.
	    */
	  }
	else // if regular plane somewhere at beginning or in the middle of stack:
	  {
            if (convergence_direct!=0) sprintf(filename, "%s/%s_ic%d/%s/%s_%s_ic%d_%04dxy_%04dr_%03dp.%s", parameters.plane_urpath, parameters.simulation_codename, ic_number, parameters.planes_folder, parameters.density_basename, parameters.simulation_codename, ic_number, nx, realization_number, plane_number, parameters.extension);
            else sprintf(filename, "%s/%s_ic%d/%s/%s_%s_ic%d_%04dxy_%04dr_%03dp.%s", parameters.plane_urpath, parameters.simulation_codename, ic_number, parameters.planes_folder, parameters.potential_basename, parameters.simulation_codename, ic_number, nx, realization_number, plane_number, parameters.extension);

	    /*
	    if (parameters.convergence_direct!=0) sprintf(filename, "%s/%s_%s_%04dxy_%04dr_%03dp.%s", parameters.plane_path, parameters.density_basename, parameters.simulation_codename, parameters.nx, parameters.realization, plane_number, parameters.extension);
	    else sprintf(filename, "%s/%s_%s_%04dxy_%04dr_%03dp.%s", parameters.plane_path, parameters.potential_basename, parameters.simulation_codename, parameters.nx, parameters.realization, plane_number, parameters.extension);
	    */
	  }

#ifdef HAVE_FITS
    readFITSpotential_singleplane(filename, potential_array);
#else
        exit(1); // currently not implemented, must have FITS.
#endif
    }
    else // preload all planes at the beginning and then do one-sided RMA MPI_Get:
    {
        // Find on which MPI process and with what target displacement the plane has been loaded into plane_storage at the beginning of the run:
        int target_rank;
	MPI_Aint target_displacement;
        
        plane_parameters_to_displacement(plane_number, realization_number, ic_number, number_of_planes, number_of_plane_realizations, first_sim_ic, nx, ny, sim_comm, &target_rank, &target_displacement); // note that plane_number starts counting at zero, while realization_number and ic_number start counting at 1 (but plane_realization and plane_sim_ic in the function plane_parameters_to_displacement() start counting at zero).
        
        // printf("Will do RMA MPI Operation with target rank %d and target displacement %d (instance %d), loadingsize %d ...\n", target_rank, target_displacement, target_displacement/nxny, nxny);
        // fflush(stdout);
        
        // Now do one-sided RMA MPI_Get:
        MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0, *plane_storage_window);
        MPI_Get(potential_array, nxny, MPI_DOUBLE, target_rank, target_displacement, nxny, MPI_DOUBLE, *plane_storage_window);
        MPI_Win_unlock(target_rank, *plane_storage_window);
        
        // printf("Done doing RMA MPI Operation.\n");
        // fflush(stdout);

    }
    
    // Shift and rotate potential plane if use only few exclusive-sliced lens planes to save disk space.
    shift_potential_singleplane(potential_array, nx, ny, mirrot, shiftx, shifty);
    
}


///////////////////////////////////////////////////////////////////////////////////
// Functions for rotating, mirroring, and shifting lens planes during ray-tracing:
// (Not necessary if have created planes from fully rotated boxes individually for each map, but that consumes a lot of disk space.)

void shift_potential_singleplane(double *potential_array, int nx, int ny, int mirrot, int shiftx, int shifty)
{
    int i, j, ii, jj;
    double *potential_shifter, value;
    
    potential_shifter=Vector0(nx*ny);
    
    /*
    double r1, m1, m2, x, y;
     
    r1=ran2(&parameters.seed); // rotate axis if > 0.5.
    m1=ran2(&parameters.seed); // mirror axis if >0.5.
    m2=ran2(&parameters.seed);
    x=ran2(&parameters.seed); // shift axis 1 by x1*nx;
    y=ran2(&parameters.seed);
    */
     
    for (i=0; i<ny; i++)
    {
        for (j=0; j<nx; j++)
        {
            value=potential_array[i*nx+j];
            // shiftrot2D(i,j, nx, ny, r1, m1, m2, x, y, &ii, &jj);
            shiftrot2D_from_file(i, j, nx, ny, mirrot, shiftx, shifty, &ii, &jj);
            potential_shifter[ii*nx+jj]=value;
        }

    }
    
    for (i=0; i<ny; i++)
    {
        for (j=0; j<nx; j++)
        {
            potential_array[i*nx+j]=potential_shifter[i*nx+j];
        }
    }
    
    free_Vector0(potential_shifter);
    
}


void shiftrot2D_from_file(int i, int j, int nx, int ny, int mirrot, int shiftx, int shifty, int *ii, int *jj)
{
    int tempint;
    double x, y;
    
    if (mirrot>=4) // swap axes:
    {
        tempint=i;
        i=j;
        j=tempint;
        mirrot-=4;
    }
    if (mirrot>=2) // mirror x-axis:
    {
        j=nx-j;
        mirrot-=2;
    }
    if (mirrot>=1) // mirror y-axis:
    {
        i=ny-i;
    }
    
    x=shiftx/((double) MAX_RANDOM_NUMBER);
    y=shifty/((double) MAX_RANDOM_NUMBER);
    
    j=j+((int) floor(x*nx));
    if (j>=nx) j=j-nx;
    i=i+((int) floor(y*ny));
    if (i>=ny) i=i-ny;
    
    *ii=i;
    *jj=j;
    
}


#undef MAX_RANDOM_NUMBER
