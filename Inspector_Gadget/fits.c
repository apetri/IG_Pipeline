/*
 *  fits.c
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 7/11/07.
 *  Copyright 2007. All rights reserved.
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
#include "fitsio.h"

#include "main.h"
#include "mathematics.h"
#include "allocation.h"
#include "io_routines.h"
#include "2D-plane_multi.h"
#include "fits.h"
#include "weak_lensing_multi.h"

#define MAXNAME 200

void writePlaneFITSimage_f(char filename[], long naxis, long naxes[], double *imagearray, int plane_number) // ATTENTION: externally fed in image imagearray must be a 1D (pseudo-2D) array.

    /******************************************************/
    /* Create a FITS primary array containing a 2-D image */
    /******************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii;
    long  fpixel, nelements; //, exposure;
	float Z_redshift, A_scale_factor, CHI_comoving_distance, CLOSE_catchment, FAR_catchment, BOXSIZE_boxsize, GRID_resolution, GRID_space_x, GRID_space_y, H_0_parameter, Omega_m_parameter, Omega_L_parameter, w0_parameter, wa_parameter, ns_parameter, sigma8_parameter, initial_condition, N_CDM, M_CDM, N_plane, N_avg, Rot_X, Rot_Y, Rot_Z, Shift_X, Shift_Y, Shift_Z;
    char MODEL_folder[MAXNAME], MODEL[MAXNAME], USER_comment[MAXNAME];

	float *writearray;

    /* initialize FITS image parameters */
    //char filename[] = "atestfil.fit";             /* name for new FITS file */
    int bitpix   =  FLOAT_IMG; /* 32-bit(?) signed float pixel values       */
    //long naxis    =   2;  /* 2-dimensional image                            */    
    //long naxes[2] = { 300, 200 };   /* image is 300 pixels wide by 200 rows */

    /* allocate memory for the whole image */ 
	writearray=(float *)malloc(naxes[0]*naxes[1]*sizeof(float));
	assert(writearray != NULL );
	

    remove(filename);               /* Delete old file if it already exists */

    status = 0;         /* initialize status before calling fitsio routines */

    if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
         printerror( status );           /* call printerror if error occurs */

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
         printerror( status );          


    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    for (ii=0; ii<nelements; ii++) writearray[ii]=((float) imagearray[ii]);

    /* write the array of unsigned integers to the FITS file */
    //    if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, array[0], &status) )
      if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, writearray, &status) )
        printerror( status );
      
       free(writearray);  /* free previously allocated memory */


    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
   
	// Write header with parameters of potential plane into FITS file:

	Z_redshift= ((float) ((1.0/Plane[plane_number].scale_factor)-1.0));
    if ( fits_update_key(fptr, TFLOAT, "Z", &Z_redshift,
         "Redshift of plane", &status) )
         printerror( status );     
	
	A_scale_factor=((float) Plane[plane_number].scale_factor);
    if ( fits_update_key(fptr, TFLOAT, "A", &A_scale_factor,
         "Scale factor at plane (norm. to a=1 today)", &status) )
         printerror( status );     
		 
	CHI_comoving_distance =((float) Plane[plane_number].comoving_distance/1000.0); // comoving distance is so far in Mpc/h
    if ( fits_update_key(fptr, TFLOAT, "CHI", &CHI_comoving_distance,
         "Comoving distance of plane to obs (in Mpc/h)", &status) )
         printerror( status );
		 
	CLOSE_catchment=((float) Plane[plane_number].catchment_close/1000.0);	// division by 1000 to convert from kpc/h (used by Gadget-2 as program units by default) to Mpc/h. 	 
    if ( fits_update_key(fptr, TFLOAT, "CLOSE", &CLOSE_catchment,
         "Catchment area close (in Mpc/h from plane)", &status) )
         printerror( status );
		 
	FAR_catchment=((float) Plane[plane_number].catchment_far/1000.0);	// division by 1000 to convert from kpc/h (used by Gadget-2 as program units by default) to Mpc/h. 	 
    if ( fits_update_key(fptr, TFLOAT, "FAR", &FAR_catchment,
         "Catchment area far (in Mpc/h from plane)", &status) )
         printerror( status );  
	
	BOXSIZE_boxsize=((float) Plane[plane_number].boxsize/1000.0);	// division by 1000 to convert from kpc/h (used by Gadget-2 as program units by default) to Mpc/h. 	 
    if ( fits_update_key(fptr, TFLOAT, "BOXSIZE", &BOXSIZE_boxsize,
         "Lateral boxsize of snapshot / plane (in Mpc/h)", &status) )
         printerror( status );
		 
	GRID_resolution=((float) parameters.nx);	 	 
    if ( fits_update_key(fptr, TFLOAT, "GRID", &GRID_resolution,
		"Resolution of grid on plane", &status) )
		printerror( status );

	GRID_space_x=((float) Plane[plane_number].binwidth_x/1000.0);	// division by 1000 to convert from kpc/h (used by Gadget-2 as program units by default) to Mpc/h. 	 
    if ( fits_update_key(fptr, TFLOAT, "RES_X", &GRID_space_x,
         "Grid point spacing x-dir (in Mpc/h)", &status) )
         printerror( status );
	
	GRID_space_y=((float) Plane[plane_number].binwidth_y/1000.0);	// division by 1000 to convert from kpc/h (used by Gadget-2 as program units by default) to Mpc/h. 	 
    if ( fits_update_key(fptr, TFLOAT, "RES_Y", &GRID_space_y,
         "Grid point spacing y-dir (in Mpc/h)", &status) )
         printerror( status );  

	H_0_parameter =((float) Plane[plane_number].H_0);
    if ( fits_update_key(fptr, TFLOAT, "H_0", &H_0_parameter,
         "Hubble constant today (in km/s/Mpc)", &status) )
         printerror( status );  

	Omega_m_parameter =((float) Plane[plane_number].Omega0);
    if ( fits_update_key(fptr, TFLOAT, "OMEGA_M", &Omega_m_parameter,
         "Omega_matter (matter density)", &status) )
         printerror( status );
	
	Omega_L_parameter =((float) Plane[plane_number].OmegaLambda);
    if ( fits_update_key(fptr, TFLOAT, "OMEGA_L", &Omega_L_parameter,
         "Omega_Lambda (dark energy density)", &status) )
         printerror( status );
    
    w0_parameter =((float) Plane[plane_number].w0);
    if ( fits_update_key(fptr, TFLOAT, "W_0", &w0_parameter,
         "w_0 (dark energy EOS)", &status) )
        printerror( status );
    
    wa_parameter =((float) Plane[plane_number].wa);
    if ( fits_update_key(fptr, TFLOAT, "W_A", &wa_parameter,
         "w_a (dark energy EOS)", &status) )
        printerror( status );
    
    ns_parameter =((float) Plane[plane_number].ns);
    if ( fits_update_key(fptr, TFLOAT, "N_S", &ns_parameter,
         "n_s (scalar spectral index)", &status) )
        printerror( status );
    
    sigma8_parameter =((float) Plane[plane_number].sigma_8);
    if ( fits_update_key(fptr, TFLOAT, "SIGMA_8", &sigma8_parameter,
         "sigma_8 (amplitude of fluct.)", &status) )
        printerror( status );
    
    initial_condition =((float) Plane[plane_number].initial_condition);
    if ( fits_update_key(fptr, TFLOAT, "IC", &initial_condition,
         "Initial Condition", &status) )
        printerror( status );

    
	N_CDM =((float) Plane[plane_number].NumPartTotal[1]);
    if ( fits_update_key(fptr, TFLOAT, "N_CDM", &N_CDM,
         "CDM: number particles in N-body sim", &status) )
         printerror( status );
		 
	M_CDM =((float) Plane[plane_number].mass[1]);
    if ( fits_update_key(fptr, TFLOAT, "M_CDM", &M_CDM,
         "CDM: mass of particle (in 10^10 solarmass/h)", &status) )
         printerror( status );

	N_plane=((float) Plane[plane_number].particles_written[parameters.seed_block]);
    if ( fits_update_key(fptr, TFLOAT, "N_PLANE", &N_plane,
         "Number of particles written in this plane", &status) )
         printerror( status );
		 
	N_avg=((float) Plane[plane_number].particles_written[parameters.seed_block]/((float) parameters.nxny));
    if ( fits_update_key(fptr, TFLOAT, "AVG_PART", &N_avg,
         "Average number of particles per plane grid cell", &status) )
         printerror( status );

    // Snapshot Box Rotations and Shifts:

    Rot_X=((float) Plane[plane_number].rRot[0][parameters.seed_block]);
    if ( fits_update_key(fptr, TFLOAT, "ROT_X", &Rot_X,
			 "New x axis (horizontal right) is old GSS axis...", &status) )
      printerror( status );

    Rot_Y=((float) Plane[plane_number].rRot[1][parameters.seed_block]);
    if ( fits_update_key(fptr, TFLOAT, "ROT_Y", &Rot_Y,
                         "New y axis (vertical up) is old GSS axis...", &status) )
      printerror( status );

    Rot_Z=((float) Plane[plane_number].rRot[2][parameters.seed_block]);
    if ( fits_update_key(fptr, TFLOAT, "ROT_Z", &Rot_Z,
                         "New z axis (away from obs) is old GSS axis...", &status) )
      printerror( status );

    Shift_X=((float) Plane[plane_number].rShift[0][parameters.seed_block]);  // in Mpc/h, has already been converted to Mpc/h earlier upon initialization.
    if ( fits_update_key(fptr, TFLOAT, "SHIFT_X", &Shift_X,
                         "Box center shift along new x axis (in Mpc/h)", &status) )
      printerror( status );

    Shift_Y=((float) Plane[plane_number].rShift[1][parameters.seed_block]);
    if ( fits_update_key(fptr, TFLOAT, "SHIFT_Y", &Shift_Y,
                         "Box center shift along new y axis (in Mpc/h)", &status) )
      printerror( status );

    Shift_Z=((float) Plane[plane_number].rShift[2][parameters.seed_block]);
    if ( fits_update_key(fptr, TFLOAT, "SHIFT_Z", &Shift_Z,
                         "Box center shift along new z axis (in Mpc/h)", &status) )
      printerror( status );



    strncpy(MODEL_folder, parameters.simulation_codename, MAXNAME);
    if ( fits_update_key(fptr, TSTRING, "M-CODE", &MODEL_folder,
			 "", &status) )
      printerror( status );

	strncpy(MODEL, parameters.modelname, MAXNAME);
	if ( fits_update_key(fptr, TSTRING, "MODEL", &MODEL,
         "", &status) )
         printerror( status );
			 	 
	strncpy(USER_comment, parameters.plane_comment, MAXNAME);
	if ( fits_update_key(fptr, TSTRING, "COMM", &USER_comment,
         "", &status) )
         printerror( status );


	// Done writing header information into FITS file.


    if ( fits_close_file(fptr, &status) )                /* close the file */
         printerror( status );           

    return;
}


void writeWLmapFITSimage_f(char filename[], long naxis, long naxes[], double *imagearray) // ATTENTION: externally fed in image imagearray must be a true 2D array.

    /******************************************************/
    /* Create a FITS primary array containing a 2-D image */
    /******************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii;
    long  fpixel, nelements; //, exposure;
	float Z_redshift, A_scale_factor, CHI_comoving_distance, BOXSIZE_boxsize, H_0_parameter, Omega_m_parameter, Omega_L_parameter, w0_parameter, wa_parameter, ns_parameter, sigma8_parameter, initial_condition, ANGLE_survey_angle, GRID_resolution, GRID_space_x, GRID_space_y, MAP_resolution, ASPP_resolution, PPAM_resolution, N_CDM, M_CDM, PLANES, CHI_plane;
	int plane_number;
	char MODEL_folder[MAXNAME], MODEL[MAXNAME], label[MAXNAME], comment[MAXNAME], USER_comment[MAXNAME];

	float *writearray;

    /* initialize FITS image parameters */
    //char filename[] = "atestfil.fit";             /* name for new FITS file */
    int bitpix   =  FLOAT_IMG; /* 32-bit(?) signed float pixel values       */
    //long naxis    =   2;  /* 2-dimensional image                            */    
    //long naxes[2] = { 300, 200 };   /* image is 300 pixels wide by 200 rows */

    /* allocate memory for the whole image */ 
	writearray=(float *)malloc(naxes[0]*naxes[1]*sizeof(float));
	assert(writearray != NULL );
	

    remove(filename);               /* Delete old file if it already exists */

    status = 0;         /* initialize status before calling fitsio routines */

    if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
         printerror( status );           /* call printerror if error occurs */

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
         printerror( status );          


    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    for (ii=0; ii<nelements; ii++) writearray[ii]=((float) imagearray[ii]);

    /* write the array of unsigned integers to the FITS file */
    if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, writearray, &status) )
        printerror( status );
      
 
    free(writearray);  /* free previously allocated memory */

    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
   

	// Write header with parameters of potential plane into FITS file:

	Z_redshift= ((float) parameters.source_redshift);	
    if ( fits_update_key(fptr, TFLOAT, "Z", &Z_redshift,
         "Redshift of sources (plane or mean)", &status) )
         printerror( status );
	
	A_scale_factor=((float) parameters.source_scale_factor);
    if ( fits_update_key(fptr, TFLOAT, "A", &A_scale_factor,
         "Scale factor at sources (norm. to a=1 today)", &status) )
         printerror( status );     
		 	 
	CHI_comoving_distance = ((float) parameters.source_comoving_distance); // comoving distance of sources (or source plane) to observer.
    if ( fits_update_key(fptr, TFLOAT, "CHI", &CHI_comoving_distance,
         "Comoving distance of sources (in Mpc/h)", &status) )
         printerror( status );
			 	 
	ANGLE_survey_angle = ((float) parameters.survey_angle); // survey angle in degrees.
    if ( fits_update_key(fptr, TFLOAT, "ANGLE", &ANGLE_survey_angle,
         "Survey angle (in degrees)", &status) )
         printerror( status );
	
	BOXSIZE_boxsize=((float) parameters.boxsize); // NOT: division by 1000 converts from kpc/h to Mpc/h.	 	 
    if ( fits_update_key(fptr, TFLOAT, "BOXSIZE", &BOXSIZE_boxsize,
         "Boxsize of snapshot cube and plane (in Mpc/h)", &status) )
         printerror( status );
		 
	H_0_parameter =((float) parameters.H_0);
    if ( fits_update_key(fptr, TFLOAT, "H_0", &H_0_parameter,
         "Hubble constant today (in km/s/Mpc)", &status) )
         printerror( status );  
	
	Omega_m_parameter =((float) parameters.Omega_m);
    if ( fits_update_key(fptr, TFLOAT, "OMEGA_M", &Omega_m_parameter,
         "Omega_matter (relative matter density)", &status) )
         printerror( status );
	
	Omega_L_parameter =((float) parameters.Omega_Lambda);
    if ( fits_update_key(fptr, TFLOAT, "OMEGA_L", &Omega_L_parameter,
         "Omega_Lambda (relative dark energy density)", &status) )
         printerror( status );

    w0_parameter =((float) Plane[plane_number].w0);
    if ( fits_update_key(fptr, TFLOAT, "W_0", &w0_parameter,
                         "w_0 (dark energy EOS)", &status) )
        printerror( status );
    
    wa_parameter =((float) Plane[plane_number].wa);
    if ( fits_update_key(fptr, TFLOAT, "W_A", &wa_parameter,
                         "w_a (dark energy EOS)", &status) )
        printerror( status );
    
    ns_parameter =((float) Plane[plane_number].ns);
    if ( fits_update_key(fptr, TFLOAT, "N_S", &ns_parameter,
                         "n_s (scalar spectral index)", &status) )
        printerror( status );
    
    sigma8_parameter =((float) Plane[plane_number].sigma_8);
    if ( fits_update_key(fptr, TFLOAT, "SIGMA_8", &sigma8_parameter,
                         "sigma_8 (amplitude of fluct.)", &status) )
        printerror( status );
    
    initial_condition =((float) Plane[plane_number].initial_condition);
    if ( fits_update_key(fptr, TFLOAT, "IC", &initial_condition,
                         "Initial Condition", &status) )
        printerror( status );

    
	GRID_resolution=((float) parameters.nx);	 	 
    if ( fits_update_key(fptr, TFLOAT, "GRID", &GRID_resolution,
         "Resolution of grid on planes", &status) )
         printerror( status );

	GRID_space_x=((float) parameters.boxsize/(float) parameters.nx); // parameters.boxsize is already in Mpc/h during the "second pass", i.e. during creating WL maps from potential planes.
    if ( fits_update_key(fptr, TFLOAT, "RES_X", &GRID_space_x,
         "Grid point spacing x-dir (in Mpc/h)", &status) )
         printerror( status );
	
	GRID_space_y=((float) parameters.boxsize/(float) parameters.ny); // parameters.boxsize is already in Mpc/h during the "second pass", i.e. during creating WL maps from potential planes.
    if ( fits_update_key(fptr, TFLOAT, "RES_Y", &GRID_space_y,
         "Grid point spacing y-dir (in Mpc/h)", &status) )
         printerror( status );  
		 
	MAP_resolution=((float) parameters.NbinsX);	 	 
    if ( fits_update_key(fptr, TFLOAT, "MAP", &MAP_resolution,
         "Pixels in WL Map", &status) )
         printerror( status );
		 
	ASPP_resolution=((float) parameters.ASPP_resolution);	 	 
    if ( fits_update_key(fptr, TFLOAT, "ASPP", &ASPP_resolution,
         "Map Resolution: arcsec per pixel", &status) )
         printerror( status );     

	PPAM_resolution=((float) parameters.PPAM_resolution);	 	 
    if ( fits_update_key(fptr, TFLOAT, "PPAM", &PPAM_resolution,
         "Map Resolution: pixels per arcminute", &status) )
         printerror( status );
		 
	N_CDM =((float) parameters.NumPartTotal[1]);
    if ( fits_update_key(fptr, TFLOAT, "N_CDM", &N_CDM,
         "CDM: number particles in N-body sim", &status) )
         printerror( status );
		 
	M_CDM =((float) parameters.mass[1]);
    if ( fits_update_key(fptr, TFLOAT, "M_CDM", &M_CDM,
         "CDM: mass of particle (in 10^10 solarmass/h)", &status) )
         printerror( status );


    strncpy(MODEL_folder, parameters.simulation_codename, MAXNAME);
    if ( fits_update_key(fptr, TSTRING, "M-CODE", &MODEL_folder,
                         "", &status) )
      printerror( status );

	strncpy(MODEL, parameters.modelname, MAXNAME);
	if ( fits_update_key(fptr, TSTRING, "MODEL", &MODEL,
         "", &status) )
         printerror( status );

	PLANES =((float) parameters.plane_before_source+1);
    if ( fits_update_key(fptr, TFLOAT, "PLANES", &PLANES,
         "Number of potential planes used", &status) )
         printerror( status );

	for (plane_number=parameters.first_plane;plane_number<=parameters.last_plane; plane_number++)
	{
		sprintf(label, "CHI_P%03d", plane_number);
		sprintf(comment, "Comoving distance of Plane %03d (in Mpc/h)", plane_number);
		CHI_plane =((float) Plane[plane_number].comoving_distance); // comoving distance is so far in Mpc/h
		if ( fits_update_key(fptr, TFLOAT, label, &CHI_plane,
			comment, &status) )
			printerror( status );
	}

	strncpy(USER_comment, parameters.WL_map_comment, MAXNAME);
	if ( fits_update_key(fptr, TSTRING, "COMM", &USER_comment,
         "", &status) )
         printerror( status );


	// Done writing header information into FITS file.


    if ( fits_close_file(fptr, &status) )                /* close the file */
         printerror( status );           

    return;
}



void writeRayPlaneFITSimage_f(char filename[], long naxis, long naxes[], double *imagearray, int plane_number) // ATTENTION: externally fed in image imagearray must be a 1D-(pseudo-2D) array.

     /******************************************************/
     /* Create a FITS primary array containing a 2-D image */
     /******************************************************/
{
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status, ii;
  long  fpixel, nelements; //, exposure;
  float Z_redshift, A_scale_factor, CHI_comoving_distance, CLOSE_catchment, FAR_catchment, ANGLE_survey_angle, BOXSIZE_boxsize, GRID_resolution, GRID_space_x, GRID_space_y, H_0_parameter, Omega_m_parameter, Omega_L_parameter, w0_parameter, wa_parameter, ns_parameter, sigma8_parameter, initial_condition, N_CDM, M_CDM, N_plane, N_avg, Rot_X, Rot_Y, Rot_Z, Shift_X, Shift_Y, Shift_Z;
  char MODEL_folder[MAXNAME], MODEL[MAXNAME], USER_comment[MAXNAME];

  float *writearray;

  nelements=1;
  for (ii=0; ii<naxis; ii++) nelements*=naxes[ii];  /* number of pixels to write */
  

  /* initialize FITS image parameters */
  //char filename[] = "atestfil.fit";             /* name for new FITS file */
  int bitpix   =  FLOAT_IMG; /* 32-bit(?) signed float pixel values       */
  //long naxis    =   2;  /* 2-dimensional image                            */    
  //long naxes[2] = { 300, 200 };   /* image is 300 pixels wide by 200 rows */

  /* allocate memory for the whole image */ 
 writearray=(float *) malloc(nelements*sizeof(float));
    assert(writearray!=NULL);


  remove(filename);               /* Delete old file if it already exists */

  status = 0;         /* initialize status before calling fitsio routines */

  if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
    printerror( status );           /* call printerror if error occurs */

  /* write the required keywords for the primary array image.     */
  /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
  /* a FITS image with BITPIX = 16 (signed short integers) with   */
  /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
  /* FITS uses to store unsigned integers.  Note that the BSCALE  */
  /* and BZERO keywords will be automatically written by cfitsio  */
  /* in this case.                                                */

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    printerror( status );          

  fpixel = 1;                               /* first pixel to write      */

  for (ii=0; ii<nelements; ii++) writearray[ii]=((float) imagearray[ii]);

  /* write the array of unsigned integers to the FITS file */
  // if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, array[0], &status) )
  if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, writearray, &status) )
    printerror( status );
      
  free(writearray);  /* free previously allocated memory */


  /* write another optional keyword to the header */
  /* Note that the ADDRESS of the value is passed in the routine */
   

  // Write header with parameters of potential plane into FITS file:

  Z_redshift= ((float) ((1.0/Plane[plane_number].scale_factor)-1.0));
  if ( fits_update_key(fptr, TFLOAT, "Z", &Z_redshift,
		       "Redshift of plane", &status) )
    printerror( status );     
  
  A_scale_factor=((float) Plane[plane_number].scale_factor);
  if ( fits_update_key(fptr, TFLOAT, "A", &A_scale_factor,
		       "Scale factor at plane (norm. to a=1 today)", &status) )
    printerror( status );     
   
  CHI_comoving_distance =((float) Plane[plane_number].comoving_distance); // comoving distance is so far in Mpc/h
  if ( fits_update_key(fptr, TFLOAT, "CHI", &CHI_comoving_distance,
		       "Comoving distance of plane to obs (in Mpc/h)", &status) )
    printerror( status );
   
  CLOSE_catchment=((float) Plane[plane_number].catchment_close); // in Mpc/h.
  if ( fits_update_key(fptr, TFLOAT, "CLOSE", &CLOSE_catchment,
		       "Catchment area close (in Mpc/h from plane)", &status) )
    printerror( status );
   
  FAR_catchment=((float) Plane[plane_number].catchment_far); // in Mpc/h.
  if ( fits_update_key(fptr, TFLOAT, "FAR", &FAR_catchment,
		       "Catchment area far (in Mpc/h from plane)", &status) )
    printerror( status );  
  
  ANGLE_survey_angle = ((float) parameters.survey_angle); // survey angle in degrees.
  if ( fits_update_key(fptr, TFLOAT, "ANGLE", &ANGLE_survey_angle,
		       "Survey angle (in degrees)", &status) )
    printerror( status );

  BOXSIZE_boxsize=((float) Plane[plane_number].boxsize); // in Mpc/h.
  if ( fits_update_key(fptr, TFLOAT, "BOXSIZE", &BOXSIZE_boxsize,
		       "Lateral boxsize of snapshot / plane (in Mpc/h)", &status) )
    printerror( status );
   
  GRID_resolution=((float) parameters.nx);  
  if ( fits_update_key(fptr, TFLOAT, "GRID", &GRID_resolution,
		       "Resolution of grid on plane", &status) )
    printerror( status );

  GRID_space_x=((float) Plane[plane_number].binwidth_x); // in Mpc/h.
  if ( fits_update_key(fptr, TFLOAT, "RES_X", &GRID_space_x,
		       "Grid point spacing x-dir (in Mpc/h)", &status) )
    printerror( status );
  
  GRID_space_y=((float) Plane[plane_number].binwidth_y); // in Mpc/h.
  if ( fits_update_key(fptr, TFLOAT, "RES_Y", &GRID_space_y,
		       "Grid point spacing y-dir (in Mpc/h)", &status) )
    printerror( status );  

  H_0_parameter =((float) parameters.H_0);
  if ( fits_update_key(fptr, TFLOAT, "H_0", &H_0_parameter,
		       "Hubble constant today (in km/s/Mpc)", &status) )
    printerror( status );  

  Omega_m_parameter =((float) parameters.Omega_m);
  if ( fits_update_key(fptr, TFLOAT, "OMEGA_M", &Omega_m_parameter,
		       "Omega_matter (matter density)", &status) )
    printerror( status );
  
  Omega_L_parameter =((float) parameters.Omega_Lambda);
  if ( fits_update_key(fptr, TFLOAT, "OMEGA_L", &Omega_L_parameter,
		       "Omega_Lambda (dark energy density)", &status) )
    printerror( status );

    w0_parameter =((float) Plane[plane_number].w0);
    if ( fits_update_key(fptr, TFLOAT, "W_0", &w0_parameter,
                         "w_0 (dark energy EOS)", &status) )
        printerror( status );
    
    wa_parameter =((float) Plane[plane_number].wa);
    if ( fits_update_key(fptr, TFLOAT, "W_A", &wa_parameter,
                         "w_a (dark energy EOS)", &status) )
        printerror( status );
    
    ns_parameter =((float) Plane[plane_number].ns);
    if ( fits_update_key(fptr, TFLOAT, "N_S", &ns_parameter,
                         "n_s (scalar spectral index)", &status) )
        printerror( status );
    
    sigma8_parameter =((float) Plane[plane_number].sigma_8);
    if ( fits_update_key(fptr, TFLOAT, "SIGMA_8", &sigma8_parameter,
                         "sigma_8 (amplitude of fluct.)", &status) )
        printerror( status );
    
    initial_condition =((float) Plane[plane_number].initial_condition);
    if ( fits_update_key(fptr, TFLOAT, "IC", &initial_condition,
                         "Initial Condition", &status) )
        printerror( status );
 
    
  N_CDM =((float) Plane[plane_number].NumPartTotal[1]);
  if ( fits_update_key(fptr, TFLOAT, "N_CDM", &N_CDM,
		       "CDM: number particles in N-body sim", &status) )
    printerror( status );
   
  M_CDM =((float) Plane[plane_number].mass[1]);
  if ( fits_update_key(fptr, TFLOAT, "M_CDM", &M_CDM,
		       "CDM: mass of particle (in 10^10 solarmass/h)", &status) )
    printerror( status );

  N_plane=((float) Plane[plane_number].particles_written[parameters.seed_block]);
  if ( fits_update_key(fptr, TFLOAT, "N_PLANE", &N_plane,
		       "Number of particles written in this plane", &status) )
    printerror( status );
   
  N_avg=((float) Plane[plane_number].particles_written[parameters.seed_block]/((float) parameters.nxny));
  if ( fits_update_key(fptr, TFLOAT, "AVG_PART", &N_avg,
		       "Average number of particles per plane grid cell", &status) )
    printerror( status );

  // Snapshot Box Rotations and Shifts:

  Rot_X=((float) Plane[plane_number].Rot[0]);
  if ( fits_update_key(fptr, TFLOAT, "ROT_X", &Rot_X,
		       "New x axis (horizontal right) is old GSS axis...", &status) )
    printerror( status );

  Rot_Y=((float) Plane[plane_number].Rot[1]);
  if ( fits_update_key(fptr, TFLOAT, "ROT_Y", &Rot_Y,
		       "New y axis (vertical up) is old GSS axis...", &status) )
    printerror( status );

  Rot_Z=((float) Plane[plane_number].Rot[2]);
  if ( fits_update_key(fptr, TFLOAT, "ROT_Z", &Rot_Z,
		       "New z axis (away from obs) is old GSS axis...", &status) )
    printerror( status );

  Shift_X=((float) Plane[plane_number].Shift[0]);  // in Mpc/h
  if ( fits_update_key(fptr, TFLOAT, "SHIFT_X", &Shift_X,
		       "Box center shift along new x axis (in Mpc/h)", &status) )
    printerror( status );

  Shift_Y=((float) Plane[plane_number].Shift[1]);
  if ( fits_update_key(fptr, TFLOAT, "SHIFT_Y", &Shift_Y,
		       "Box center shift along new y axis (in Mpc/h)", &status) )
    printerror( status );

  Shift_Z=((float) Plane[plane_number].Shift[2]);
  if ( fits_update_key(fptr, TFLOAT, "SHIFT_Z", &Shift_Z,
		       "Box center shift along new z axis (in Mpc/h)", &status) )
    printerror( status );

  strncpy(MODEL_folder, parameters.simulation_codename, MAXNAME);
  if ( fits_update_key(fptr, TSTRING, "M-CODE", &MODEL_folder,
		       "", &status) )
    printerror( status );

  strncpy(MODEL, parameters.modelname, MAXNAME);
  if ( fits_update_key(fptr, TSTRING, "MODEL", &MODEL,
		       "", &status) )
    printerror( status );
   
  strncpy(USER_comment, parameters.plane_comment, MAXNAME);
  if ( fits_update_key(fptr, TSTRING, "COMM", &USER_comment,
		       "", &status) )
    printerror( status );


  // Done writing header information into FITS file.


  if ( fits_close_file(fptr, &status) )                /* close the file */
    printerror( status );           

  return;
}





void readFITSheader (char filename[], int plane_number, struct fitsheader *FITSheader, struct plane_2D *Plane)

    /**********************************************************************/
    /* Print out all the header keywords in all extensions of a FITS file */
    /**********************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */

    int status, nkeys, keypos, hdutype, ii, jj;
    //char filename[]  = "atestfil.fit";     /* name of existing FITS file   */
    char card[FLEN_CARD];   /* standard string lengths defined in fitsioc.h */

    status = 0;

    if ( fits_open_file(&fptr, filename, READONLY, &status) ) 
         printerror( status );

    // NOTE: This really should be done better than having code write out and actual file containing a keyword card, which is then read back in again right away.
	FILE *tempfile;
	char tempname[200];
#if defined(MPI_COMPILE)
	sprintf(tempname, "tempfile_%d.txt", parameters.process_number); // necessary to avoid processors writing into same directory to overwrite each others tempfile, which would lead to corrupted header read in.
#else
	sprintf(tempname, "tempfile.txt"); // with condor, each process is copied to its own local folder, so there is no danger of interference.
#endif
	tempfile=fopen(tempname, "w");

    /* attempt to move to next HDU, until we get an EOF error */
    for (ii = 1; !(fits_movabs_hdu(fptr, ii, &hdutype, &status) ); ii++) 
    {
        /* get no. of keywords */
        if (fits_get_hdrpos(fptr, &nkeys, &keypos, &status) )
            printerror( status );

        printf("Header listing for HDU #%d:\n", ii);
        for (jj = 1; jj <= nkeys; jj++)  {
            if ( fits_read_record(fptr, jj, card, &status) )
                 printerror( status );

            printf("%s\n", card); /* print the keyword card */
			fprintf(tempfile, "%s\n", card);
        }
        printf("END\n\n");  /* terminate listing with END */
    }

    if (status == END_OF_FILE)   /* status values are defined in fitsioc.h */
        status = 0;              /* got the expected EOF error; reset = 0  */
    else
       printerror( status );     /* got an unexpected error                */

    if ( fits_close_file(fptr, &status) )
         printerror( status );
		 
	fclose(tempfile);

	/////////////////////////////////
	// Read out header parameters:

	int NAXIS, NAXIS1, NAXIS2;
	double Z, A, CHI, BOXSIZE, CATCHMENT_close, CATCHMENT_far, H_0, Omega_m, Omega_Lambda, w0, wa, ns, sigma_8, initial_condition, NumPartTotal[6], mass[6], Rot_X, Rot_Y, Rot_Z, Shift_X, Shift_Y, Shift_Z;
	NAXIS=0; NAXIS1=0; NAXIS2=0; Z=0; A=0; CHI=0; BOXSIZE=0; CATCHMENT_close=0; CATCHMENT_far=0; H_0=0; Omega_m=0; Omega_Lambda=0; w0=0; wa=0; ns=0; sigma_8=0; initial_condition=0; Rot_X=0; Rot_Y=0; Rot_Z=0; Shift_X=0; Shift_Y=0; Shift_Z=0;
	char line[1000], parameter_name[1000], equalsign[100];
	int line_length;
	float parameter_value;
	int i;
	
	for(i=0;i<6;i++)
	{
		NumPartTotal[i]=0.0;
		mass[i]=0.0;
	}

	// NOTE: This is the readin of keyword card that should be done better without file:
	FILE *temp_file;
	temp_file=fopen(tempname, "r");
	
	while ((line_length=fgetline(temp_file, line, sizeof(line))) > 0)
	{
		if (sscanf(line, "%s %s %e", parameter_name, equalsign, &parameter_value) == 3)
		{

		if (strcmp(parameter_name,"NAXIS")==0 && strcmp(equalsign,"=")==0)
			{
				NAXIS=(int) parameter_value;
				if (feedback>2) printf("NAXIS detected: %d\n", NAXIS);
			}	
			
		if (strcmp(parameter_name,"NAXIS1")==0 && strcmp(equalsign,"=")==0)
			{
				NAXIS1=(int) parameter_value;
				if (feedback>2) printf("NAXIS1 detected: %d\n", NAXIS1);
			}				

		if (strcmp(parameter_name,"NAXIS2")==0 && strcmp(equalsign,"=")==0)
			{
				NAXIS2=(int) parameter_value;
				if (feedback>2) printf("NAXIS2 detected: %d\n", NAXIS2);
			}	
			
		if (strcmp(parameter_name,"Z")==0 && strcmp(equalsign,"=")==0)
			{
				Z=(double) parameter_value;
				if (feedback>2) printf("Z detected: %e\n", Z);
			}	

		if (strcmp(parameter_name,"A")==0 && strcmp(equalsign,"=")==0)
			{
				A=(double) parameter_value;
				if (feedback>2) printf("A detected: %e\n", A);
			}	

		if (strcmp(parameter_name,"CHI")==0 && strcmp(equalsign,"=")==0)
			{
				CHI=(double) parameter_value;
				if (feedback>2) printf("CHI detected: %e\n", CHI);
			}

		if (strcmp(parameter_name,"BOXSIZE")==0 && strcmp(equalsign,"=")==0)
			{
				BOXSIZE=(double) parameter_value;
				if (feedback>2) printf("BOXSIZE detected: %e\n", BOXSIZE);
			}

                if (strcmp(parameter_name,"CLOSE")==0 && strcmp(equalsign,"=")==0)
		  {
		    CATCHMENT_close=(double) parameter_value;
		    if (feedback>2) printf("Catchment CLOSE detected: %e\n", CATCHMENT_close);
		  }

                if (strcmp(parameter_name,"FAR")==0 && strcmp(equalsign,"=")==0)
		  {
		    CATCHMENT_far=(double) parameter_value;
		    if (feedback>2) printf("Catchment FAR detected: %e\n", CATCHMENT_far);
		  }
		
		if (strcmp(parameter_name,"H_0")==0 && strcmp(equalsign,"=")==0)
			{
				H_0=(double) parameter_value;
				if (feedback>2) printf("H_0 detected: %e\n", H_0);
			}

		if (strcmp(parameter_name,"OMEGA_M")==0 && strcmp(equalsign,"=")==0)
			{
				Omega_m=(double) parameter_value;
				if (feedback>2) printf("Omega_m detected: %e\n", Omega_m);
			}
			
		if (strcmp(parameter_name,"OMEGA_L")==0 && strcmp(equalsign,"=")==0)
			{
				Omega_Lambda=(double) parameter_value;
				if (feedback>2) printf("Omega_Lambda detected: %e\n", Omega_Lambda);
			}

			
        if (strcmp(parameter_name,"W_0")==0 && strcmp(equalsign,"=")==0)
			{
				w0=(double) parameter_value;
				if (feedback>2) printf("w_0 detected: %e\n", w0);
			}
			
        if (strcmp(parameter_name,"W_A")==0 && strcmp(equalsign,"=")==0)
			{
				wa=(double) parameter_value;
				if (feedback>2) printf("w_0 detected: %e\n", wa);
			}

        if (strcmp(parameter_name,"N_S")==0 && strcmp(equalsign,"=")==0)
			{
				ns=(double) parameter_value;
				if (feedback>2) printf("n_s detected: %e\n", ns);
			}

        if (strcmp(parameter_name,"SIGMA_8")==0 && strcmp(equalsign,"=")==0)
			{
				sigma_8=(double) parameter_value;
				if (feedback>2) printf("sigma_8 detected: %e\n", sigma_8);
			}

        if (strcmp(parameter_name,"IC")==0 && strcmp(equalsign,"=")==0)
			{
				initial_condition=(double) parameter_value;
				if (feedback>2) printf("Initial Condition (IC) detected: %e\n", initial_condition);
			}
            
		if (strcmp(parameter_name,"N_CDM")==0 && strcmp(equalsign,"=")==0)
			{
				NumPartTotal[1]=(double) parameter_value;
				if (feedback>2) printf("Number of CDM particles (N_CDM) detected: %e\n", NumPartTotal[1]);
			}

		if (strcmp(parameter_name,"M_CDM")==0 && strcmp(equalsign,"=")==0)
			{
				mass[1]=(double) parameter_value;
				if (feedback>2) printf("Mass of CDM particles (M_CDM) detected (in 10^10 solar masses/h): %e\n", mass[1]);
			}

                if (strcmp(parameter_name,"ROT_X")==0 && strcmp(equalsign,"=")==0)
		  {
		    Rot_X=(double) parameter_value;
		    if (feedback>2) printf("Rot_X detected: %e\n", Rot_X);
		  }

                if (strcmp(parameter_name,"ROT_Y")==0 && strcmp(equalsign,"=")==0)
                  {
                    Rot_Y=(double) parameter_value;
                    if (feedback>2) printf("Rot_Y detected: %e\n", Rot_Y);
                  }

                if (strcmp(parameter_name,"ROT_Z")==0 && strcmp(equalsign,"=")==0)
                  {
                    Rot_Z=(double) parameter_value;
                    if (feedback>2) printf("Rot_Z detected: %e\n", Rot_Z);
                  }

                if (strcmp(parameter_name,"SHIFT_X")==0 && strcmp(equalsign,"=")==0)
                  {
                    Shift_X=(double) parameter_value;
                    if (feedback>2) printf("Shift_X detected: %e\n", Shift_X);
                  }

                if (strcmp(parameter_name,"SHIFT_Y")==0 && strcmp(equalsign,"=")==0)
                  {
                    Shift_Y=(double) parameter_value;
                    if (feedback>2) printf("Shift_Y detected: %e\n", Shift_Y);
                  }

                if (strcmp(parameter_name,"SHIFT_Z")==0 && strcmp(equalsign,"=")==0)
                  {
                    Shift_Z=(double) parameter_value;
                    if (feedback>2) printf("Shift_Z detected: %e\n", Shift_Z);
                  }

		
		}
				
	}
	fclose(temp_file);
	
	FITSheader[plane_number].NAXIS=NAXIS;
	FITSheader[plane_number].NAXIS1=NAXIS1;
	FITSheader[plane_number].NAXIS2=NAXIS2;
	FITSheader[plane_number].Z=Z;
	FITSheader[plane_number].CHI=CHI;
	FITSheader[plane_number].A=A;
	FITSheader[plane_number].BOXSIZE=BOXSIZE;
	FITSheader[plane_number].H_0=H_0;
	FITSheader[plane_number].Omega_m=Omega_m;
	FITSheader[plane_number].Omega_Lambda=Omega_Lambda;
    FITSheader[plane_number].w0=w0;
    FITSheader[plane_number].wa=wa;
    FITSheader[plane_number].ns=ns;
    FITSheader[plane_number].sigma_8=sigma_8;
    FITSheader[plane_number].initial_condition=initial_condition;

	Plane[plane_number].catchment_close=CATCHMENT_close;
	Plane[plane_number].catchment_far=CATCHMENT_far;
	Plane[plane_number].Rot[0]=Rot_X;
	Plane[plane_number].Rot[1]=Rot_Y;
	Plane[plane_number].Rot[2]=Rot_Z;
	Plane[plane_number].Shift[0]=Shift_X;
	Plane[plane_number].Shift[1]=Shift_Y;
	Plane[plane_number].Shift[2]=Shift_Z;

	for (i=0;i<6;i++)
	{
		FITSheader[plane_number].NumPartTotal[i]=NumPartTotal[i];
		FITSheader[plane_number].mass[i]=mass[i];
	}

	printf("Header corrupt? Filename, H_0, Omega_m Omega_Lambda, w0, wa, ns sigma_8, initial_condition:\n%s:\n%e %e %e %e %e %e %e %e %e %e (FITSheader, local)\n", filename,  FITSheader[plane_number].H_0, FITSheader[plane_number].Omega_m, H_0, Omega_m, Omega_Lambda, w0, wa, ns, sigma_8, initial_condition);


    return;
}





void readFITSpotential_singleplane(char filename[], double *potential_array)

     /************************************************************************/
     /* Read a FITS image and determine the minimum and maximum pixel values */
     /************************************************************************/
{
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status,  nfound, anynull;
  long naxes[2], fpixel, nbuffer, npixels, ii;
  
  double *image;
  int i, j, k;

#define buffsize 1000
  float datamin, datamax, nullval, buffer[buffsize];
  //char filename[]  = "atestfil.fit";     /* name of existing FITS file   */

  status = 0;

  if ( fits_open_file(&fptr, filename, READONLY, &status) )
    printerror( status );

  /* read the NAXIS1 and NAXIS2 keyword to get image size */
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    printerror( status );

  npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
  image=(double *) malloc(npixels*sizeof(double));
    assert(image!=NULL);
  fpixel   = 1;
  nullval  = 0;                /* don't check for null values in the image */
  datamin  = 1.0E30;
  datamax  = -1.0E30;

  while (npixels > 0)
    {
      nbuffer = npixels;
      if (npixels > buffsize)
        nbuffer = buffsize;     /* read as many pixels as will fit in buffer */

      /* Note that even though the FITS images contains unsigned integer */
      /* pixel values (or more accurately, signed integer pixels with    */
      /* a bias of 32768),  this routine is reading the values into a    */
      /* float array.   Cfitsio automatically performs the datatype      */
      /* conversion in cases like this.                                  */

      if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
			 buffer, &anynull, &status) )
	printerror( status );

      for (ii = 0; ii < nbuffer; ii++) {
        
	image[fpixel-1+ii]=buffer[ii];
	
	/*
	  if ( buffer[ii] < datamin )
            datamin = buffer[ii];

        if ( buffer[ii] > datamax )
            datamax = buffer[ii];
	*/
	
	//printf("Buffer %e \n", buffer[ii]);
      }
      npixels -= nbuffer;    /* increment remaining number of pixels */
      fpixel  += nbuffer;    /* next pixel to be read in image */
    }

  //printf("\nMin and max image pixels =  %.0f, %.0f\n", datamin, datamax);
  //printf("\nMin and max image pixels =  %e, %e\n", datamin, datamax);

  if ( fits_close_file(fptr, &status) )
    printerror( status );
  
/*
  double random_number1, random_number2;
  int iii, jjj;
  
  if (parameters.plane_shift==1)
    {
      random_number1=ran2(&parameters.seed);
      random_number2=ran2(&parameters.seed);
      if (feedback>2) printf("Random Numbers: %e %e \n", random_number1, random_number2);

      
      for (i=0; i<parameters.ny; i++)
	{
	  for (j=0; j<parameters.nx; j++)
	    {
		  k=i*parameters.nx+j;
	      iii=i+((int) (random_number1*parameters.ny));
	      jjj=j+((int) (random_number2*parameters.nx));
	      if (iii>=parameters.ny) iii=iii-parameters.ny;
	      if (jjj>=parameters.nx) jjj=jjj-parameters.nx;
	      potential_array[iii*parameters.nx+jjj]=image[k];
	    }
	}
    }
  else
    {  
      for (k=0; k<parameters.nxny; k++)
	{
	      potential_array[k]=image[k];
	}
    }
 */
    
    // for (k=0; k<parameters.nxny; k++)
    for (k=0; k<naxes[0]*naxes[1]; k++)
    {
        potential_array[k]=image[k];
    }
 
  free(image);
  return;
}




void printerror(int status)
{
	printf("ERROR! Status value: %d.\n", status);
	exit(1);
}
