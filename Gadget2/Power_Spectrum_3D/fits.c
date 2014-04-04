/*
 *  fits.c
 *  FITS-Modifier
 *
 *  Created by Jan Michael Kratochvil on 1/26/09 at Columbia University.
 *  Copyright 2009 Jan Michael Kratochvil. All rights reserved.
 *
 */





#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
//#include <assert.h>

//#include <complex.h>
//#include <fftw3.h>
#include "fitsio.h"
//#include "main.h"
#include "fits.h"



/*
int main (int argc, const char * argv[]) {
    // insert code here...


	char filename[200], filename2[200];
	double *imagearray, min, max, nu;
	int nx, ny, nxny, i;
	// Image dimensions:
	nx=2048; ny=2048; // The sample map which is read in and written out by this example program is 2048x2048 pixels. In general, it is recommended to read in the header of the FITS file first and determine the pixel dimensions of the image automatically. This simple sample program does not have header reading capability though.
	long mapdims[2]={2048, 2048};
	nxny=nx*ny;
	imagearray=malloc(nxny*sizeof(double));
	
	nu=0.05; // Threshold above which to color area (for Minkowski functionals). This is just used as an example here: the code reads in a weak lensing map (FITS file), sets all pixels above this threshold to zero and the ones below to 1, and the writes out the modified map into a new FITS file.

	//sprintf(filename, "WL-conv_j512w08_4_2048xy_0003r_0054p_0200z_og.fit");
	sprintf(filename, "WL-conv_j512w10_4_2048xy_0003r_0054p_0200z_0060s.fit");
	readFITSimage_f(filename, imagearray);

	/////////////////////////////////////////////////////////
	// Modify image here by modifying elements of imagearray:
	
	min=1e100; max=-1e100;

	for (i=0; i<nxny-1; i++)
	{
		if (min>imagearray[i]) min=imagearray[i];
		if (max<imagearray[i]) max=imagearray[i];
	}
	printf("Min and max pixel values in FITS image: %e %e\n", min, max);

	for (i=0; i<nxny-1; i++)
	{
		if(imagearray[i]>=nu) imagearray[i]=0.0;
		else imagearray[i]=1.0;
	}

	
	/////////////////////////////////////////////////////////

	sprintf(filename2, "Modified_%s", filename);
	printf("%s\n", filename2);
	
	writeFITSimage_f(filename2, 2, mapdims, imagearray);



    printf("Done with FITS image modification. Program ran successfully.\n");
    return 0;
}
*/






void readFITSimage_f(char filename[], double *imagearray)

    /************************************************************************/
    /* Read a FITS image and determine the minimum and maximum pixel values */
    /************************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status,  nfound, anynull;
    long naxes[2], fpixel, nbuffer, npixels, ii;
	int npixels2;
	
	double *imagecore;
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
	npixels2  = ((int) naxes[0] * naxes[1]);
	printf("npixels %d\n", npixels);
	imagecore=(double *) malloc(npixels*sizeof(double));
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

      for (ii = 0; ii < nbuffer; ii++)  {
        
		imagecore[fpixel-1+ii]=buffer[ii];
		
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
		
	for (i=0; i<npixels2; i++)
	{
		imagearray[i]=imagecore[i];
	}
	
	free(imagecore);
    return;
}



void writeFITSimage_f(char filename[], long naxis, long naxes[], double *imagearray) // ATTENTION: externally fed in image imagearray must be a 1D-(pseudo-2D) array.

    /******************************************************/
    /* Create a FITS primary array containing a 2-D image */
    /******************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj;
    long  fpixel, nelements;
    long exposure;
	float Z_redshift, A_scale_factor, CHI_comoving_distance, BOXSIZE_boxsize;
    //unsigned short *array[200];
	// float **array;
	// array = (float **) malloc(naxes[1]*sizeof (float));

	float *writearray;
	writearray = (float *) malloc(naxes[0]*naxes[1]*sizeof(float));

    /* initialize FITS image parameters */
    //char filename[] = "atestfil.fit";             /* name for new FITS file */
    int bitpix   =  FLOAT_IMG; /* 32-bit(?) signed float pixel values       */
    //long naxis    =   2;  /* 2-dimensional image                            */    
    //long naxes[2] = { 300, 200 };   /* image is 300 pixels wide by 200 rows */

    /* allocate memory for the whole image */ 
    // array[0] = (float *)malloc( naxes[0] * naxes[1]
    //                             * sizeof( float ) );

    /* initialize pointers to the start of each row of the image */
    //for( ii=1; ii<naxes[1]; ii++ )
    //  array[ii] = array[ii-1] + naxes[0];

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

    /* initialize the values in the image with a linear ramp function */
    /*
    for (jj = 0; jj < naxes[1]; jj++)
    {   for (ii = 0; ii < naxes[0]; ii++)
        {
            // array[jj][ii] = imagearray[ii][naxes[1]-1-jj];  // here is where one reads in an external file <JMK>.
			// array[jj][ii] = imagearray[ii][jj];  // here is where one reads in an external file <JMK>.
			array[jj][ii]=imagearray[jj*naxes[1]+ii];
        }
    }
    */

    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

	for (ii=0; ii<nelements; ii++)
	{
		writearray[ii]=((float) imagearray[ii]);
	}

    /* write the array of unsigned integers to the FITS file */
	    if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, writearray, &status) )
	      //    if ( fits_write_img(fptr, TFLOAT, fpixel, nelements, array[0], &status) )
        printerror( status );
      
	    //    free( array[0] );  /* free previously allocated memory */
	    free(writearray);

    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
      
	 exposure = 1500.;
    if ( fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
         "Total Exposure Time", &status) )
         printerror( status );           

	// Done writing header information into FITS file.


    if ( fits_close_file(fptr, &status) )                /* close the file */
         printerror( status );           

    return;
}


void printerror(int status)
{
	printf("ERROR! Status value: %d.\n", status);
	exit(1);
}
