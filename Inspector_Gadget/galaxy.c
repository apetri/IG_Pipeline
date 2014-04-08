/*
 *  galaxy.c
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at the University of KwaZulu-Natal on 02/09/2014.
 *  Copyright 2014. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <time.h>

#include "fitsio.h"

#include "main.h"
#include "io_routines.h"
#include "galaxy.h"
#include "2D-plane_multi.h"
#include "fits.h"


#define MAXLINE 1000


int get_number_of_galaxies(char galaxy_input_filename[], int *number_of_redshifts_per_galaxy)
// This function returns the number of galaxies in the catalogue file, as well as the largest number of redshifts for any galaxy in the file (number_of_redshifts_per_galaxy is an output variable).
{
    FILE *galaxy_input_file;
    int number_of_galaxies, galaxy_read;
    struct galaxy_struct galaxy;
    int redshifts_per_galaxy=0;
    *number_of_redshifts_per_galaxy=0;
    
    galaxy_input_file=fopen(galaxy_input_filename, "r");
    
    galaxy_read=1;
    number_of_galaxies=0;
    
    while (galaxy_read==1)
    {
        galaxy_read=get_galaxy_parameters(galaxy_input_file, &galaxy, &redshifts_per_galaxy);
        if (galaxy_read==1) number_of_galaxies++;
        if (redshifts_per_galaxy>*number_of_redshifts_per_galaxy) *number_of_redshifts_per_galaxy=redshifts_per_galaxy;
    }
    
    fclose(galaxy_input_file);
    
    return number_of_galaxies;
}



int get_galaxy_parameters(FILE *galaxy_input_file, struct galaxy_struct *galaxy, int *number_of_redshifts_per_galaxy)
{
    char line[MAXLINE];
    int line_length;
    float parameter_value0, parameter_value1, parameter_value2, parameter_value3;
    float parameter_value[100];
    int i;
    int feedback=0;
    
    int galaxy_read=0;
    
    
        while ((line_length=fgetline(galaxy_input_file, line, sizeof(line))) > 0)
        {
            if (feedback >=2) printf("Reading a line from input file:\nLine length: %d \nLine Content: %s", line_length, line);
				
        
            // WARNING: In the lines below, only up to five redshifts per galaxy have been implemented. Add additional cases if want more redshifts. (This could be implemented better, more universally).
            
            ///////////////////////////////////////////////////////////
            // Reverse order here is important, because the 1-redshift read is satisfied even if there are 2 or more redshifts, etc., so need to start with the largest prospective number of galaxies.
            // For 5 redshifts per galaxy:
            if (sscanf(line, "%e %e %e %e %e %e %e", &parameter_value[1], &parameter_value[2], &parameter_value[3], &parameter_value[4], &parameter_value[5], &parameter_value[6], &parameter_value[7]) == 7)
            {
                *number_of_redshifts_per_galaxy=5;
                galaxy_read=1; // set to 1 if have read a galaxy (regardless of number of redshifts).
            }
            // For 4 redshifts per galaxy:
            else if (sscanf(line, "%e %e %e %e %e %e", &parameter_value[1], &parameter_value[2], &parameter_value[3], &parameter_value[4], &parameter_value[5], &parameter_value[6]) == 6)
            {
                *number_of_redshifts_per_galaxy=4;
                galaxy_read=1; // set to 1 if have read a galaxy (regardless of number of redshifts).
            }
            // For 3 redshifts per galaxy:
            else if (sscanf(line, "%e %e %e %e %e", &parameter_value[1], &parameter_value[2], &parameter_value[3], &parameter_value[4], &parameter_value[5]) == 5)
            {
                *number_of_redshifts_per_galaxy=3;
                galaxy_read=1; // set to 1 if have read a galaxy (regardless of number of redshifts).
            }
            // For 2 redshift per galaxy:
            else if (sscanf(line, "%e %e %e %e", &parameter_value[1], &parameter_value[2], &parameter_value[3], &parameter_value[4]) == 4)
            {
                *number_of_redshifts_per_galaxy=2;
                galaxy_read=1; // set to 1 if have read a galaxy (regardless of number of redshifts).
            }
            // For 1 redshift per galaxy:
            else if (sscanf(line, "%e %e %e", &parameter_value[1], &parameter_value[2], &parameter_value[3]) == 3)
            {
                *number_of_redshifts_per_galaxy=1;
                galaxy_read=1; // set to 1 if have read a galaxy (regardless of number of redshifts).
            }
            ///////////////////////////////////////////////////////////
            
            if (galaxy_read==1)
            {
                (*galaxy).theta[0]=(double) parameter_value[1];
                (*galaxy).theta[1]=(double) parameter_value[2];
                for (i=0; i<*number_of_redshifts_per_galaxy; i++) (*galaxy).redshift[i]=(double) parameter_value[i+3];
        
                
                galaxy_read=1; // set to 1 if have read a galaxy.
                return galaxy_read; // returns value 1 if a galaxy has been read successfully.
            }
        }
    
    // If no galaxy read, zero out part of galaxy structure:
    (*galaxy).theta[0]=0.0;
    (*galaxy).theta[1]=0.0;
    (*galaxy).redshift[0]=0.0; // WARNING: Note that other redshift entries not zeroed out, only the first one for this galaxy.
    *number_of_redshifts_per_galaxy=0;
    
    return galaxy_read; // returns value 0 if no galaxy has been read (also means end of file is reached and there are no more galaxies in the file).
}


/*
void save_WL_galaxy_catalogue(char galaxy_input_filename[], char galaxy_output_filename[], double **theta1, double **theta2, double **theta1_ini, double **theta2_ini, double **A11, double **A12, double **A21, double **A22, double *writeout_array, int NbinsX, int NbinsY)
{
    // NOTE: writeout_array not used in this function.
    
    FILE *galaxy_input_file;
    FILE *galaxy_output_file;
    
    
    int Xbin, Ybin;
    int galaxy_number, galaxy_read;
    struct galaxy_struct galaxy;
    
    galaxy_number=0;
    galaxy_read=1;
    
    
    galaxy_input_file=fopen(galaxy_input_filename, "r");
    galaxy_output_file=fopen(galaxy_output_filename, "w");
    
    
    for (Ybin=0; Ybin<NbinsY; Ybin++)
    {
        for (Xbin=0; Xbin<NbinsX; Xbin++)
        {
            if (galaxy_read>0)
            {
                galaxy_read=get_galaxy_parameters(galaxy_input_file, &galaxy);
    
                if (galaxy_read>0)
                {
                    galaxy_number++; // galaxy numbering in output catalogue starts at 1.
                    
                    // kappa:
                    galaxy.kappa=-(A11[Ybin][Xbin]+A22[Ybin][Xbin])/2.0+1.0;
                    // gamma1:
                    galaxy.gamma[0]=-(A11[Ybin][Xbin]-A22[Ybin][Xbin])/2.0;
                    // gamma2:
                    galaxy.gamma[1]=-(A12[Ybin][Xbin]+A21[Ybin][Xbin])/2.0;
                    // omega:
                    galaxy.omega=-(A12[Ybin][Xbin]-A21[Ybin][Xbin])/2.0;
            
 
                    fprintf(galaxy_output_file, "%d %e %e %e %e %e %e\n", galaxy_number, galaxy.theta[0], galaxy.theta[1], galaxy.redshift, galaxy.kappa, galaxy.gamma[0], galaxy.gamma[1]);
                    
                }
            }
        }
    }
    
    fclose(galaxy_input_file);
    fclose(galaxy_output_file);
            
}
 
 */

void save_WL_galaxy_catalogue_output_only(char galaxy_output_filename[], double **theta1, double **theta2, double **theta1_ini, double **theta2_ini, double **A11, double **A12, double **A21, double **A22, double *writeout_array, int last_plane_number, int NbinsX, int NbinsY, int number_of_galaxies, int number_of_redshifts_per_galaxy)
{
    // NOTE: writeout_array not used in this function.
    

    FILE *galaxy_output_file;

    int i, j;
    int Xbin, Ybin;
    int galaxy_number, galaxy_read;
    struct galaxy_struct galaxy;
    double *output_array;
    long mapdims[2];
    
    mapdims[0]=number_of_redshifts_per_galaxy*4;
    mapdims[1]=number_of_galaxies;
    
    output_array=malloc(mapdims[0]*mapdims[1]*sizeof(double));
    
    galaxy_number=0;
    i=0;
    j=0;
    
    for (Ybin=0; Ybin<NbinsY; Ybin++)
    {
        for (Xbin=0; Xbin<NbinsX; Xbin++)
        {
            
            // kappa:
            galaxy.kappa[i]=-(A11[Ybin][Xbin]+A22[Ybin][Xbin])/2.0+1.0;
            // gamma1:
            galaxy.gamma1[i]=-(A11[Ybin][Xbin]-A22[Ybin][Xbin])/2.0;
            // gamma2:
            galaxy.gamma2[i]=-(A12[Ybin][Xbin]+A21[Ybin][Xbin])/2.0;
            // omega:
            galaxy.omega[i]=-(A12[Ybin][Xbin]-A21[Ybin][Xbin])/2.0;
            
            if (galaxy_number<number_of_galaxies)
            {
                output_array[j]=galaxy.kappa[i];
                output_array[j+1]=galaxy.gamma1[i];
                output_array[j+2]=galaxy.gamma2[i];
                output_array[j+3]=galaxy.omega[i];
                // printf("Output: %e %e %e %e\n", output_array[j], output_array[j+1], output_array[j+2], output_array[j+3]);
                j=j+4;
                i++;
                if (i==number_of_redshifts_per_galaxy)
                {
                    galaxy_number++;
                    i=0;
                }
            }
            
            
         }
    }
    
    writeWLmapFITSimage_f(galaxy_output_filename, 2, mapdims, output_array,last_plane_number);
    
}



#undef MAXLINE
