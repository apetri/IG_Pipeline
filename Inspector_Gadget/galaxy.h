/*
 *  galaxy.h
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at the University of KwaZulu-Natal on 02/09/2014.
 *  Copyright 2014. All rights reserved.
 *
 */




struct galaxy_struct
{
    double theta[2]; // anglular position in map (flat sky), must be in radians (not degrees).
    double redshift[100]; // redshift of galaxy (is an array, because it is coded such that each galaxy can have multiple redshifts, drawn from a distribution).
    double kappa[100]; // convergence (is an array for the same reasons as above; each redshift wil produce a different value of convergence).
    double gamma1[100]; // shear (first component)
    double gamma2[100]; // shear (second component)
    double omega[100]; // rotation (a rarely used weak lensing quantity, included mostly for completeness).
    
    // WARNING: Make sure the above array lengths are longer than the number_of_galaxies_per_redshift in the parameter file.
};


int get_number_of_galaxies(char galaxy_input_filename[], int *number_of_galaxies_per_redshift);
int get_galaxy_parameters(FILE *galaxy_input_file, struct galaxy_struct *galaxy, int *number_of_galaxies_per_redshift);
// void save_WL_galaxy_catalogue(char galaxy_input_filename[], char galaxy_output_filename[], double **theta1, double **theta2, double **theta1_ini, double **theta2_ini, double **A11, double **A12, double **A21, double **A22, double *writeout_array, int NbinsX, int NbinsY);
void save_WL_galaxy_catalogue_output_only(char galaxy_output_filename[], double **theta1, double **theta2, double **theta1_ini, double **theta2_ini, double **A11, double **A12, double **A21, double **A22, double *writeout_array, int last_plane_number, int NbinsX, int NbinsY, int number_of_galaxies, int number_of_redshifts_per_galaxy);

