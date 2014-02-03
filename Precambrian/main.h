/*
 *  main.h
 *  Precambrian (Cosmological N-body Simulation Initial Conditions Organizer) 
 *
 *  Created by Jan Michael Kratochvil - on 01/05/10 - at the University of Miami, Coral Gables, Florida, USA. 
 *  Copyright 2010 Jan Michael Kratochvil. All rights reserved.
 *
 */


// Must be included in all .c files that are not main.c
// Contains all references to global variables and functions which are declared in main.c

extern double vel_prefac_lam;


// Global path and folder variables:
extern char CAMB_dir[1000], NGenIC_dir_L[1000], NGenIC_dir_P[1000], Gadget_dir_L[1000], Gadget_dir_P[1000]; // full paths to applications (CAMB on power spectrum generation computer (LSST), N-GenIC and Gadget (on Blue Gene/L and Blue Gene/P, two separate versions) on simulation computer). For Gadget there can be several alternatives: with w, Octopus-w, noGSL-w, and Octopus-noGSL-w.
extern char CAMB_exec[200], NGenIC_exec[200], Gadget_exec[200]; // names of executables at above paths.              
extern char ics_dir_speccomp[1000], ics_dir_simcomp[1000], ics_data_dir_simcomp[1000], Gadget_output_dir[1000]; // full master paths on power spectrum generation computer cluster (LSST) and simulation computer cluster (NYBlue, private parameter and public data path), as well as location of simulation output (snapshots) of Gadget-2 on simulation computer cluster (NYBlue, public).
extern char series[200], series_folder[200]; // N-body simulation series name (short identifier, typically one letter, e.g. "m") and folder name for this simulation series (e.g. "m-series").                                    
extern char CAMB_folder[200], NGenIC_folder[200], Gadget_folder[200]; // names of folders for supportive data for CAMB, N-GenIC, and Gadget-2 (parameter files, input and output data, job description files for Condor and NYBlue, and log files from runs.
extern char CAMB_param_folder[200], CAMB_data_folder[200], CAMB_jobs_folder[200], CAMB_logs_folder[200]; // folder for CAMB parameter files, output power spectra, Condor submission script, and logs from Condor run.            
extern char NGenIC_param_folder[200], NGenIC_data_folder[200], NGenIC_jobs_folder[200], NGenIC_logs_folder[200]; // folder for N-GenIC parameter files, input power spectra (converted from CAMB output), NYBlue job submission scripts, and logs from runs on NYBlue (<number>.out and <number>.err).
extern char Gadget_param_folder[200], Gadget_data_folder[200], Gadget_jobs_folder[200], Gadget_logs_folder[200]; // folder for Gadget-2 parameter files, initial condition data files (output of N-GenIC), NYBlue job submission scripts, and logs from runs on NYBlue.
extern char IG_output_dir[1000], IG_planes_folder[200], IG_maps_folder[200], IG_products_nonoise_folder[200], IG_products_noise_folder[200];
extern char home_path[1000], repository_path[1000], mass_storage_path[1000];


// Functions:
void convert_CAMB_power_spectrum(char power_spectrum_filename[], char converted_power_spectrum_filename[]);
void write_CAMB_parameter_file(char CAMB_param_filename[], char filebase[], char power_end[], double OBh2, double OCh2, double OM, double OL, double OK, double w0, double wa, double ns, double As, double h, double z);
void write_CAMB_condor_job_description(FILE* CAMB_condor_file, char CAMB_param_filename[], char filebase[]);
void write_NGenIC_parameter_file(char converted_power_spectrum_filename[], char NGenIC_param_filename[], char filebase2[], int part, double boxsize, double OBh2, double OCh2, double OM, double OL, double OK, double w0, double wa, double ns, double s8, double h, double z, int seed, int power_spectrum_at_zini, double Dplus, double vel_prefac_lam);
void write_Gadget_parameter_file(char Gadget_param_filename[], char simulation_codename[], char filebase2[], char output_list_filename[], double boxsize, double OBh2, double OM, double OL, double w0, double wa, double h, double z, double soft);
void write_BGL_description(char NGenIC_param_filename[], char Gadget_param_filename[], char jobstamm[], int job_nr);
void write_BGP_description(char **BGP_NGenIC_param_filename, char **BGP_Gadget_param_filename, char jobstamm[], int job_nr, int octo_counter);
void write_submission_script_L(int jobs, char jobstamm[], char BG_type[]);
void write_submission_script_P(int jobs, int Pjobs, char jobstamm[], char BG_type[]);
void write_directory_script(FILE *script_file, char simulation_codename[]);
void write_IG_directory_script(FILE *script_file, char simulation_codename[]);
void write_simple_list(FILE *script_file, char simulation_codename[]);


