/*
*
* Precambrian. Copyright 2010 Jan Michael Kratochvil at the University of Miami.
* Sets up parameter files for CAMB, N-GenIC, and Gadget2, as part of the Inspector Gadget Weak Gravitational Lensing Simulation Pipeline.
* This code calls (includes) Lam Hui's FORTRAN 77 growth factor code, which computes the growth factor for w(z) for the modified N-GenIC initial conditions generator parameter file.
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "main.h"
#include "darkenergy.h"
#include "interface.h"

#include "options.h"
#include "ini.h"



double vel_prefac_lam; // global variable for velocity prefactor for particles, written into by interface function to Lam Hui's growth factor code (the code returns the growth factor as its return value).

// Global path and folder variables:
char CAMB_dir[1000], NGenIC_dir_L[1000], NGenIC_dir_P[1000], Gadget_dir_L[1000], Gadget_dir_P[1000]; // full paths to applications (CAMB on power spectrum generation computer (LSST), N-GenIC and Gadget (on Blue Gene/L and Blue Gene/P, two separate versions) on simulation computer). For Gadget there can be several alternatives: with w, Octopus-w, noGSL-w, and Octopus-noGSL-w. 
char CAMB_exec[200], NGenIC_exec[200], Gadget_exec[200]; // names of executables at above paths.
char ics_dir_speccomp[1000], ics_dir_simcomp[1000], ics_data_dir_simcomp[1000], Gadget_output_dir[1000]; // full master paths on power spectrum generation computer cluster (LSST) and simulation computer cluster (NYBlue, private parameter and public data path), as well as location of simulation output (snapshots) of Gadget-2 on simulation computer cluster (NYBlue, public).
char series[200], series_folder[200]; // N-body simulation series name (short identifier, typically one letter, e.g. "m") and folder name for this simulation series (e.g. "m-series").
char CAMB_folder[200], NGenIC_folder[200], Gadget_folder[200]; // names of folders for supportive data for CAMB, N-GenIC, and Gadget-2 (parameter files, input and output data, job description files for Condor and NYBlue, and log files from runs.
char CAMB_param_folder[200], CAMB_data_folder[200], CAMB_jobs_folder[200], CAMB_logs_folder[200]; // folder for CAMB parameter files, output power spectra, Condor submission script, and logs from Condor run.
char NGenIC_param_folder[200], NGenIC_data_folder[200], NGenIC_jobs_folder[200], NGenIC_logs_folder[200]; // folder for N-GenIC parameter files, input power spectra (converted from CAMB output), NYBlue job submission scripts, and logs from runs on NYBlue (<number>.out and <number>.err).
char Gadget_param_folder[200], Gadget_data_folder[200], Gadget_jobs_folder[200], Gadget_logs_folder[200]; // folder for Gadget-2 parameter files, initial condition data files (output of N-GenIC), NYBlue job submission scripts, and logs from runs on NYBlue.
char IG_output_dir[1000], IG_planes_folder[200], IG_maps_folder[200], IG_products_nonoise_folder[200], IG_products_noise_folder[200];

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
int fgetline(FILE *stream, char s[], int lim);



int main (int argc, const char * argv[]) {

char CAMB_param_filename[200], NGenIC_param_filename[200], Gadget_param_filename[200], power_end[200], power_front[200], ics_front[200], filerawbase[200], filebase[200], simulation_codename[200], filebase2[200], power_spectrum_filename[200], converted_power_spectrum_filename[200];	
// Strings above contain filenames and file name fragments:
// filerawbase: parameters necessary for power spectrum (missing number of particles, boxsize, sigma_8).
// filebase: filerawbase with series name prepended (used for power spectrum files for CAMB and N-GenIC, where normalization, random seed, particle number and boxsize do not matter yet).
// simulation_codename: Full simulation identifying name (with all required parameters): series name, number of particles, boxsize, power spectrum parameters (from filerawbase), sigma_8 and random seed index added.
// filebase2: ics_front prepended to simulation codename to identify that it is the initial conditions for that simulation (used for parameter file and output of N-GenIC). 
// ics_front: ("ics" by default) identifies that file is initial conditions or parameter file for N-GenIC.
// power_front: prepended string in front of converted power spectrum file name to identify that it has been converted for input into N-GenIC ("MPS_N-GenIC").
// power_end: ending of power spectrum output file of CAMB ("... .dat").
// CAMB_param_filename: fully constructed parameter filename for CAMB ("param... .ini").
// NGenIC_param_filename: fully constructed parameter filename for N-GenIC ("ics_... .param").
int i;
char **BGP_NGenIC_param_filename, **BGP_Gadget_param_filename; // same as above but for Octopus on Blue Gene/P (array of 8 names).
BGP_NGenIC_param_filename=(char **)malloc(8*sizeof(char *));
BGP_Gadget_param_filename=(char **)malloc(8*sizeof(char *));
for (i=0; i<8; i++)
{
	BGP_NGenIC_param_filename[i]=(char *)malloc(200*sizeof(char));
	BGP_Gadget_param_filename[i]=(char *)malloc(200*sizeof(char));	
}
char output_list_filename[200]; // name of input file for Gadget-2 which contains the scale factor (time) values at which simulation snapshot outputs are to be generated.

int Nobh2, Nom, Nol, Nw0, Nwa, Nns, Nas, Ns8, Nh, Nz, Nseed; // Number of parameter values for each of these parameters. 
double *OBh2, *OM, *OL, *w0, *wa, *ns, *As, *s8, *h, *z; // Arrays containing the different parameter values for the cosmological parameters.
int *seed; // Array of random number seeds for N-GenIC to realize different random initial conditions for a simulation with same power spectrum and cosmological parameters.
double OK, OCh2, redshift, Dplus; // Derived cosmological parameters (not user seletable, as they follow from the above and are automatically computed).
int i_OBh2, i_OM, i_OL, i_w0, i_wa, i_ns, i_As, i_s8, i_h, i_z, i_seed; // Index counting through values of parameters.
int mode, submission_style, power_spectrum_at_zini, flat_universe, remove; // Flags to select different modes and options.
int part, Nboxsize, i_boxsize; // simulation specific parameters.
double *boxsize; // Array containing box sizes for simulations (more than one size if want to do telescopic weak lensing study). 
 double soft; // gravitational softening length for CDM simulation (needed for Gadget-2 parameter file).
// for Condor:
int process_number, jobs, Pjobs, octo_counter; // counts through number of processes (not used much).
char command[200], job_description_filename[1000], jobstammL[200], jobstammP[200], BG_type[200], directory_script_filename[1000], IG_directory_script_filename[1000], simple_list_filename[1000]; // string containing Condor job description filename and NYBlue job description filenames.
FILE *CAMB_condor_file, *directory_script_file, *IG_directory_script_file, *simple_list_file; // pointer to Condor job description file.
char home_path[1000], repository_path[1000], mass_storage_path[1000];

// Set these parameters to span full suite of N-body simulations. Every combination will be evaluated.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Usage Instructions:
// 0. Compile Precambrian (this code). Stand-alone file is with "gcc main.c -o Precambrian -lm", but if want to include growth factor, need to use Makefile and version which has Lam Hui's FORTRAN 77 growth factor code incorporated.
// 1. Run mode=1 (prepares CAMB).
// 2. Run Condor script which runs CAMB on LSST cluster (this is obsolete, LSST/Astro cluster at BNL does not use Condor anymore).
// 3. Run mode=2 (converts CAMB output for ingestion by N-GenIC).
// 4. Run mode=3 (generates parameter files for N-GenIC and Gadget-2 (all cosmological parameter combinations)
// 5. Run mode=4 (selects which combinations will be run, and generates Blue Gene job description files).
// 6. Submit N-GenIC to Blue Gene
// 7. Submit Gadget-2 to Blue Gene
// 8. Optionally repeat steps 4.--6. if want to add different box sizes later (which were not included in first pass).
// Note: Steps 1--3 need to be rerun only if cosmological parameters change. If power spectrum is being scaled back from z=0, for changes in sigma_8 and z_ini, it is sufficient to rerun Steps 4--6.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Execution modes are directly chosen by providing an argument to the executable. Additional parameters are selectable in this file below.


if (argc < 3)
{
	printf("Run make_directories.py before running Precambrian!!!\n\n");
	printf("Usage: ./Precambrian <ini_options_file> <mode>\nwhere <mode> can equal:\n\n");
	printf("1: Generate parameter files for CAMB and Condor job description file for CAMB execution.\n");
	printf("2: Convert CAMB matter power spectra to N-GenIC power spectra.\n");
	printf("3: Generate N-GenIC and Gadget-2 parameter files (all combinations).\n");
	printf("4: Generate selected NYBlue job description files and submission shell scripts;\n   this option also takes into account the submission_style variable in the code.\n");
	printf("Aborting. Rerun Precambrian with one of the above modes as its argument.\n\n");
	exit(1);
}

//Parse options
sys_options *options = malloc(sizeof(sys_options));

if(ini_parse(argv[1],handler,options)<0){
	fprintf(stderr,"ini options file %s not found\n",argv[1]);
	exit(1);
}

// Mode: 
mode=atoi(argv[2]); // First (and only) argument of Precambrian is mode; select 1, 2, 3 or 4.
printf("\nRunning Precambrian in Mode=%d\n", mode);

// Safety for NYBlue:
/*
if (mode==1)
{
  printf("Mode 1 is not available on NYBlue, because CAMB not set up to run here.\n");
  printf("Use preferrably only Mode 4 here, to select which jobs to run, all the rest should have been set up already.\n"); 
  printf("But Modes 2 (convert power spectrum) and Mode 3 (set up parameter files for N-GenIC and Gadget-2), are also available.\n");
  exit(1);
}
*/
sleep(1);

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// Adjustable Parameters:
/////////////////////////////////////////////////////////////////////

// Submission style (NYBlue job submission files are generated only for desired jobs and numbered sequentially for easy submission later);
submission_style=options->submission_style; // select 1, 2 (recommended default), 3, 4, 5, or 6 here. Affects only mode=4 (submission script generation) and mode=6 (simple simulation codename list output).
// 1: all, 2: cross jobs only, 3: fiducial model only (many random seeds), 4: anti-cross (complement of 2), 5: cross-anti-fiducial (like 2 but without 3), 6: anti-fiducial (complement of 3).
// Standard submissions are: NYBlue/L: single, nonGSL, long jobs on 128 nodes in VN mode (all CPUs), 256 CPUs, wall time 72 hours.
//                           NYBlue/P: Octopus, GSL, normal jobs on 512 nodes in VN mode (all CPUs), 8 simulations in parallel, 256 CPUs each, wall time 48 hours.

///////////////////////////////                                                
// Home, Repository, and Mass Storage Paths (set here depending on Blue Gene type):
// sprintf(home_path, "/gpfs/home2/jank"); // for Blue Gene/L and /P
sprintf(home_path, "%s",options->home_path); // for Blue Gene/Q
sprintf(repository_path, "%s%s", home_path,options->repository_relative_path); // path of Inspector Gadget pipeline repository in home directory.
sprintf(mass_storage_path, "%s",options->mass_storage_path); // path on mass storage disk (for storage of bulky stuff like simulations).  

// N-body simulation specs:
sprintf(series,"%s",options->series_name); // simulation series name (typically one lower case letter).
part=options->num_particles_side; // Number of particles in one dimension, N-body simulation has part^3 particles.
// List of Box sizes (in Mpc/h) for N-body simulation:
Nboxsize=options->Nboxsize; // number of different boxsizes to be evaluated (more can be added later easily).
boxsize=(double *)malloc(Nboxsize*sizeof(double));
for(i=0;i<Nboxsize;i++){
	boxsize[i]=options->boxsize[i];
}
/*
boxsize[0]=480.0;
boxsize[1]=400.0;
boxsize[2]=320.0;
boxsize[3]=240.0;
boxsize[4]=160.0;
boxsize[5]=80.0;
*/

// Settings:
power_spectrum_at_zini=options->power_spectrum_at_zini; // Set !=0 if want to generate power spectrum at initial redshift rather than scaling it back in N-GenIC.
flat_universe=options->flat_universe; // Set to !=0 if want universe to be flat (no curvature); ignores settings of OL (Omega_Lambda) and creates flat universe based on OM (Omega_matter).
remove=options->remove_old; // Set !=0 if want to remove old job files (old ones will be overwritten, but there may be stray superfluous ones so it's recommended even though job submission shell scripts will be updated such that they ignore the superfluous ones).

/////////////////////////////////////////////////////////////////////////////
// COSMOLOGICAL PARAMETERS: (arrays, can run various combinations)
////////////////////////////

// OMEGA BARYON: Fractional baryon density * h^2:
Nobh2=options->Nobh2; // Number of different parameter values of OMEGA BARYON to be investigated; make sure number equals number of different parameter values below.
OBh2=(double *)malloc(Nobh2*sizeof(double));
for(i=0;i<Nobh2;i++){
	OBh2[i]=options->OBh2[i];
}  // OB ~ 0.042;   // OB=0.0437885802469 

// OMEGA MATTER: Fractional total matter density today (CDM + baryons, has to add up to 1 with OL below for flat universe):
Nom=options->Nom; // Number of different parameter values of OMEGA MATTER to be investigated; make sure number equals number of different parameter values below.
OM=(double *)malloc(Nom*sizeof(double));
for(i=0;i<Nom;i++){
	OM[i]=options->OM[i];
}

// OMEGA DARK ENERGY: Fractional dark energy density today (has to add up to 1 with OM above for flat universe):
Nol=options->Nol; // make sure number equals number of different parameter values below.
OL=(double *)malloc(Nol*sizeof(double));
for(i=0;i<Nol;i++){
	OL[i]=options->OL[i];
} 
// will be reset to value to make universe flat if flat_universe flag above is set.

/////////////////////////////////
// DARK ENERGY EQUATION OF STATE:
/////////////////////////////////
// Dark Energy Model: w(z)=w_0+(z/(z+1))*w_a.
// Currently only w_0 works (constant w) with CAMB.
Nw0=options->Nw0; // make sure number equals number of different parameter values below.
w0=(double *)malloc(Nw0*sizeof(double));
for(i=0;i<Nw0;i++){
	w0[i]=options->w0[i];
}

Nwa=options->Nwa; // make sure number equals number of different parameter values below.
wa=(double *)malloc(Nwa*sizeof(double));
for(i=0;i<Nwa;i++){
	wa[i]=options->wa[i];
}
/////////////////////////////////

// Scalar spectral index n_s:
Nns=options->Nns; // make sure number equals number of different parameter values below.
ns=(double *)malloc(Nns*sizeof(double));
for(i=0;i<Nns;i++){
	ns[i]=options->ns[i];
}
/*
ns[1]=0.92;
ns[2]=1.00;
*/

// Primordial amplitude of density perturbations A_s (Note: depends on pivot scale, set in CAMB parameter file):
// (This parameter is reset by sigma_8 normalization in postprocessing.)
Nas=options->Nas; // make sure number equals number of different parameter values below.
As=(double *)malloc(Nas*sizeof(double));
for(i=0;i<Nas;i++){
	As[i]=options->as[i];
}

// sigma_8:
Ns8=options->Ns8; // make sure number equals number of different parameter values below.
s8=(double *)malloc(Ns8*sizeof(double));
for(i=0;i<Ns8;i++){
	s8[i]=options->s8[i];  // m-series was: 0.798;     // 0.79841924
}

// Hubble parameter h: H_0 = 100 * h km/s/Mpc.
Nh=options->Nh; // make sure number equals number of different parameter values below.
h=(double *)malloc(Nh*sizeof(double));
for(i=0;i<Nh;i++){
	h[i]=options->h[i];
}

// Starting Redshift of N-body simulations:
// (Power spectrum redshift can be set via power_spectrum_at_zini flag above to initial redshift of simulation or to z=0 and normalized to sigma_8 today by modified N-GenIC, which is capable of scaling back with dark energy with w(z).)
Nz=options->Nz; // make sure number equals number of different parameter values below.
z=(double *)malloc(Nz*sizeof(double));
for(i=0;i<Nz;i++){
	z[i]=options->z[i];
}

// Random number seed for N-GenIC:
Nseed=options->Nseed; // make sure number equals number of different parameter values below.
seed=(int *)malloc(Nseed*sizeof(int));
for(i=0;i<Nseed;i++){
	seed[i]=options->seed[i];
}

/*
seed[5]=194340;
seed[6]=705031;
seed[7]=674951;
seed[8]=495306;
seed[9]=105884;
seed[10]=932155;
seed[11]=758004;
seed[12]=521745;
seed[13]=610028;
seed[14]=493763;
seed[15]=817290;
seed[16]=715483;
seed[17]=219569;
seed[18]=563201;
seed[19]=695264;
seed[20]=124485;
seed[21]=545083;
seed[22]=814251;
seed[23]=173353;
seed[24]=614519;
seed[25]=840351;
seed[26]=308614;
seed[27]=247316;
seed[28]=782067;
seed[29]=651665;
seed[30]=956032;
seed[31]=943607;
seed[32]=546898;
seed[33]=570575;
seed[34]=871385;
seed[35]=600617;
seed[36]=381579;
seed[37]=225985;
seed[38]=493429;
seed[39]=207875;
seed[40]=370995;
seed[41]=481012;
seed[42]=368944;
seed[43]=662792;
seed[44]=362064;
seed[45]=647263;
seed[46]=195591;
seed[47]=482732;
seed[48]=521713;
seed[49]=855138;
*/

free(options);

/////////////////////////////////////////////////////////////
// Paths, folders, and names (for description of all the variables here, see their declaration above):
// Application Paths:
sprintf(CAMB_dir, "%s/camb", repository_path);
sprintf(NGenIC_dir_L, "%s/N-GenIC", repository_path); // path to version for Blue Gene/L
sprintf(NGenIC_dir_P, "%s/N-GenIC", repository_path); // path to version for Blue Gene/P
sprintf(Gadget_dir_L, "%s/Gadget2", repository_path);
sprintf(Gadget_dir_P, "%s/Gadget2", repository_path);
// Executables:
sprintf(CAMB_exec, "camb");
sprintf(NGenIC_exec, "N-GenIC");
sprintf(Gadget_exec, "Gadget2");
// Master directories for initial conditions:
sprintf(ics_data_dir_simcomp, "%s/Storage/sims/ics", mass_storage_path); // public path for large data.
sprintf(ics_dir_simcomp, "%s/localStorage/ics", repository_path); // private path for parameters (small data).
sprintf(ics_dir_speccomp, "%s", ics_dir_simcomp); // For LSST cluster: sprintf(ics_dir_speccomp, "/data/jank/Storage/ics");
sprintf(Gadget_output_dir, "%s/Storage/sims/snapshots", mass_storage_path);

sprintf(IG_output_dir, "%s/Storage/wl/IG", mass_storage_path); // output directory for Inspector Gadget (weak lensing map generation).
// Folders:
sprintf(series_folder, "%s-series", series);
sprintf(CAMB_folder, "data_CAMB");
sprintf(CAMB_param_folder, "Parameters");
sprintf(CAMB_data_folder, "Output_Data");
sprintf(CAMB_jobs_folder, "Jobs");
sprintf(CAMB_logs_folder, "Logs");
sprintf(NGenIC_folder, "data_N-GenIC");
sprintf(NGenIC_param_folder, "Parameters");
sprintf(NGenIC_data_folder, "Power_Spectra");
sprintf(NGenIC_jobs_folder, "Jobs");
sprintf(NGenIC_logs_folder, "Logs");
sprintf(Gadget_folder, "data_Gadget");
sprintf(Gadget_param_folder, "Parameters");
sprintf(Gadget_data_folder, "IC_Files");
sprintf(Gadget_jobs_folder, "Jobs");
sprintf(Gadget_logs_folder, "Logs");

sprintf(IG_planes_folder, "Planes");
sprintf(IG_maps_folder, "Maps");
sprintf(IG_products_nonoise_folder, "Products_nonoise");
sprintf(IG_products_noise_folder, "Products_noise");
// File name elements:
sprintf(output_list_filename, "outputs_%s-series.txt", series); // filename containing scale factor (time) list when Gadget-2 snapshots outputs are to be made. 
sprintf(power_end, "CAMB-MPS.dat"); // Filename ending for CAMB total matter power spectrum file (will be prepended by filebase).
sprintf(power_front, "MPS-N-GenIC"); // Filename beginning for converted matter power spectrum (input file for N-GenIC).
sprintf(ics_front, "ics");

sprintf(jobstammL, "jobsubmitL_%s", series_folder);
sprintf(jobstammP, "jobsubmitP_%s", series_folder);

/////////////////////////////////////////

printf("Read in list of parameters.\n");
fflush(stdout);

////////////////////////
// Do not modify anything below here unless you have specific reason to.
////////////////////////

//////////////////////////////////////////////////
// Loop over parameters and write parameter files for CAMB:
process_number=0;
jobs=0;
Pjobs=0;
octo_counter=0;

if (mode==3 && remove!=0) // remove old Blue Gene job description files before writing new ones.
{
	printf("Removing old Blue Gene job description files.\n");
	sprintf(command, "rm %s/%s/%s/%s/*.ll", ics_dir_speccomp, series_folder, NGenIC_folder, NGenIC_jobs_folder);
	system(command);
	sprintf(command, "rm %s/%s/%s/%s/*.ll", ics_dir_speccomp, series_folder, Gadget_folder, Gadget_jobs_folder);
	system(command);
}

if (mode==4)
{
  sprintf(directory_script_filename, "%s/%s/directory_creation_%s.sh", Gadget_output_dir, series_folder, series_folder);
  printf("Opening file (directory generation script) for writing:\n%s\n", directory_script_filename);
  fflush(stdout);

  directory_script_file=fopen(directory_script_filename, "w");
  fprintf(directory_script_file, "#!/bin/bash\n\n");

}

if (mode==5)
{
  sprintf(IG_directory_script_filename, "%s/%s/IG_directory_creation_%s.sh", IG_output_dir, series_folder, series_folder);
  printf("Opening file (IG directory generation script) for writing:\n%s\n", IG_directory_script_filename);
  fflush(stdout);

  IG_directory_script_file=fopen(IG_directory_script_filename, "w");
  fprintf(IG_directory_script_file, "#!/bin/bash\n\n");

}

if (mode==6)
{
  sprintf(simple_list_filename, "./simple_simulation_code_list.txt");
  printf("Opening file (IG directory generation script) for writing:\n%s\n", simple_list_filename);
  fflush(stdout);

  simple_list_file=fopen(simple_list_filename, "w");

}



if (mode==1) // open condor job description file for writing and write header
{
	sprintf(job_description_filename, "%s/%s/%s/%s/CAMB_job_desc", ics_dir_speccomp, series_folder, CAMB_folder, CAMB_jobs_folder);
	CAMB_condor_file=fopen(job_description_filename, "w");
	
	fprintf(CAMB_condor_file,  "# Condor Job Description File for CAMB for generation of power spectra for N-body simulations\n"); 
	fprintf(CAMB_condor_file, "Executable = %s/%s\n", CAMB_dir, CAMB_exec); 
	fprintf(CAMB_condor_file, "Universe = vanilla\n");
	fprintf(CAMB_condor_file, "Requirements = (CPU_Experiment == \"lsst\")\n");
	// fprintf(CAMB_condor_file, "requirements = Memory >= 1200\n");
	fprintf(CAMB_condor_file, "notification = Complete\n");
	fprintf(CAMB_condor_file, "notify_user = jank@astro.columbia.edu\n");
	// fprintf(CAMB_condor_file, "image_size = 1200000\n"
	fprintf(CAMB_condor_file, "initialdir = %s/%s/%s/%s\n", ics_dir_speccomp, series_folder, CAMB_folder, CAMB_data_folder);
	fprintf(CAMB_condor_file, "\n"); 
	fprintf(CAMB_condor_file, "should_transfer_files = YES\n");
	fprintf(CAMB_condor_file, "stream_output = TRUE\n");
	fprintf(CAMB_condor_file, "when_to_transfer_output = ON_EXIT\n\n");
}

for (i_OBh2=0; i_OBh2<Nobh2; i_OBh2++)
{
for (i_OM=0; i_OM<Nom; i_OM++)
{
for (i_OL=0; i_OL<Nol; i_OL++)
{
for (i_w0=0; i_w0<Nw0; i_w0++)
{
for (i_wa=0; i_wa<Nwa; i_wa++)
{
for (i_ns=0; i_ns<Nns; i_ns++)
{
for (i_As=0; i_As<Nas; i_As++)
{
for (i_h=0; i_h<Nh; i_h++)
{
for (i_z=0; i_z<Nz; i_z++)
{
	// Calculate derived quantities:
	if (flat_universe!=0) OL[i_OL]=1-OM[i_OM]; // override to make flat universe if flat_universe flag is set;
	OK=1-OM[i_OM]-OL[i_OL]; // curvature
	OCh2=OM[i_OM]*h[i_h]*h[i_h]-OBh2[i_OBh2]; // Fractional CDM density * h^2
	Dplus=0; 
	vel_prefac_lam=0.0; // initialization to zero.
	sprintf(filerawbase, "Om%1.3f_Ol%1.3f_w%1.3f_ns%1.3f", OM[i_OM], OL[i_OL], w0[i_w0], ns[i_ns]); // base of filenames for a particular parameter combination (must contain distinction based on parameter varied values).
	sprintf(filebase, "%s-%s", series, filerawbase);
	sprintf(CAMB_param_filename, "params_%s.ini", filebase); // construct name of CAMB parameter file (input file for CAMB).
	sprintf(power_spectrum_filename, "%s_%s", filebase, power_end); // construct filename of file which will contain output power spectrum of CAMB (output of CAMB).
	sprintf(converted_power_spectrum_filename, "%s_%s.txt", power_front, filebase); // This file contains the power spectrum in final format and will be input file for N-GenIC. MPS stands for (total) Matter Power Spectrum 

	//printf("Looping: %d %d %d\n", i_OM, i_OL, i_w0);
	//fflush(stdout);

	if (mode==1) // write CAMB parameter file
	{
		if (power_spectrum_at_zini==0)
		{
			printf("NOTE: Setting up power spectrum to be generated at z=0 and scaled back later by N-GenIC.\n");
			redshift=0;
		}
		else redshift=z[i_z];
		write_CAMB_parameter_file(CAMB_param_filename, filebase, power_end, OBh2[i_OBh2], OCh2, OM[i_OM], OL[i_OL], OK, w0[i_w0], wa[i_wa], ns[i_ns], As[i_As], h[i_h], redshift);
		write_CAMB_condor_job_description(CAMB_condor_file, CAMB_param_filename, filebase);
		process_number++;
	}
	else if (mode==2)
	{
	  /*
		// Simulate power spectrum file: (turn this off for actual production run, is here only for code test)
		FILE *sim_spec_file;
		char sim_spec_filename[1000];
		int ii;
		sprintf(sim_spec_filename, "%s/%s/%s/%s/%s", ics_dir_speccomp, series_folder, CAMB_folder, CAMB_data_folder, power_spectrum_filename);
		sim_spec_file=fopen(sim_spec_filename, "w");
		for (ii=0; ii<50; ii++) fprintf(sim_spec_file, "%e %e\n", ((double) ii+10), ((double) ii*ii+100));
		fclose(sim_spec_file);
		// End of power spectrum file simulation.
	   */
		// Convert power spectrum file:
		convert_CAMB_power_spectrum(power_spectrum_filename, converted_power_spectrum_filename);
	}
	else if (mode==3 || mode==4 || mode==5 || mode==6) // generate all N-GenIC and Gadget-2 parameter files, and selected Blue Gene job description files and submission shell scripts. 
	{
		
		if (mode==3)
		{
		  
			// Calculate Growth Factor and Velocity prefactor (depends on cosmological parameters other than sigma_8, does not depend on boxsize or random seed):
		        Dplus=0.0;
	         	vel_prefac_lam=0.0;
	     		Dplus=Dplus_Interface(z[i_z], 0.0, OM[i_OM], OL[i_OL], 0.0, w0[i_w0], wa[i_wa], h[i_h]); // redshift today always set to z_today=0.0, same for Omega Cosmological constant (additional cosmological constant in addition to dark energy model).
			// vel_prefac_lam (realized as a global variable for simplicity) is also set by this function.
	
			// Make a safety check, because the above is C calling FORTRAN 77:
			printf("Growth factor and velocity prefactor from Lam Hui's growth factor code: Dplus=%f, vel_prefac_lam=%f\n", Dplus, vel_prefac_lam);
			if (fabs(Dplus)<0.001 || fabs(vel_prefac_lam)<0.001)
			{
				printf("ERROR: Growth factor or velocity prefactor returned incorrectly from Lam Hui's growth factor code. Check for compilation errors (this is C calling FORTRAN 77 for this part) and or warnings/limitations of parameter ranges in growth factor routine package.\n");
				exit(2);
			}
		  
			// exit(1111);
            
		}

		for (i_s8=0; i_s8<Ns8; i_s8++)
		{
		for (i_seed=0; i_seed<Nseed; i_seed++)
		{
		for (i_boxsize=0; i_boxsize<Nboxsize; i_boxsize++)
		{
		  // Caluclate gravitational softening length based on boxsize and number of particles:
		  soft=7.5*512.0*boxsize[i_boxsize]/(200.0*part); // This is scaled formula from standard simulation (similar softening to Millennium simulation, Bologna simulation, etc.).
		  
			//printf("Seed: seed, iseed: %d %d\n", seed[i_seed], i_seed);
			sprintf(simulation_codename, "%s-%db%d_%s_si%1.3f_ic%d", series, part, ((int) boxsize[i_boxsize]), filerawbase, s8[i_s8], i_seed+1);
			sprintf(filebase2, "%s_%s", ics_front, simulation_codename);
			sprintf(NGenIC_param_filename, "%s.param", filebase2); // This file contains the power spectrum in final format and will be input file for N-GenIC. MPS stands for (total) Matter Power Spectrum 
			sprintf(Gadget_param_filename, "%s.param", simulation_codename);
			
			// write N-GenIC and Gadget-2 parameter files for all possible runs.
			if (mode==3)
			{
			        write_NGenIC_parameter_file(converted_power_spectrum_filename, NGenIC_param_filename, filebase2, part, boxsize[i_boxsize], OBh2[i_OBh2], OCh2, OM[i_OM], OL[i_OL], OK, w0[i_w0], wa[i_wa], 1.0, s8[i_s8], h[i_h], z[i_z], seed[i_seed], power_spectrum_at_zini, Dplus, vel_prefac_lam); // ns must always be 1.0, because this is an _additional_ tilt N-GenIC applies. Variable A_s (primordial amplitude) now replaced by s8 (sigma_8).
				write_Gadget_parameter_file(Gadget_param_filename, simulation_codename, filebase2, output_list_filename, boxsize[i_boxsize], OBh2[i_OBh2], OM[i_OM], OL[i_OL], w0[i_w0], wa[i_wa], h[i_h], z[i_h], soft); // associated parameter file for Gadget-2 N-body run with these initial conditions (contains run parameter optimization).
			}
			
			// Create job description files for Blue Gene/L if this job is to be run for the selected submission_style:
			if (mode==4)
			{
				if (submission_style==1 || (submission_style==2 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)<=1) || (submission_style==3 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)==0) || (submission_style==4 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)>1) || (submission_style==5 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)==1) || (submission_style==6 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)!=0))
				{
					jobs++;
					write_BGL_description(NGenIC_param_filename, Gadget_param_filename, jobstammL, jobs);
					octo_counter++;
					sprintf(BGP_NGenIC_param_filename[octo_counter-1], "%s", NGenIC_param_filename);
					sprintf(BGP_Gadget_param_filename[octo_counter-1], "%s", Gadget_param_filename);
					if (octo_counter==8)
					{
						Pjobs++;
						write_BGP_description(BGP_NGenIC_param_filename, BGP_Gadget_param_filename, jobstammP, Pjobs, octo_counter);
						octo_counter=0;
					}

					// Add directory to directory creation script:
					write_directory_script(directory_script_file, simulation_codename);

				}
			}

			if (mode==5)
			  {

			    if (submission_style==1 || (submission_style==2 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)<=1) || (submission_style==3 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)==0) || (submission_style==4 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)>1) || (submission_style==5 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)==1) || (submission_style==6 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)!=0))
			      {
				// Add directory to IG directory creation script:
				write_IG_directory_script(IG_directory_script_file, simulation_codename);
			      }
			  }

				
		} // end loop over i_boxsize.
		} // end loop over i_seed.

		if (mode==6)
		  {
		    // Add cosmology to simple list file:
		    if (submission_style==1 || (submission_style==2 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)<=1) || (submission_style==3 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)==0) || (submission_style==4 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)>1) || (submission_style==5 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)==1) || (submission_style==6 && (i_OBh2!=0)+(i_OM!=0)+(i_OL!=0)+(i_w0!=0)+(i_wa!=0)+(i_ns!=0)+(i_s8!=0)!=0))
		      {
			write_simple_list(simple_list_file, simulation_codename);
		      }
		  }

		} // end loop over i_s8.
	}
	else
	{
		printf("ERROR: Selected mode does not exist. Select a valid mode instead.\n");
		exit(1);
	}
	
	
	// printf("Done with loop.\n");
	// fflush(stdout);

} // end loop over i_z.
} // end loop over i_h.
} // end loop over i_As.
} // end loop over i_ns.
} // end loop over i_wa.
} // end loop over i_w0.
} // end loop over i_OL.
} // end loop over i_OM.
} // end loop over i_OBh2.
	
	
	if (mode == 1) fclose(CAMB_condor_file);
	if (mode == 4)
	{
		// write shell script which submits all NYBlue/L jobs upon execution:
		sprintf(BG_type, "L");
		write_submission_script_L(jobs, jobstammL, BG_type);
		
		if (octo_counter!=0) // write rest of Octopus jobs which don't add up to final eight.
		{
			Pjobs++;
			write_BGP_description(BGP_NGenIC_param_filename, BGP_Gadget_param_filename, jobstammP, Pjobs, octo_counter);
		}
		sprintf(BG_type, "P");
		write_submission_script_P(jobs, Pjobs, jobstammP, BG_type);

		fclose(directory_script_file);

		//set appropriate execution permissions 
		chmod(directory_script_filename,S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	}
        if (mode == 5) fclose(IG_directory_script_file);
        if (mode == 6) fclose(simple_list_file);
	
    printf("Done!!\n\n");
    return 0;
}


// Convert power spectrum from CAMB output format to N-GenIC input format: {k, P} --> {ln(k)/ln(10), ln(k^3*P)/(2*pi^2)/ln(10)}, where k is in (Mpc/h)^(-1).
void convert_CAMB_power_spectrum(char power_spectrum_filename[], char converted_power_spectrum_filename[])
{
	FILE *power_spectrum_file, *converted_power_spectrum_file;
	char filename_in[1000], filename_out[1000];
	char line[1000];
	int line_length;
	float k, p;
	int line_counter, stopwriting;
	line_counter=0; stopwriting=0;
	double logk, logp, prevk, prevp, pi;
	pi=3.14159265359;
	k=0; p=0; logk=0; logp=0; prevk=0; prevp=0;
	
	sprintf(filename_in, "%s/%s/%s/%s/%s", ics_dir_speccomp, series_folder, CAMB_folder, CAMB_data_folder, power_spectrum_filename);
	power_spectrum_file=fopen(filename_in, "r");
	sprintf(filename_out, "%s/%s/%s/%s/%s", ics_dir_speccomp, series_folder, NGenIC_folder, NGenIC_data_folder, converted_power_spectrum_filename);
	converted_power_spectrum_file=fopen(filename_out, "w");
	// printf("Converting CAMB power spectrum:\n%s\ninto N-GenIC power spectrum:\n%s\n", power_spectrum_filename, converted_power_spectrum_filename);
	printf("Converting CAMB power spectrum:\n%s\ninto N-GenIC power spectrum:\n%s\n", filename_in, filename_out);
	fflush(stdout);

	// Read each line, convert and write. Need fgetline-while structure.
	
	while ((line_length=fgetline(power_spectrum_file, line, sizeof(line))) > 0)
	{
	  // printf("Reading a line from input file:\nLine length: %d \nLine Content: %s\n", line_length, line);
		if (sscanf(line, "%e %e", &k, &p) == 2)
		{
			line_counter++;
			
			// printf("Test 1: %e %e\n", k, p);
			// fflush(stdout);

			logk=log(k)/log(10.0); // log() in C is natural logarithm ln().
			logp=log(k*k*k*p/(2.0*pi*pi))/log(10.0);

			// printf("Test 2: %e %e\n", logk, logp);
			// fflush(stdout);
			
			if (line_counter>1)
			{
				if(fabs(logp-prevp)<0.1 && stopwriting==0) // catches outliers in CAMB file (happen towards the right (large-k) end of power spectrum around k~1000, which is ca. as high as we need to go for the N-body simulations.
				{
			    	      	fprintf(converted_power_spectrum_file, "%e %e\n", prevk, prevp);
					prevk=logk;
					prevp=logp;
				}
				else
				{
					printf("WARNING: Rejected power spectrum data point: log(k)=%e, log(Delta)=%e; previous point was log(k)=%e, log(Delta)=%e\n", logk, logp, prevk, prevp);
					stopwriting=1;
					if (logk<2.9)
					  {
					    printf("ERROR: Accuracy problem in matter power spectrum from CAMB, first outlier to the right occurs at smaller k than selected threshold log(k)=2.9.\nCheck files and simulations if really such large k is needed, and rerun.\n");
					    exit(1);
					  }
				}
			}
			else
			{
				prevk=logk;
				prevp=logp;
			}
		}
	}
	
	fclose(converted_power_spectrum_file);
	fclose(power_spectrum_file);

	printf("Done with a power spectrum conversion.\n");
	fflush(stdout);
}


void write_CAMB_parameter_file(char CAMB_param_filename[], char filebase[], char power_end[], double OBh2, double OCh2, double OM, double OL, double OK, double w0, double wa, double ns, double As, double h, double z)
{
	FILE* param_file;
	char filename[1000];
	sprintf(filename, "%s/%s/%s/%s/%s", ics_dir_speccomp, series_folder, CAMB_folder, CAMB_param_folder, CAMB_param_filename);
	// Write parameter file for CAMB:
	param_file=fopen(filename, "w");

	printf("Writing CAMB parameter file:\n%s\n", filename);
	fflush(stdout);

	/////////////////////////////////////////////
	/////////////////////////////////////////////
	
	fprintf(param_file, "#Parameters for CAMB\n\n");

	fprintf(param_file, "#output_root is prefixed to output file names\n");
	fprintf(param_file, "output_root = %s\n", filebase); // no path needed here, Condor just has to have the proper working directory (CAMB data output folder).

	fprintf(param_file,"#What to do\n");
	fprintf(param_file, "get_scalar_cls = T\n");
	fprintf(param_file, "get_vector_cls = F\n");
	fprintf(param_file, "get_tensor_cls = F\n");
	fprintf(param_file, "get_transfer   = T\n\n");

	fprintf(param_file, "#if do_lensing then scalar_output_file contains additional columns of l^4 C_l^{pp} and l^3 C_l^{pT}\n");
	fprintf(param_file, "#where p is the projected potential. Output lensed CMB Cls (without tensors) are in lensed_output_file below.\n");
	fprintf(param_file, "do_lensing     = F\n\n");

	fprintf(param_file, "# 0: linear, 1: non-linear matter power (HALOFIT), 2: non-linear CMB lensing (HALOFIT)\n");
	fprintf(param_file, "do_nonlinear = 0\n\n");

	fprintf(param_file, "#Maximum multipole and k*eta.\n"); 
	fprintf(param_file, "#  Note that C_ls near l_max are inaccurate (about 5%%), go to 50 more than you need\n");
	fprintf(param_file, "#  Lensed power spectra are computed to l_max_scalar-250 where accurate at %%-level\n");
	fprintf(param_file, "#  For high accuracy lensed spectra set l_max_scalar = (l you need) + 500\n");
	fprintf(param_file, "#  To get accurate lensed BB need to have l_max_scalar>2000, k_eta_max_scalar > 10000\n");
	fprintf(param_file, "#  Otherwise k_eta_max_scalar=2*l_max_scalar usually suffices\n");
	fprintf(param_file, "l_max_scalar      = 8000\n"); // 8000 for final run
	fprintf(param_file, "k_eta_max_scalar  = 16000\n"); // 16000
	fprintf(param_file, "#l_max_scalar      = 2000\n");
	fprintf(param_file, "#k_eta_max_scalar  = 4000\n\n");

	fprintf(param_file, "#  Tensor settings should be less than or equal to the above\n");
	fprintf(param_file, "l_max_tensor      = 1500\n");
	fprintf(param_file, "k_eta_max_tensor  = 3000\n\n");

	fprintf(param_file, "#Main cosmological parameters, neutrino masses are assumed degenerate\n");
	fprintf(param_file, "# If use_phyical set phyiscal densities in baryone, CDM and neutrinos + Omega_k\n");
	fprintf(param_file, "use_physical   = T\n");
	fprintf(param_file, "ombh2          = %f\n", OBh2); //  0.0227
	fprintf(param_file, "omch2          = %f\n", OCh2); //  0.112084
	fprintf(param_file, "omnuh2         = 0\n"); // no neutrinos implemented so far.
	fprintf(param_file, "omk            = %f\n", OK); //  0
	fprintf(param_file, "hubble         = %f\n", h*100.0); //  72
	fprintf(param_file, "#effective equation of state parameter for dark energy, assumed constant\n");
	fprintf(param_file, "w              = %f\n", w0); //  -1
	fprintf(param_file, "#constant comoving sound speed of the dark energy (1=quintessence)\n");
	fprintf(param_file, "cs2_lam        = 1\n\n");

	fprintf(param_file, "#if use_physical = F set parameters as here\n");
	fprintf(param_file, "#omega_baryon   = 0.0462\n");
	fprintf(param_file, "#omega_cdm      = 0.2538\n");
	fprintf(param_file, "#omega_lambda   = 0.7\n");
	fprintf(param_file, "#omega_neutrino = 0\n\n");

	fprintf(param_file, "#massless_neutrinos is the effective number (for QED + non-instantaneous decoupling)\n");
	fprintf(param_file, "temp_cmb           = 2.726\n");
	fprintf(param_file, "helium_fraction    = 0.24\n");
	fprintf(param_file, "massless_neutrinos = 3.04\n");
	fprintf(param_file, "massive_neutrinos  = 0\n\n");

	fprintf(param_file, "#Neutrino mass splittings\n");
	fprintf(param_file, "nu_mass_eigenstates = 1\n");
	fprintf(param_file, "#nu_mass_degeneracies = 0 sets nu_mass_degeneracies = massive_neutrinos\n");
	fprintf(param_file, "#otherwise should be an array\n");
	fprintf(param_file, "#e.g. for 3 neutrinos with 2 non-degenerate eigenstates, nu_mass_degeneracies = 2 1\n");
	fprintf(param_file, "nu_mass_degeneracies = 0\n");  
	fprintf(param_file, "#Fraction of total omega_nu h^2 accounted for by each eigenstate, eg. 0.5 0.5\n");
	fprintf(param_file, "nu_mass_fractions = 1\n\n");

	fprintf(param_file, "#Initial power spectrum, amplitude, spectral index and running. Pivot k in Mpc^{-1}.\n");
	fprintf(param_file, "initial_power_num         = 1\n");
	fprintf(param_file, "pivot_scalar              = 0.002\n");
	fprintf(param_file, "# Alternative popular value for pivot_scalar: 0.05\n");
	fprintf(param_file, "pivot_tensor              = 0.002\n");
	fprintf(param_file, "scalar_amp(1)             = 2.41e-9\n");
	fprintf(param_file, "scalar_spectral_index(1)  = %f\n", ns); // 0.96
	fprintf(param_file, "scalar_nrun(1)            = 0\n");
	fprintf(param_file, "tensor_spectral_index(1)  = 0\n");
	fprintf(param_file, "#ratio is that of the initial tens/scal power spectrum amplitudes\n");
	fprintf(param_file, "initial_ratio(1)          = 0\n");
	fprintf(param_file, "#initial_ratio(1)          = 1\n");
	fprintf(param_file, "#note vector modes use the scalar settings above\n\n\n");


	fprintf(param_file, "#Reionization, ignored unless reionization = T, re_redshift measures where x_e=0.5\n");
	fprintf(param_file, "reionization         = T\n\n");

	fprintf(param_file, "re_use_optical_depth = T\n");
	fprintf(param_file, "re_optical_depth     = 0.087\n"); // currently not implemented yet as variable parameter.
	fprintf(param_file, "#If re_use_optical_depth = F then use following, otherwise ignored\n");
	fprintf(param_file, "re_redshift          = 11\n");
	fprintf(param_file, "#width of reionization transition. CMBFAST model was similar to re_delta_redshift~0.5.\n");
	fprintf(param_file, "re_delta_redshift    = 0.5\n");
	fprintf(param_file, "#re_delta_redshift    = 1.5\n");
	fprintf(param_file, "#re_ionization_frac=-1 sets to become fully ionized using YHe to get helium contribution\n");
	fprintf(param_file, "#Otherwise x_e varies from 0 to re_ionization_frac\n");
	fprintf(param_file, "re_ionization_frac   = -1\n\n\n");


	fprintf(param_file, "#RECFAST 1.4 recombination parameters\n");
	fprintf(param_file, "RECFAST_fudge = 1.14\n");
	fprintf(param_file, "RECFAST_fudge_He = 0.86\n");
	fprintf(param_file, "RECFAST_Heswitch = 6\n\n\n");


	fprintf(param_file, "#Initial scalar perturbation mode (adiabatic=1, CDM iso=2, Baryon iso=3, \n");
	fprintf(param_file, "# neutrino density iso =4, neutrino velocity iso = 5) \n");
	fprintf(param_file, "initial_condition   = 1\n");
	fprintf(param_file, "#If above is zero, use modes in the following (totally correlated) proportions\n");
	fprintf(param_file, "#Note: we assume all modes have the same initial power spectrum\n");
	fprintf(param_file, "initial_vector = -1 0 0 0 0\n\n");

	fprintf(param_file, "#For vector modes: 0 for regular (neutrino vorticity mode), 1 for magnetic\n");
	fprintf(param_file, "vector_mode = 0\n\n");

	fprintf(param_file, "#Normalization\n");
	fprintf(param_file, "COBE_normalize = F\n");
	fprintf(param_file, "##CMB_outputscale scales the output Cls\n");
	fprintf(param_file, "#To get MuK^2 set realistic initial amplitude (e.g. scalar_amp(1) = 2.3e-9 above) and\n");
	fprintf(param_file, "#otherwise for dimensionless transfer functions set scalar_amp(1)=1 and use\n");
	fprintf(param_file, "#CMB_outputscale = 1\n");
	fprintf(param_file, "CMB_outputscale = 7.4311e12\n\n");

	fprintf(param_file, "#Transfer function settings, transfer_kmax=0.5 is enough for sigma_8\n");
	fprintf(param_file, "#transfer_k_per_logint=0 sets sensible non-even sampling; \n");
	fprintf(param_file, "#transfer_k_per_logint=5 samples fixed spacing in log-k\n");
	fprintf(param_file, "#transfer_interp_matterpower =T produces matter power in regular interpolated grid in log k; \n");
	fprintf(param_file, "# use transfer_interp_matterpower =F to output calculated values (e.g. for later interpolation)\n");
	fprintf(param_file, "transfer_high_precision = T\n");
	fprintf(param_file, "transfer_kmax           = 1000\n"); // 1000 for final run
	fprintf(param_file, "transfer_k_per_logint   = 100\n");
	fprintf(param_file, "transfer_num_redshifts  = 1\n");
	fprintf(param_file, "transfer_interp_matterpower = T\n");
	fprintf(param_file, "transfer_redshift(1)    = %f\n", z);
	fprintf(param_file, "transfer_filename(1)    = CAMB-TFR.dat\n");  //   transfer_out_z0.da
	fprintf(param_file, "#Matter power spectrum output against k/h in units of h^{-3} Mpc^3\n");
	fprintf(param_file, "transfer_matterpower(1) = %s\n", power_end);  //   matterpower_z0.dat
	fprintf(param_file, "transfer_power_var = 7\n\n"); // 7 sets the matter power spectrum to be for total matter.

	fprintf(param_file, "#Output files not produced if blank. make camb_fits to use use the FITS setting.\n");
	fprintf(param_file, "scalar_output_file = scalCls.dat\n");
	fprintf(param_file, "vector_output_file = vecCls.dat\n");
	fprintf(param_file, "tensor_output_file = tensCls.dat\n");
	fprintf(param_file, "total_output_file  = totCls.dat\n");
	fprintf(param_file, "lensed_output_file = lensedCls.dat\n");
	fprintf(param_file, "lensed_total_output_file  =lensedtotCls.dat\n");
	fprintf(param_file, "FITS_filename      = scalCls.fits\n\n");

	fprintf(param_file, "##Optional parameters to control the computation speed,accuracy and feedback\n\n");

	fprintf(param_file, "#If feedback_level > 0 print out useful information computed about the model\n");
	fprintf(param_file, "feedback_level = 1\n\n");

	fprintf(param_file, "# 1: curved correlation function, 2: flat correlation function, 3: inaccurate harmonic method\n");
	fprintf(param_file, "lensing_method = 1\n");
	fprintf(param_file, "accurate_BB = T\n\n\n");


	fprintf(param_file, "#massive_nu_approx: 0 - integrate distribution function\n");
	fprintf(param_file, "#                   1 - switch to series in velocity weight once non-relativistic\n");
	fprintf(param_file, "#                   2 - use fast approximate scheme (CMB only- accurate for light neutrinos)\n");
	fprintf(param_file, "#                   3 - intelligently use the best accurate method\n");
	fprintf(param_file, "massive_nu_approx = 3\n\n");

	fprintf(param_file, "#Whether you are bothered about polarization. \n");
	fprintf(param_file, "accurate_polarization   = T\n\n");

	fprintf(param_file, "#Whether you are bothered about percent accuracy on EE from reionization\n");
	fprintf(param_file, "accurate_reionization   = T\n\n");

	fprintf(param_file, "#whether or not to include neutrinos in the tensor evolution equations\n");
	fprintf(param_file, "do_tensor_neutrinos     = F\n\n");

	fprintf(param_file, "#Whether to turn off small-scale late time radiation hierarchies (save time,v. accurate)\n");
	fprintf(param_file, "do_late_rad_truncation   = F\n\n");

	fprintf(param_file, "#Computation parameters\n");
	fprintf(param_file, "#if number_of_threads=0 assigned automatically\n");
	fprintf(param_file, "number_of_threads       = 0\n\n");

	fprintf(param_file, "#Default scalar accuracy is about 0.3%% (except lensed BB). \n");
	fprintf(param_file, "#For 0.1%%-level try accuracy_boost=2, l_accuracy_boost=2.\n\n");

	fprintf(param_file, "#Increase accuracy_boost to decrease time steps, use more k values,  etc.\n");
	fprintf(param_file, "#Decrease to speed up at cost of worse accuracy. Suggest 0.8 to 3.\n");
	fprintf(param_file, "accuracy_boost          = 3\n\n"); // 3 for final run.

	fprintf(param_file, "#Larger to keep more terms in the hierarchy evolution. \n");
	fprintf(param_file, "l_accuracy_boost        = 3\n\n"); // 3

	fprintf(param_file, "#Increase to use more C_l values for interpolation.\n");
	fprintf(param_file, "#Increasing a bit will improve the polarization accuracy at l up to 200 -\n");
	fprintf(param_file, "#interpolation errors may be up to 3%%\n");
	fprintf(param_file, "#Decrease to speed up non-flat models a bit\n");
	fprintf(param_file, "l_sample_boost          = 3\n"); // 3

	
	/////////////////////////////////////////////
	/////////////////////////////////////////////
	

	fclose(param_file);
}


void write_CAMB_condor_job_description(FILE* CAMB_condor_file, char CAMB_param_filename[], char filebase[])
{
	 fprintf(CAMB_condor_file, "transfer_input_files = %s/%s/%s/%s/%s\n", ics_dir_speccomp, series_folder, CAMB_folder, CAMB_param_folder, CAMB_param_filename);
     fprintf(CAMB_condor_file, "Error           = %s/%s/%s/%s/err_CAMB.%s\n", ics_dir_speccomp, series_folder, CAMB_folder, CAMB_logs_folder, filebase);
     fprintf(CAMB_condor_file, "Output          = %s/%s/%s/%s/out_CAMB.%s\n", ics_dir_speccomp, series_folder, CAMB_folder, CAMB_logs_folder, filebase);
     fprintf(CAMB_condor_file, "Log             = %s/%s/%s/%s/log_CAMB.%s\n", ics_dir_speccomp, series_folder, CAMB_folder, CAMB_logs_folder, filebase);
     fprintf(CAMB_condor_file, "Arguments = %s\n", CAMB_param_filename); // local path is enough here, because it gets copied by Condor first to local job execution folder. 
     fprintf(CAMB_condor_file, "Queue\n\n\n");
}



void write_NGenIC_parameter_file(char converted_power_spectrum_filename[], char NGenIC_param_filename[], char filebase2[], int part, double boxsize, double OBh2, double OCh2, double OM, double OL, double OK, double w0, double wa, double ns, double s8, double h, double z, int seed, int power_spectrum_at_zini, double Dplus, double vel_prefac_lam)
{
	FILE* param_file;
	char filename[1000];
	int renormalize_spectrum;
	
	if (power_spectrum_at_zini==0) renormalize_spectrum=1;
	else renormalize_spectrum=0;
	
	sprintf(filename, "%s/%s/%s/%s/%s", ics_dir_speccomp, series_folder, NGenIC_folder, NGenIC_param_folder, NGenIC_param_filename);
	// Write parameter file for N-GenIC:
	param_file=fopen(filename, "w");
	
	printf("Writing N-GenIC parameter file:\n%s\n", filename);
	fflush(stdout);
	
	printf("WARNING: ns in N-GenIC parameter file is an additional spectral index for additional tilt. Leave at ns=1 (unless want to tilt a previously untilted power spectrum.\n");
	
	/////////////////////////////////
	/////////////////////////////////
	
	fprintf(param_file, "%% Automatically generated parameter file for w(z)-enhanced version of public N-GenIC. Do not modify cosmological parameters by hand with the exception of sigma_8!\n");
	fprintf(param_file, "%% Use standard parameter file and regular version of N-GenIC if want to do hand-modifications. Limitation then: no w(z) dark energy models, only cosmological constant.\n");
	
	fprintf(param_file, "Nmesh            %d        %% This is the size of the FFT grid used to \n", 2*part); // 1024
	

	
	fprintf(param_file, "                           %% compute the displacement field. One\n");
	fprintf(param_file, "                           %% should have Nmesh >= Nsample.\n\n");

	fprintf(param_file, "Nsample          %d        %% sets the maximum k that the code uses,\n", part); // 512
	fprintf(param_file, "                           %% i.e. this effectively determines the\n");
	fprintf(param_file, "                           %% Nyquist frequency that the code assumes,\n");
	fprintf(param_file, "                           %% k_Nyquist = 2*PI/Box * Nsample/2\n");
	fprintf(param_file, "                           %% Normally, one chooses Nsample such that\n");
	fprintf(param_file, "                           %% Ntot =  Nsample^3, where Ntot is the\n");
	fprintf(param_file, "                           %% total number of particles\n\n\n");
 


	fprintf(param_file, "Box              %f   %% Periodic box size of simulation\n\n", boxsize*1000.0); // boxsize here must be in kpc/h. // 200000.0

	fprintf(param_file, "FileBase         %s                 %% Base-filename of output files\n", filebase2); // ics_i512w10_3
	fprintf(param_file, "OutputDir        %s/%s/%s/%s                 %%/home/volker/ics/   %% Directory for output\n\n", ics_data_dir_simcomp, series_folder, Gadget_folder, Gadget_data_folder); // /gpfs/scratch3/jank/Storage/sims/ics/i-series/

	fprintf(param_file, "GlassFile        %s/dummy_glass_big_endian.dat  %% Grid-File or Glass-File \n\n", NGenIC_dir_L); // select between "grid_little_endian.dat" and "grid_bid_endian.dat", depending on endianness of machine on which N-GenIC and Gadget-2 will be run. 

	fprintf(param_file, "TileFac     %d                %% Number of times the glass file is\n", part/16); // above grid files consist of 16^3 particles, so they need to be repoduced (number of particles in one dim / 16)-times. If want number of particles not divisible by 16, must create new grid files.
	fprintf(param_file, "                                  %% tiled in each dimension (must be\n");
	fprintf(param_file, "                                  %% an integer)\n\n");

	fprintf(param_file, "                                  %% NOTE: The dark matter particle\n");
	fprintf(param_file, "                                  %% number \"Ntot\" will be equal to the number\n");
	fprintf(param_file, "                                  %% of particles in the file dummy_glass.dat\n");
	fprintf(param_file, "                                  %% multiplied with GlassTileFac^3.\n\n\n\n");



	fprintf(param_file, "Omega            %f       %% Total matter density  (at z=0)\n", OM);
	fprintf(param_file, "OmegaLambda      %f       %% Cosmological constant (at z=0)\n", OL);
	fprintf(param_file, "OmegaBaryon      %f       %% Baryon density        (at z=0)\n", OBh2/(h*h));
	fprintf(param_file, "HubbleParam      %f       %% Hubble paramater (may be used for power spec parameterization)\n\n", h);

	fprintf(param_file, "Redshift         %f        %% Starting redshift\n", z);
	fprintf(param_file, "w0               %f        %% Dark energy model w(z)=w0+(z/(1.0+z))*wa\n", w0);
	fprintf(param_file, "wa               %f        %% Dark energy model parameter\n", wa);
	fprintf(param_file, "GrowthFactor     %f        %% Growth Factor (Dplus). Depends on Starting Redshift and Cosmological Parameters!\n", Dplus);
	fprintf(param_file, "                           %% Encodes dark energy parameters for w(z) model which are not passed to N-GenIC otherwise.\n");
	fprintf(param_file, "                           %% Needs to be calculated specifically for enhanced version of N-GenIC, do not modify by hand!\n");
	fprintf(param_file, "VelocityPrefactor  %f      %% Velocity Prefactor (vel_prefac). Depends on Starting Redshift and Cosmological Parameters!\n", vel_prefac_lam);
	fprintf(param_file, "                           %% Needs to be calculated specifically for enhanced version of N-GenIC, do not modify by hand!\n");
	fprintf(param_file, "Sigma8           %f        %% power spectrum normalization\n\n\n", s8);



	fprintf(param_file, "SphereMode       1         %% if \"1\" only modes with |k| < k_Nyquist are\n");
	fprintf(param_file, "                           %% used (i.e. a sphere in k-space), otherwise modes with\n");
	fprintf(param_file, "                           %% |k_x|,|k_y|,|k_z| < k_Nyquist are used\n");
	fprintf(param_file, "                           %% (i.e. a cube in k-space)\n\n\n");
          

	fprintf(param_file, "WhichSpectrum    2         %% \"1\" selects Eisenstein & Hu spectrum,\n");
	fprintf(param_file, "                           %% \"2\" selects a tabulated power spectrum in\n");
	fprintf(param_file, "                           %% the file 'FileWithInputSpectrum'\n");
	fprintf(param_file, "                           %% otherwise, Efstathiou parametrization is used\n\n\n");




	fprintf(param_file, "FileWithInputSpectrum   %s/%s/%s/%s/%s\n                            %% filename of tabulated input\n", ics_dir_simcomp, series_folder, NGenIC_folder, NGenIC_data_folder, converted_power_spectrum_filename);
	fprintf(param_file, "						                           %% spectrum (if used)\n");
	fprintf(param_file, "InputSpectrum_UnitLength_in_cm  3.085678e24 %% defines length unit of tabulated\n");
	fprintf(param_file, "                                            %% input spectrum in cm/h. \n");
	fprintf(param_file, "                                            %% Note: This can be chosen different from UnitLength_in_cm\n\n");

	fprintf(param_file, "ReNormalizeInputSpectrum   %d                %% if set to zero, the\n", renormalize_spectrum);  // 0 (but 1 now newly).
	fprintf(param_file, "                                            %% tabulated spectrum is\n");
	fprintf(param_file, "                                            %% assumed to be normalized\n");
	fprintf(param_file, "                                            %% already in its amplitude to\n");
	fprintf(param_file, "                                            %% the starting redshift,\n");
	fprintf(param_file, "                                            %% otherwise this is recomputed\n");
	fprintf(param_file, "                                            %% based on the specified sigma8\n\n\n");



	fprintf(param_file, "ShapeGamma       0.201     %% only needed for Efstathiou power spectrum \n");
	fprintf(param_file, "PrimordialIndex  %f       %% may be used to tilt the primordial index\n\n\n", ns);


		  
	fprintf(param_file, "Seed             %d      %% seed for IC-generator\n\n\n", seed);  // % missing


	fprintf(param_file, "NumFilesWrittenInParallel 4  %% limits the number of files that are\n");
	fprintf(param_file, "                             %% written in parallel when outputting\n\n\n");




	fprintf(param_file, "UnitLength_in_cm          3.085678e21   %% defines length unit of output (in cm/h) \n");
	fprintf(param_file, "UnitMass_in_g             1.989e43      %% defines mass unit of output (in g/cm)\n");
	fprintf(param_file, "UnitVelocity_in_cm_per_s  1e5           %% defines velocity unit of output (in cm/sec)\n\n\n\n");


	fclose(param_file);

}


void write_Gadget_parameter_file(char Gadget_param_filename[], char simulation_codename[], char filebase2[], char output_list_filename[], double boxsize, double OBh2, double OM, double OL, double w0, double wa, double h, double z, double soft)
// Writes parameter file for Gadget-2.
{
	FILE* param_file;
	char filename[1000];
 
	sprintf(filename, "%s/%s/%s/%s/%s", ics_dir_speccomp, series_folder, Gadget_folder, Gadget_param_folder, Gadget_param_filename); 
	param_file=fopen(filename, "w");
		
		
	fprintf(param_file, "%%  Relevant files\n\n");

	fprintf(param_file, "InitCondFile  	   %s/%s/%s/%s/%s\n", ics_data_dir_simcomp, series_folder, Gadget_folder, Gadget_data_folder, filebase2);    
	fprintf(param_file, "OutputDir          %s/%s/%s/\n\n", Gadget_output_dir, series_folder, simulation_codename);

	fprintf(param_file, "EnergyFile         energy.txt\n");
	fprintf(param_file, "InfoFile           info.txt\n");
	fprintf(param_file, "TimingsFile        timings.txt\n");
	fprintf(param_file, "CpuFile            cpu.txt\n\n");

	fprintf(param_file, "RestartFile        restart\n"); // change this to something better.
	fprintf(param_file, "SnapshotFileBase   snapshot\n\n"); // change this to something better.

	fprintf(param_file, "OutputListFilename %s/%s/%s/%s/%s\n\n", ics_dir_simcomp, series_folder, Gadget_folder, Gadget_param_folder, output_list_filename);

	fprintf(param_file, "%% CPU time -limit\n\n");

	fprintf(param_file, "TimeLimitCPU      345600  %% = 4 days\n");
	fprintf(param_file, "ResubmitOn        0\n");
	fprintf(param_file, "ResubmitCommand   my-scriptfile\n");  


	fprintf(param_file, "%% Code options\n\n\n");


	fprintf(param_file, "ICFormat                 1\n");
	fprintf(param_file, "SnapFormat               1\n");
	fprintf(param_file, "ComovingIntegrationOn    1\n\n");

	fprintf(param_file, "TypeOfTimestepCriterion  0\n");
	fprintf(param_file, "OutputListOn             1\n");
	fprintf(param_file, "PeriodicBoundariesOn     1\n\n");

	fprintf(param_file, "%%  Caracteristics of run\n\n");

	fprintf(param_file, "TimeBegin           %f  %% Begin of the simulation (scale factor time, i.e. 1/(z+1))\n", 1.0/(z+1));
	fprintf(param_file, "TimeMax	            1.0\n\n");

	fprintf(param_file, "Omega0	              %f\n", OM);
	fprintf(param_file, "OmegaLambda           %f\n", OL);
	fprintf(param_file, "OmegaBaryon           %f\n", OBh2/(h*h));
	fprintf(param_file, "HubbleParam           %f\n", h);
	fprintf(param_file, "BoxSize               %f\n", boxsize*1000.0);
	fprintf(param_file, "w0                    %f\n", w0);
	fprintf(param_file, "wa                    %f\n\n", wa);

	fprintf(param_file, "%% Output frequency\n\n");

	fprintf(param_file, "TimeBetSnapshot        0.5\n");
	fprintf(param_file, "TimeOfFirstSnapshot    0\n\n");

	fprintf(param_file, "CpuTimeBetRestartFile     45000.0  %% = 12.5 hours   ; here in seconds\n");
	fprintf(param_file, "TimeBetStatistics         0.05\n\n");

	fprintf(param_file, "NumFilesPerSnapshot       16\n");
	fprintf(param_file, "NumFilesWrittenInParallel 8\n\n\n\n");



	fprintf(param_file, "%% Accuracy of time integration\n\n");

	fprintf(param_file, "ErrTolIntAccuracy      0.025 \n\n");

	fprintf(param_file, "MaxRMSDisplacementFac  0.2\n\n");

	fprintf(param_file, "CourantFac             0.15     \n\n");

	fprintf(param_file, "MaxSizeTimestep       0.02\n");
	fprintf(param_file, "MinSizeTimestep       0.0\n\n\n\n\n");




	fprintf(param_file, "%% Tree algorithm, force accuracy, domain update frequency\n\n");

	fprintf(param_file, "ErrTolTheta            0.45\n");            
	fprintf(param_file, "TypeOfOpeningCriterion 1\n");
	fprintf(param_file, "ErrTolForceAcc         0.005\n\n\n");


	fprintf(param_file, "TreeDomainUpdateFrequency    0.025\n\n\n");


	fprintf(param_file, "%%  Further parameters of SPH\n\n");

	fprintf(param_file, "DesNumNgb              33\n");
	fprintf(param_file, "MaxNumNgbDeviation     2\n");
	fprintf(param_file, "ArtBulkViscConst       0.8\n");
	fprintf(param_file, "InitGasTemp            1000.0        %% always ignored if set to 0 \n");
	fprintf(param_file, "MinGasTemp             50.0    \n\n\n");


	fprintf(param_file, "%% Memory allocation\n\n");

	fprintf(param_file, "PartAllocFactor       1.3\n");
	fprintf(param_file, "TreeAllocFactor       0.7\n");
	fprintf(param_file, "BufferSize            20          %% in MByte\n\n\n");


	fprintf(param_file, "%% System of units\n\n");

	fprintf(param_file, "UnitLength_in_cm         3.085678e21        ;  1.0 kpc \n");
	fprintf(param_file, "UnitMass_in_g             1.989e43    ;  1.0e10 solar masses \n"); // obsolete from Edwin Sirko's IC generator: 3.94677661169e34 %%correcting for ic-generator
	fprintf(param_file, "UnitVelocity_in_cm_per_s 1e5                ;  1 km/sec \n");
	fprintf(param_file, "GravityConstantInternal  0\n\n\n");
 

	fprintf(param_file, "%% Softening lengths\n\n");

	fprintf(param_file, "MinGasHsmlFractional 0.25\n\n");

	fprintf(param_file, "SofteningGas       0\n");
	fprintf(param_file, "SofteningHalo      %f\n", soft);
	fprintf(param_file, "SofteningDisk      0\n");
	fprintf(param_file, "SofteningBulge     0\n");           
	fprintf(param_file, "SofteningStars     0\n");
	fprintf(param_file, "SofteningBndry     0\n\n");

	fprintf(param_file, "SofteningGasMaxPhys       0\n");
	fprintf(param_file, "SofteningHaloMaxPhys      %f\n", soft);
	fprintf(param_file, "SofteningDiskMaxPhys      0\n");
	fprintf(param_file, "SofteningBulgeMaxPhys     0\n");      
	fprintf(param_file, "SofteningStarsMaxPhys     0\n");
	fprintf(param_file, "SofteningBndryMaxPhys     0\n");

	fclose(param_file);
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// NYBlue job description and submission:

void write_BGL_description(char NGenIC_param_filename[], char Gadget_param_filename[], char jobstamm[], int job_nr)
{
	FILE *description_file;
	char jobname[1000];
	
	sprintf(jobname, "%s/%s/%s/%s/%s_%d.ll", ics_dir_speccomp, series_folder, NGenIC_folder, NGenIC_jobs_folder, jobstamm, job_nr);
	description_file=fopen(jobname, "w");
	
		// Write job description file for N-GenIC for NYBlue/L (Blue Gene/L at BNL):
	////////////////////////////////////////
	
	fprintf(description_file, "# @ job_type = bluegene\n");
	fprintf(description_file, "# @ class = long\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Classes above are currently: short (1024, 2048, 3072, or 4096 nodes, 24h) , normaldyn: (512 nodes, 48h), long (32 or 128 nodes, 72h).\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# The executable that will run your parallel application should always be specified as per the next line.\n");
	fprintf(description_file, "# @ executable = /bgl/BlueLight/ppcfloor/bglsys/bin/mpirun\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Run on 512 nodes using a dynamic partition.\n");
	fprintf(description_file, "# Specify partition size using the following statement. This statement is the only way a partition size \n");
	fprintf(description_file, "# should ever be specified in a LoadLeveler job control file, i.e. use of bg_partition\n");
	fprintf(description_file, "# has been eliminated.\n");
	fprintf(description_file, "# It is possible to run on fewer processors than those afforded by the partition size, see the \n");
	fprintf(description_file, "# description of -np in the \"Notes for Sample\" below. BUT NOTE THAT you must use # @ bg_size \n");
	fprintf(description_file, "# and not -np to specify the partition size.  In other words, you must use # @ bg_size\n");
	fprintf(description_file, "# to allocate the partition. \n");
	fprintf(description_file, "# Then you can optionally use -np to run on fewer processors if this is necessary \n");
	fprintf(description_file, "# for your run -- but you will charged for the entire partition that you allocated.\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# @ bg_size = 128\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# initialdir (see the next active line) will be used as the working directory for this batch job. \n");
	fprintf(description_file, "# @ initialdir = %s/%s/%s/%s \n", ics_dir_simcomp, series_folder, NGenIC_folder, NGenIC_logs_folder);
	fprintf(description_file, "#\n");
	fprintf(description_file, "# If for example your jobid is 82, your output and error will be written in\n");
	fprintf(description_file, "# directory /home/johndoe/app1/runs, to files 82.out and 82.err respectively.\n");
	fprintf(description_file, "# @ input = /dev/null\n");
	fprintf(description_file, "# @ output = $(jobid).out\n");
	fprintf(description_file, "# @ error = $(jobid).err\n");
	fprintf(description_file, "#\n"); 
	fprintf(description_file, "# Maximum wall clock time for job will be 20 minutes.\n");
	fprintf(description_file, "# @ wall_clock_limit = 00:20:00\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Send email to johndoe@bnl.gov when job has completed.\n");
	fprintf(description_file, "# @ notification = complete\n");
	fprintf(description_file, "# @ notify_user = jank@astro.columbia.edu\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Specify executable for your parallel application, and arguments to that executable.\n");
	fprintf(description_file, "# Note that the arguments to specify for the executable will vary depending upon the executable.\n");
	fprintf(description_file, "# Specify any special environment variables for your application that need to be in the environment \n");
	fprintf(description_file, "# presented to the job on the compute nodes, they will vary\n");
	fprintf(description_file, "# depending upon the application -  some applications will not require any - so delete or modify the\n");
	fprintf(description_file, "# -env specification below.\n");
	fprintf(description_file, "# @ arguments = -np 256 -exe %s/%s \\ \n", NGenIC_dir_L, NGenIC_exec);
	fprintf(description_file, "-cwd %s/%s/%s/%s \\ \n", ics_dir_simcomp, series_folder, NGenIC_folder, NGenIC_param_folder);
	fprintf(description_file, "-mode VN \\ \n");
	fprintf(description_file, "-args \"%s\" \n", NGenIC_param_filename);
	fprintf(description_file, "#\n");
	fprintf(description_file, "# The next statement marks the end of the job step. This example is a one-job-step batch job,\n");
	fprintf(description_file, "# so this is equivalent to saying that the next statement marks the end of the batch job.\n");
	fprintf(description_file, "# @ queue\n");

		
	////////////////////////////////////////
	
	fclose(description_file);
		
	sprintf(jobname, "%s/%s/%s/%s/%s_%d.ll", ics_dir_speccomp, series_folder, Gadget_folder, Gadget_jobs_folder, jobstamm, job_nr);
	description_file=fopen(jobname, "w");
	
		// Write job description file for Gadget-2 NYBlue/L (Blue Gene/L at BNL):
		
	////////////////////////////////////////
	
	fprintf(description_file, "# @ job_type = bluegene\n");
	fprintf(description_file, "# @ class = long\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Classes above are currently: short (1024, 2048, 3072, or 4096 nodes, 24h) , normaldyn: (512 nodes, 48h), long (32 or 128 nodes, 72h).\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# The executable that will run your parallel application should always be specified as per the next line.\n");
	fprintf(description_file, "# @ executable = /bgl/BlueLight/ppcfloor/bglsys/bin/mpirun\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Run on 512 nodes using a dynamic partition.\n");
	fprintf(description_file, "# Specify partition size using the following statement. This statement is the only way a partition size \n");
	fprintf(description_file, "# should ever be specified in a LoadLeveler job control file, i.e. use of bg_partition\n");
	fprintf(description_file, "# has been eliminated.\n");
	fprintf(description_file, "# It is possible to run on fewer processors than those afforded by the partition size, see the \n");
	fprintf(description_file, "# description of -np in the \"Notes for Sample\" below. BUT NOTE THAT you must use # @ bg_size \n");
	fprintf(description_file, "# and not -np to specify the partition size.  In other words, you must use # @ bg_size\n");
	fprintf(description_file, "# to allocate the partition. \n");
	fprintf(description_file, "# Then you can optionally use -np to run on fewer processors if this is necessary \n");
	fprintf(description_file, "# for your run -- but you will charged for the entire partition that you allocated.\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# @ bg_size = 128\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# initialdir (see the next active line) will be used as the working directory for this batch job. \n");
	//	printf(description_file, "# @ initialdir = /gpfs/home2/jank/Documents/Gadget/noGSL-w-Gadget-2.0.3/Gadget2\n");
	fprintf(description_file, "# @ initialdir = %s/%s/%s/%s \n", ics_dir_simcomp, series_folder, Gadget_folder, Gadget_logs_folder);
	fprintf(description_file, "#\n");
	fprintf(description_file, "# If for example your jobid is 82, your output and error will be written in\n");
	fprintf(description_file, "# directory /home/johndoe/app1/runs, to files 82.out and 82.err respectively.\n");
	fprintf(description_file, "# @ input = /dev/null\n");
	fprintf(description_file, "# @ output = $(jobid).out\n");
	fprintf(description_file, "# @ error = $(jobid).err\n");
	fprintf(description_file, "#\n"); 
	fprintf(description_file, "# Maximum wall clock time for job will be 72 hours.\n");
	fprintf(description_file, "# @ wall_clock_limit = 72:00:00\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Send email to johndoe@bnl.gov when job has completed.\n");
	fprintf(description_file, "# @ notification = complete\n");
	fprintf(description_file, "# @ notify_user = jank@astro.columbia.edu\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Specify executable for your parallel application, and arguments to that executable.\n");
	fprintf(description_file, "# Note that the arguments to specify for the executable will vary depending upon the executable.\n");
	fprintf(description_file, "# Specify any special environment variables for your application that need to be in the environment \n");
	fprintf(description_file, "# presented to the job on the compute nodes, they will vary\n");
	fprintf(description_file, "# depending upon the application -  some applications will not require any - so delete or modify the\n");
	fprintf(description_file, "# -env specification below.\n");
	fprintf(description_file, "# @ arguments = -np 256 -exe %s/%s \\ \n", Gadget_dir_L, Gadget_exec);
	// printf(description_file, "-cwd /gpfs/home2/jank/Documents/Gadget/noGSL-w-Gadget-2.0.3/Gadget2 \ \n");
	fprintf(description_file, "-cwd %s/%s/%s/%s \\ \n", ics_dir_simcomp, series_folder, Gadget_folder, Gadget_param_folder);
	fprintf(description_file, "-mode VN \\ \n");
	// printf(description_file, "-args \"%s/%s/%s/%s/%s\" \n", ics_dir_simcomp, series_folder, Gadget_folder, Gadget_param_folder, Gadget_param_filename);
	fprintf(description_file, "-args \"%s\" \n", Gadget_param_filename);
	fprintf(description_file, "#\n");
	fprintf(description_file, "# The next statement marks the end of the job step. This example is a one-job-step batch job,\n");
	fprintf(description_file, "# so this is equivalent to saying that the next statement marks the end of the batch job.\n");
	fprintf(description_file, "# @ queue\n");

	////////////////////////////////////////

	fclose(description_file);

}

void write_BGP_description(char **BGP_NGenIC_param_filename, char **BGP_Gadget_param_filename, char jobstamm[], int job_nr, int octo_counter)
{

	FILE *description_file;
	char jobname[1000], arguments[1000], temp[1000];
	int i;
	
	for (i=0; i<octo_counter; i++) // make individual jobs for N-GenIC, not Octopus.
	{
		sprintf(jobname, "%s/%s/%s/%s/%s_%d.ll", ics_dir_speccomp, series_folder, NGenIC_folder, NGenIC_jobs_folder, jobstamm, (job_nr-1)*8+i);
		description_file=fopen(jobname, "w");
	
		// Write job description file for N-GenIC for NYBlue/L (Blue Gene/L at BNL):
		////////////////////////////////////////

	fprintf(description_file, "# @ job_type = bluegene\n");
	fprintf(description_file, "# @ class = normal\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# The executable that will run your parallel application should always be specified as per the next line.\n");
	fprintf(description_file, "# @ executable = /usr/bin/mpirun\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Run job on 512 compute nodes. LoadLeveler will dynamically allocate the partition.\n");
	fprintf(description_file, "# @ bg_size = 128\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# initialdir will be the initial directory. LoadLeveler changes to this\n");
	fprintf(description_file, "# directory before running the job. If you don't specify this, it defaults to your current working directory\n");
	fprintf(description_file, "# at the time the job was submitted.\n");
	fprintf(description_file, "# File names mentioned in the batch script that do not begin with a slash ( / ) are relative to the initial\n");
	fprintf(description_file, "# directory.\n");
	fprintf(description_file, "# The initial directory must exist on both the fen and the compute nodes.\n");
	fprintf(description_file, "# @ initialdir = %s/%s/%s/%s\n", ics_dir_simcomp, series_folder, NGenIC_folder, NGenIC_logs_folder);
	fprintf(description_file, "#\n");
	fprintf(description_file, "# If for example your jobid is 82, your output and error will be written in\n");
	fprintf(description_file, "# directory /home/johndoe/app1/runs, to files 82.out and 82.err respectively.\n");
	fprintf(description_file, "# @ input = /dev/null\n");
	fprintf(description_file, "# @ output = $(jobid).out\n");
	fprintf(description_file, "# @ error = $(jobid).err\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Maximum wall clock time for job will be 20 minutes.\n");
	fprintf(description_file, "# @ wall_clock_limit = 00:20:00\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Send email to johndoe@bnl.gov when job has completed.\n");
	fprintf(description_file, "# @ notification = complete\n");
	fprintf(description_file, "# @ notify_user = jank@astro.columbia.edu\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Specify executable for your parallel application, and arguments to that executable.\n");
	fprintf(description_file, "# Note that the arguments to specify for the executable will vary depending upon the executable.\n");
	fprintf(description_file, "# Specify any special environment variables for your application that need to be in the environment \n");
	fprintf(description_file, "# presented to the job on the compute nodes, they will vary\n");
	fprintf(description_file, "# depending upon the application -  some applications will not require any - so delete or modify the\n");
	fprintf(description_file, "# -env specification below.\n");
	fprintf(description_file, "# @ arguments = -np 256 -exe %s/%s \\ \n", NGenIC_dir_P, NGenIC_exec);
	fprintf(description_file, "-cwd %s/%s/%s/%s \\ \n", ics_dir_simcomp, series_folder, NGenIC_folder, NGenIC_param_folder);
	fprintf(description_file, "-mode DUAL \\ \n");
	fprintf(description_file, "-args \"%s\" \n", BGP_NGenIC_param_filename[i]);
	fprintf(description_file, "#\n");
	fprintf(description_file, "# The next statement marks the end of the job step. This example is a one-job-step batch job,\n");
	fprintf(description_file, "# so this is equivalent to saying that the next statement marks the end of the batch job.\n");
	fprintf(description_file, "# @ queue\n");

		////////////////////////////////////////
	
		fclose(description_file);
	}
	
					
	sprintf(jobname, "%s/%s/%s/%s/%s_%d.ll", ics_dir_speccomp, series_folder, Gadget_folder, Gadget_jobs_folder, jobstamm, job_nr);
	description_file=fopen(jobname, "w");
	
		// Write job description file for Gadget-2 NYBlue/L (Blue Gene/L at BNL):
		
	////////////////////////////////////////
	
	fprintf(description_file, "# @ job_type = bluegene\n");
	fprintf(description_file, "# @ class = normal\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# The executable that will run your parallel application should always be specified as per the next line.\n");
	fprintf(description_file, "# @ executable = /usr/bin/mpirun\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Run job on 512 compute nodes. LoadLeveler will dynamically allocate the partition.\n");
	fprintf(description_file, "# @ bg_size = %d\n", 64*octo_counter); // each simulation runs on 64 Blue Gene/P nodes (on all of their 256 CPUs)
	fprintf(description_file, "#\n");
	fprintf(description_file, "# initialdir will be the initial directory. LoadLeveler changes to this\n");
	fprintf(description_file, "# directory before running the job. If you don't specify this, it defaults to your current working directory\n");
	fprintf(description_file, "# at the time the job was submitted.\n");
	fprintf(description_file, "# File names mentioned in the batch script that do not begin with a slash ( / ) are relative to the initial\n");
	fprintf(description_file, "# directory.\n");
	fprintf(description_file, "# The initial directory must exist on both the fen and the compute nodes.\n");
	fprintf(description_file, "# @ initialdir = %s/%s/%s/%s\n", ics_dir_simcomp, series_folder, Gadget_folder, Gadget_logs_folder);
	fprintf(description_file, "#\n");
	fprintf(description_file, "# If for example your jobid is 82, your output and error will be written in\n");
	fprintf(description_file, "# directory /home/johndoe/app1/runs, to files 82.out and 82.err respectively.\n");
	fprintf(description_file, "# @ input = /dev/null\n");
	fprintf(description_file, "# @ output = $(jobid).out\n");
	fprintf(description_file, "# @ error = $(jobid).err\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Maximum wall clock time for job will be 48 hours.\n");
	fprintf(description_file, "# @ wall_clock_limit = 48:00:00\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Send email to johndoe@bnl.gov when job has completed.\n");
	fprintf(description_file, "# @ notification = complete\n");
	fprintf(description_file, "# @ notify_user = jank@astro.columbia.edu\n");
	fprintf(description_file, "#\n");
	fprintf(description_file, "# Specify executable for your parallel application, and arguments to that executable.\n");
	fprintf(description_file, "# Note that the arguments to specify for the executable will vary depending upon the executable.\n");
	fprintf(description_file, "# Specify any special environment variables for your application that need to be in the environment \n");
	fprintf(description_file, "# presented to the job on the compute nodes, they will vary\n");
	fprintf(description_file, "# depending upon the application -  some applications will not require any - so delete or modify the\n");
	fprintf(description_file, "# -env specification below.\n");
	fprintf(description_file, "# @ arguments = -np %d -exe %s/%s \\ \n", 256*octo_counter, Gadget_dir_P, Gadget_exec);
	fprintf(description_file, "-cwd %s/%s/%s/%s \\ \n", ics_dir_simcomp, series_folder, Gadget_folder, Gadget_param_folder);
	fprintf(description_file, "-mode VN \\ \n");
	sprintf(arguments, "%d 256", octo_counter); // octo_counter-many simulations (typically 8), each on 256 processors. 
	for (i=0; i<octo_counter; i++)
	{
		sprintf(temp, "%s %s", arguments, BGP_Gadget_param_filename[i]);
		sprintf(arguments, "%s", temp);
	}
	fprintf(description_file, "-args \"%s\" \n", arguments); // 8 simulations on 256 processors each.
	fprintf(description_file, "#\n");
	fprintf(description_file, "# The next statement marks the end of the job step. This example is a one-job-step batch job,\n");
	fprintf(description_file, "# so this is equivalent to saying that the next statement marks the end of the batch job.\n");
	fprintf(description_file, "# @ queue\n");
	
	////////////////////////////////////////

	fclose(description_file);

}
				
void write_submission_script_L(int jobs, char jobstamm[], char BG_type[])
{
	FILE *script_file;
	char filename[1000], jobname[200];
	int i, app;

	for (app=0; app<=1; app++)
	{
		if (app==0) sprintf(filename, "%s/%s/%s/subm-script_N-GenIC_%s_%s.sh", ics_dir_speccomp, series_folder, NGenIC_folder, series_folder, BG_type);
		else sprintf(filename, "%s/%s/%s/subm-script_Gadget_%s_%s.sh", ics_dir_speccomp, series_folder, Gadget_folder, series_folder, BG_type);
		script_file=fopen(filename, "w");
	
		fprintf(script_file, "#!/bin/bash\n\n");
		fprintf(script_file, "# Shell script for submission of initial condition generation for simulation series on Blue Gene/%s.\n\n\n", BG_type); 
		for (i=1; i<=jobs; i++)
		{
			sprintf(jobname, "%s_%d.ll", jobstamm, i);
			fprintf(script_file, "echo \"Submitting job: %s ...\"\n", jobname);
			if (app==0) fprintf(script_file, "llsubmit %s/%s\n", NGenIC_jobs_folder, jobname);
			else fprintf(script_file, "llsubmit %s/%s\n", Gadget_jobs_folder, jobname);
			fprintf(script_file, "sleep 5\n");
		}
	
		fprintf(script_file, "echo \"Done submitting all jobs.\"\n");
	
// Scripting example to wait for process to finish (include such structure if necessary):
// chessboard >/dev/null & EPID=$!
// wait $EPID
// echo "Chessboard has finished"
// echo "Process ID was: Pid = $EPID"
	
		fclose(script_file);
	}
}

void write_submission_script_P(int jobs, int Pjobs, char jobstamm[], char BG_type[])
{
	FILE *script_file;
	char filename[1000], jobname[200];
	int i, app, n;

	for (app=0; app<=1; app++)
	{
		if (app==0) sprintf(filename, "%s/%s/%s/subm-script_N-GenIC_%s_%s.sh", ics_dir_speccomp, series_folder, NGenIC_folder, series_folder, BG_type);
		else sprintf(filename, "%s/%s/%s/subm-script_Gadget_%s_%s.sh", ics_dir_speccomp, series_folder, Gadget_folder, series_folder, BG_type);
		script_file=fopen(filename, "w");
	
		fprintf(script_file, "#!/bin/bash\n\n");
		fprintf(script_file, "# Shell script for submission of initial condition generation for simulation series on Blue Gene/%s.\n\n\n", BG_type); 
		if (app==0) n=jobs;
		else n=Pjobs;
		for (i=1; i<=n; i++)
		{
			sprintf(jobname, "%s_%d.ll", jobstamm, i);
			fprintf(script_file, "echo \"Submitting job: %s ...\"\n", jobname);
			if (app==0) fprintf(script_file, "llsubmit %s/%s\n", NGenIC_jobs_folder, jobname);
			else fprintf(script_file, "llsubmit %s/%s\n", Gadget_jobs_folder, jobname);
			fprintf(script_file, "sleep 5\n");
		}
	
		fprintf(script_file, "echo \"Done submitting all jobs.\"\n");
	
		fclose(script_file);
	}
}

void write_directory_script(FILE *script_file, char simulation_codename[])
{
  fprintf(script_file, "   if [ ! -d %s/%s/%s ] ; then\n       mkdir %s/%s/%s\n   fi\n\n", Gadget_output_dir, series_folder, simulation_codename, Gadget_output_dir, series_folder, simulation_codename);
}


void write_IG_directory_script(FILE *script_file, char simulation_codename[])
{
  fprintf(script_file, "   if [ ! -d %s/%s/%s ] ; then\n       mkdir %s/%s/%s\n   fi\n\n", IG_output_dir, series_folder, simulation_codename, IG_output_dir, series_folder, simulation_codename);
  fprintf(script_file, "   if [ ! -d %s/%s/%s/%s ] ; then\n       mkdir %s/%s/%s/%s\n   fi\n\n", IG_output_dir, series_folder, simulation_codename, IG_planes_folder, IG_output_dir, series_folder, simulation_codename, IG_planes_folder);
  fprintf(script_file, "   if [ ! -d %s/%s/%s/%s ] ; then\n       mkdir %s/%s/%s/%s\n   fi\n\n", IG_output_dir, series_folder, simulation_codename, IG_maps_folder, IG_output_dir, series_folder, simulation_codename, IG_maps_folder);
  fprintf(script_file, "   if [ ! -d %s/%s/%s/%s ] ; then\n       mkdir %s/%s/%s/%s\n   fi\n\n", IG_output_dir, series_folder, simulation_codename, IG_products_nonoise_folder, IG_output_dir, series_folder, simulation_codename, IG_products_nonoise_folder);
  fprintf(script_file, "   if [ ! -d %s/%s/%s/%s ] ; then\n       mkdir %s/%s/%s/%s\n   fi\n\n\n", IG_output_dir, series_folder, simulation_codename, IG_products_noise_folder, IG_output_dir, series_folder, simulation_codename, IG_products_noise_folder);

}

void write_simple_list(FILE *script_file, char simulation_codename[])
{
  fprintf(script_file, "%s\n", simulation_codename);
}



int fgetline(FILE *stream, char s[], int lim)
{
	int c, i, j;

	for (i=0;i<lim-1 && (c=getc(stream))!=EOF && c!='\n' && c!='#'; ++i) s[i]=c;
	if (c == '\n' || c == '#')
	{
		s[i] = '\n';
		i++;
	}
	if (c == '#') {for (j=0;j<lim-1-i && (c=getc(stream))!=EOF && c!='\n'; ++j);};
	s[i]='\0';
	return i;
}
