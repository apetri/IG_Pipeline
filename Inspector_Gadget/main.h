/*
 *  main.h
 *  Inspector Gadget, Version 5.3
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 5/8/2007.
 *  Later improved at the University of Miami and the University of KwaZulu-Natal.
 *  Last updated on April 4, 2014.
 *  Copyright 2007-14 Jan Michael Kratochvil. All rights reserved.
 *  This is not a public release.
 *
 */

///////////////////////////////////////////
// MANUAL SELECTION:
///////////////////////////////////////////

// Use if want MPI compilation (e.g. for Blue Gene/P), comment out if want serial compilation for submission with Condor scheduling system (e.g. on LSST cluster):
#define MPI_COMPILE 1
// NOTE: The current code must have MPI enables. Running it without MPI is a depreciated option.

// Use this line if want to have OpenMP shared memory threading in addition to distributed memory MPI parallelization:
//#define HAVE_OpenMP
// The number of OpenMP threads used by the code will depend on the setting of the OpenMP environment variable OMP_NUM_THREADS on the system on which the code is run.

// Use this line to define the file output format. Currently only FITS files (cfitsion library) are supported:
#define HAVE_FITS

// For best performance, use all lines above (do not outcomment any) - except for OpenMP for Mode 1 (potential plane generation), where is may cause a (potentially serious) performance degradation.

///////////////////////////////////////////
///////////////////////////////////////////

#include <time.h>


#define pi 3.141592653589793
#define Gravi_const 1.0
// Gravitational constant

// Maximal allowable line length for read-in from input file:
#define MAXLINE 2000
// Maximal allowable path name length for paths from input file:
#define MAXPATHNAME 2000
// Maximal allowable name length for parameters and file names in input file:
#define MAXNAME 200
// The above values can be modified as needed for longer line and name lengths.
// WARNING: Make sure you also change the corresponding string array lengths in the struct analysis_parameters parameters in the file main.h!


#if defined(MPI_COMPILE)
#include <mpi.h>
extern MPI_Comm sim_comm;
extern MPI_Comm snapshot_comm; // USED??
extern MPI_Comm map_comm;
#endif
extern int superrank;
extern int supersize;

extern int **cosmology_sampler;

extern double UnitMass_in_g;

/*
// Old classic version of Gadget-2 snapshot header:
struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  // fills to 256 Bytes
} ;
*/


struct io_header_1
{
  int npart[6];                        /*!< number of particles of each type in this file */
  double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int flag_sfr;                        /*!< flags whether the simulation was including star formation */
  int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                    /*!< flags whether cooling was included  */
  int num_files;                       /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  // <JMK>:                                                                                                 
  double w0;                           /*!< Dark energy equation of state parameters (not part of standard Gadget-2 release) */
  double wa;
  double comoving_distance;            /*!< Comoving distance of snapshot from observer at z=0 (a=1) (in h^-1 kpc) */
  // </JMK>.              
  // Variables below are not used in Inspector Gadget:                                                                                  
  int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  // <JMK>:                                                                                                 
  int nothing;                         /*!< Meaningless safety integer for better cross-platform byte alignment of structure (not defined in C standard) */
  char fill[32];  // original: char fill[60];                  /*!< fills to 256 Bytes */                   
  // </JMK>.                                                                                                
} ;
							           /*!< holds header for snapshot files */
extern struct io_header_1 header1;



struct particle_data 
{
  float  Pos[3];
  //  float  Vel[3];
  float  Mass;
  int    Type;

  //  float  Rho, U, Temp, Ne;
} ;

extern struct particle_data *P;


struct snapshot_control
{
	struct particle_data *P;
	int *Id;
	int NumPart;
	int Ngas;
	int Nhalo;
	double mass[6];
	double Time;
	double Redshift;
	int NumPartTotal[6];
	int files;
	double Boxsize;
	double Omega0;
	double OmegaLambda;
	double H_0;
	
	double w0;
	double wa;
	double ns;
	double sigma_8;

    double comoving_distance;
    double initial_condition;

    int set; // Flag if values in this structure have been set properly for a snapshot number (set to 1 if yes).
	
} ;

extern struct snapshot_control *Snapshot;


struct particle_analysis 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;
  int	 write;
  float  Color[3];
  int    Bin;
  //float  Rho, U, Temp, Ne;
} ;



struct analysis_parameters
{ 
  int feedback;
  int species;
  int NbinsX;		// number of weak lensing evaluation bins (traced lightrays) on sky, in x direction (horizontal). 
  int NbinsY;		// number of weak lensing evaluation bins (traced lightrays) on sky, in x direction (vertical). 
  int nx;			// number of grid points on 2D lens planes for Fourier transforms (Poisson equation), in x direction. Reasonable, fast number is 512, can be higher though.
  int ny;			// number of grid points on 2D lens planes for Fourier transforms (Poisson equation), in y direction. Reasonable, fast number is 512, can be higher though.
  int nxny;			// nxny = nx*ny, is a derived parameter from the above. 
  int scramble_mode;
  int plane_assignment_mode;
  int averaging_mode;
  int ray_tracing;
  float survey_angle;
  float survey_angle_in_rad;
  double plane_padding;
  int snapshots; // corrected from float
  int number_of_planes;
  float	boxsize;
  int generate_planes; // if 0, will read in pre-fabricated potential planes into ray-tracing code, rather than analyzing a box-ensemble from an N-body simulation.
  int generate_maps; // if 1, will generate weak lensing maps.
  int plane_shift;
  int seed;
  int cell_embedding;
  int raypoint_averaging;
  int first_plane; // plane_number of first (closest) plane to be used (depending on mode).
  int last_plane; // plane_number of last (farthest) plane to be used (depending on mode).
  
  int parallel;
  int process_number;
  int number_of_seeds;
  int seed_block;
  int realization;
  
  int first_snapshot; // file number of first snapshot to be processed by one process (on one CPU).
  int last_snapshot; // file number of last snapshot to be processed by one process (on one CPU).
  int number_of_snapshots; // this is the number of snapshots processed locally by one process (on one CPU).
  int global_first_snapshot; // file number of first snapshot processed by any process (under parallel).
  int global_last_snapshot; // file number of last snapshot processed by any process (under parallel).
  // the number of global snapshots is called "parameters.snapshots" and declared further above already.
  int first_realization; // realization number of first realization processed by one process (on one CPU).
  int last_realization; // realization number of last realization processed by one process (on one CPU).
  int number_of_realizations; // number of realizations processed by one process (on one CPU).
  // for realizations, the global range doesn't need to be indicated (unlike with snapshots), as they are completely independent of each other.
  int global_first_realization; // new since IG X-3.0 for Mode=2: first realization number processed by any process.
  int global_last_realization; // new since IG X-3.0 for Mode=2: first realization number processed by any process.
  int fiducial; // new since IG X-3.0 for Mode=2: selects if fiducial cosmology (influences how many different ICs for the sims are taken into account).

  // Cosmoogical parameters:
  double H_0; // Hubble constant today (in km/s/Mpc).
  double h; // Hubble parameter (H_0=100*h km/s/Mpc).
  double Omega_m;
  double Omega_Lambda;
  double w0;
  double wa;
  double ns;
  double sigma_8;
  double initial_condition;

  // Survey parameters (written into WL map FITS headers for easier reference and post-analysis):
  double source_redshift; // redshift of sources (assuming delta distribution or mean)
  double source_scale_factor;
  double source_comoving_distance;
  double ASPP_resolution; // Resolution of WL map: arcseconds/pixel.
  double PPAM_resolution; // Resolution of WL map: pixels/arcminute.
  
  double NumPartTotal[6]; // Total number of particles of each species in each simulation box of simulation.
  double mass[6]; // masses of particle species.
  
  // Path names (full or relative paths):
  char snapshot_path[2000]; // includes folder in which snapshots of the simulation are stored.
  char parameter_path[2000];
  char plane_path[2000];
  char plane_output_path[2000];
  char map_output_path[2000];
  char plane_urpath[2000]; // since IG X-3.0 for Mode=2: plane path without simulation folder (since the name of the latter is assembled on an ad hoc basis.
  
  // Parameter file names:
  char parameterfile1[200];
  char parameterfile2[200];
  
  // Read and write out file name components:
  char modelname[200]; // description of cosmological model
  char simulation_codename[200]; // this code will be written in every output file associated with the simulation.
  char snapshot_name[200];
  char density_basename[200];
  char potential_basename[200];
  
  char convergence_basename[200];
  char shear1_basename[200];
  char shear2_basename[200];
  char omega_basename[200];
  char shear_modulus_basename[200];
  char shear_angle_basename[200];
  char deflection1_basename[200];
  char deflection2_basename[200];
  char deflection_total_basename[200];
  char deflection_winkel_basename[200]; // WARNING: This is not the "deflection angle", this is the angle in the 2D map of the deflection direction (deflection_total) giving its magnitude).
  
  char plane_comment[200];
  char WL_map_comment[200];
  
  //clock_t start_time; // Start time (for timing of run).
  time_t start_time;
  float runtime; // stores time of timing routine call to next call (to compute time difference between calls).

  // Condor paths:
  char condor_script_path[1000];
  char condor_script_name[200];
  int condor_script_mode;
  char storage_node_name[200];
  char storage_Gadget_path[2000];
  char local_Gadget_path[2000];
  char storage_IG_path[2000];
  char local_IG_path[2000];
  char map_basename[200];
  char comoving_file[200];
  char simfolder[200];
  char planes_folder[200];
  char maps_folder[200];

  int snapskip; // contains how many snapshots to be skipped between planes (no skip for snapskip=1).
  int number_of_boxcenters;

  // For source planes:
  int plane_before_source; // Is set to -1 (or generally <0 ) most of the time, only to calculate the plane before source plane, terminating at the source plane, it is set to the plane number of the plane before source. This is done by the condor script and sent to Inspector Gadget as one of the input parameters.

  int convergence_direct;

  // Parameters needed only for MPI version:
  int mode;  // Mode of Inspector Gadget: 1: Generate Potential Planes, 2: Generate WL Maps.
  int number_of_source_planes; // Number of source planes from which to choose to create WL maps.

  // MPI-specific parameters (to control and check MPI run):
  int ThisTask; // Process number for MPI run (same as parameters.process_number)
  int NTasks; // Number of requested MPI processes (DIFFERENT from parameters.parallel, which gives the number of processes NEEDED to achieve requested submission task; consequently NTasks>=parallel has to hold for any MPI submission, otherwise error).

  // Endianness:
  int byteswap;   // if set to !=0, bytes are swapped upon reading in of Gadget-2 snapshots (and other binaries where necessary) to compensate for different endianness.
  int endianness_set; // records if endiannes has been determined.

  int diagnostic;

  // For MPI with multiple cosmologies in parallel:
  int number_of_cosmologies;
  int cosmology_number;
  int number_of_processes;

  // Plane and Map file extension:
  char extension[20];

  int darkenergy_initialized;
  int chi_initialized;
  int snapshot_allocated;
    
  int pass; // temporary.
    
    // For galaxy catalogues (like CFHT, LSST, etc.):
    int galaxy_catalogue_type; // 0: regular maps, 1 (and higher): actual galaxy catalogues.
    char galaxy_catalogue_path[2000];
    char galaxy_catalogue_output_path[2000];
    char galaxy_catalogue_basename[2000];
    int galaxy_subfield; // subfield number, counts through galaxy subfields during ray tracing.
    int first_galaxy_subfield; // first galaxy subfield to be processed, input parameter.
    int last_galaxy_subfield; // last galaxy subfield to be processed.
    
    // Files containing random numbers determining which plane realization and ic to choose during ray tracing (and also how to shift and rotate this plane before ray tracing if desired):
    char plane_randomizer_file_general[2000];
    char plane_randomizer_file_fiducial[2000];
    // File containing rotations and shifts for N-body simulation snapshot boxes during lens plane generation:
    char snapshot_rotation_randomizer_file[2000];
    
    int max_realizations; // independent block in random number file.
    
    // For one-sided RMA MPI communication (can be later removed and automatized by inferring from scan of random number file):
    int preload_planes;
    int number_of_plane_realizations;
    int number_of_sim_ics;
    int first_sim_ic;

    
    
} ; 

extern struct analysis_parameters parameters;


extern int feedback;

void allocate_snapshot_memory(int NumPart, int snapshot_number);
void free_snapshot_arrays(int snapshot_number);
void free_all_snapshots(void);

struct scramble
{
	int x1;
	int x2;
	int x3;
	int signx1;
	int signx2;
	int signx3;
	int centeroffset1;
	int centeroffset2;
	int centeroffset3;
	float crosshairs_x2;
	float crosshairs_x3;
} ;


void load_snapshot_multi(char *fname, int file_number, int snapshot_number, int plane_number); // Note: plane number is not being used anymore in this function (could be removed from list of arguments).
void load_snapshot_multi_header(char *fname, int file_number, int snapshot_number);

// Unused functions (included for possible later use):
void unit_conversion(int snapshot_number);
void reordering(int snapshot_number);

extern FILE *input_file_ran;


