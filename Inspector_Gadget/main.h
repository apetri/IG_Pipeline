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
//<AP> inserted appropriate compile flag in the Makefile
//#define HAVE_OpenMP
//</AP>

// The number of OpenMP threads used by the code will depend on the setting of the OpenMP environment variable OMP_NUM_THREADS on the system on which the code is run.

// Use this line to define the file output format. Currently only FITS files (cfitsion library) are supported:
#define HAVE_FITS

// For best performance, use all lines above (do not outcomment any) - except for OpenMP for Mode 1 (potential plane generation), where is may cause a (potentially serious) performance degradation.

///////////////////////////////////////////
///////////////////////////////////////////

#include <time.h>
#include "analysis_parameters.h"


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


