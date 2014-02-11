/*
 *  snapshot_read.c
 *  Gadget-2 Snapshot File Reader (low memory version with buffer)
 *
 *  Created by Jan Michael Kratochvil - on 11/03/2012.
 *  Copyright 2012 Jan Michael Kratochvil. All rights reserved.
 *
 */


#define pi 3.141592653589793
#define Gravi_const 1.0
// Gravitational constant

// Maximal allowable line length for read-in from input file:
#define MAXLINE 1000
// Maximal allowable path name length for paths from input file:
#define MAXPATHNAME 1000
// Maximal allowable name length for parameters and file names in input file:
#define MAXNAME 200
// The above values can be modified as needed for longer line and name lengths.
// WARNING: Make sure you also change the corresponding string array lengths in the struct analysis_parameters parameters in the file main.h!


extern double UnitMass_in_g;

/*
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
  // <JMK change from original sample header:>
  int flag_stellarage;                 //!< flags whether the file contains formation times of star particles 
  int flag_metals;                     //!< flags whether the file contains metallicity values for gas and star particles 
  unsigned int npartTotalHighWord[6];  //!< High word of the total number of particles of each type 
  int  flag_entropy_instead_u;
  // The above ones are unused by public Gadget but still written.
  double   w0;
  double   wa;
  double   comoving_distance;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 - 4*8 - 3*8];  // fills to 256 Bytes 
  // char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  // fills to 256 Bytes 
  // </JMK>
} ;
*/


struct io_header_1
{
	int npart[6];                        /*!< number of particles of each type in this file */
	double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are expli\
										  citly                                                                                                       
										  stored in the mass-block of the snapshot file, otherwise they a\
										  re omitted */
	double time;                         /*!< time of snapshot file */
	double redshift;                     /*!< redshift of snapshot file */
	int flag_sfr;                        /*!< flags whether the simulation was including star formation */
	int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
	unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This c\
										  an be                                                                                                       
										  different from npart if one is dealing with a multi-file snapsh\
										  ot. */
	int flag_cooling;                    /*!< flags whether cooling was included  */
	int num_files;                       /*!< number of files in multi-file snapshot */
	double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
	double Omega0;                       /*!< matter density in units of critical density */
	double OmegaLambda;                  /*!< cosmological constant parameter */
	double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
	// <JMK>:                                                                                                 
	double w0;                           /*!< Dark energy equation of state parameters (not part of standard \
										  Gadget-2 release) */
	double wa;
	double comoving_distance;            /*!< Comoving distance of snapshot from observer at z=0 (a=1) (in h^\
										  -1 kpc) */
	// </JMK>.              
	// Variables below are not used in Inspector Gadget:                                                                                  
	int flag_stellarage;                 /*!< flags whether the file contains formation times of star particl\
										  es */
	int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and \
										  star particles */
	unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
	int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
	// <JMK>:                                                                                                 
	int nothing;                         /*!< Meaningless safety integer for better cross-platform byte align\
										  ment of structure (not defined in C standard) */
	char fill[32];  // original: char fill[60];                  /*!< fills to 256 Bytes */                   
	// </JMK>.                                                                                                
} ;
/*!< holds header for snapshot files */


// extern struct io_header_1 header1;

int open_snapshot_multi(FILE **fd, char *fname, struct io_header_1 *header1);
void close_snapshot_multi(FILE **fd);
int read_particle_positions(FILE **fd, float *particle_position_buffer, int *particles_read_from_file, int *particles_read_this_time, int particles_in_file, int particle_buffer_length);

void check_snapshot_multi(float *buffer, int buffer_length, int particles_read_this_time, double boxsize);
void check_snapshot_complete(int particles_read_from_file, int particles_in_file);

void show_header(struct io_header_1 header1);
void show_particle_position_buffer(float *particle_position_buffer, int number_of_particles);


