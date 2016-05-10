#include <drfftw_mpi.h>
#include "darkenergy.h"

// <JMK>:
#define PREFAB_GROWTH
// Comment out above line if do not have prefabricated growth factor in parameter file,
// and instead want to have it calculated by Lam Hui's growth factor code which was,
// rigged into this special version of N-GenIC. The code is in FORTRAN 77 though and
// needs C to call FORTRAN (with -lg2c flag for linker in Makefile).
// The Makefile has been adapted for this, but it may not work on all systems.
// Therefore, the above line makes it possible to pre-calculate the growth factor separately.
// N-GenIC also has its own way of calculating the growth factor, but it does not work
// for w(z) dark energy models. While the code is still present in this version and could be
// reactivated, it has been disabled to avoid mistakes.
// </JMK>

#define  PI          3.14159265358979323846 
#define  GRAVITY     6.672e-8
#define  HUBBLE      3.2407789e-18   /* in h/sec */


double PowerSpec(double kmag);
double GrowthFactor(double astart, double aend);
double F_Omega(double a);
int    read_parameter_file(char *fname);
double PowerSpec_EH(double k);
double PowerSpec_Efstathiou(double k);


#ifdef T3E
typedef short int int4byte;	/* Note: int has 8 Bytes on the T3E ! */
typedef unsigned short int uint4byte;	/* Note: int has 8 Bytes on the T3E ! */
#else
typedef int int4byte;
typedef unsigned int uint4byte;
#endif



extern struct io_header_1
{
  uint4byte npart[6];      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
  double mass[6];          /*!< mass[1] gives the particle mass */
  double time;             /*!< time (=cosmological scale factor) of snapshot */
  double redshift;         /*!< redshift of snapshot */
  int4byte flag_sfr;       /*!< flags whether star formation is used (not available in L-Gadget2) */
  int4byte flag_feedback;  /*!< flags whether feedback from star formation is included */
  uint4byte npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
                                the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
  int4byte flag_cooling;   /*!< flags whether radiative cooling is included */
  int4byte num_files;      /*!< determines the number of files that are used for a snapshot */
  double BoxSize;          /*!< Simulation box size (in code units) */
  double Omega0;           /*!< matter density */
  double OmegaLambda;      /*!< vacuum energy density */
  double HubbleParam;      /*!< little 'h' */
  // <JMK>:
  double w0;
  double wa;
  double comoving_distance;
  // </JMK>
  int4byte flag_stellarage;     /*!< flags whether the age of newly formed stars is recorded and saved */
  int4byte flag_metals;         /*!< flags whether metal enrichment is included */

  // <JMK>:
  // int4byte hashtabsize;         /*!< gives the size of the hashtable belonging to this snapshot file */
  int4byte nothing[6]; // new fillers to keep structure the same size across different platforms (not always possible).
  int4byte nothing1;
  int4byte nothing2;
  char fill[32]; // decreased to accomodate additional variables, original was: char fill[84];		/*!< fills to 256 Bytes */
    // </JMK>  (put all custom additional variables at end, such that file still can be read by standard version of Gadget-2 (it will simply miss those additional variables), as header struct is written and read as one 256-byte chunk, where its components appear in order.

}
header, header1;

extern DECosmo de_cosmo;


extern int      Nglass;
extern int      *Local_nx_table;
extern int      WhichSpectrum;


extern FILE     *FdTmp, *FdTmpInput;

extern int      Nmesh, Nsample;

extern int      SphereMode;

extern long long IDStart;


extern char     GlassFile[500]; 
extern char     FileWithInputSpectrum[500];

extern int      TileFac; 

extern double   Box;
extern int Seed;

extern long long TotNumPart;

extern int      NumPart;

extern int      NTaskWithN;


extern int      *Slab_to_task;


extern struct part_data 
{
  float Pos[3];
  float Vel[3];
  long long ID;
} *P;


extern double InitTime;
extern double Redshift;
extern double MassTable[6];


// <JMK>:
// extern char OutputDir[100], FileBase[100];
extern char OutputDir[1000], FileBase[1000];
// </JMK>
extern int  NumFilesWrittenInParallel;


extern int      ThisTask, NTask;

extern int      Local_nx, Local_x_start;

extern int  IdStart;

extern rfftwnd_mpi_plan Inverse_plan;
extern fftw_real        *Disp, *Workspace;
extern fftw_complex     *Cdata;


extern double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
extern double InputSpectrum_UnitLength_in_cm;
extern double G, Hubble;
extern double RhoCrit;

extern double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
extern double OmegaBaryon, HubbleParam;
extern double PrimordialIndex;
extern double ShapeGamma;

extern double Dplus; /* growth factor */


extern int    ReNormalizeInputSpectrum;


// <JMK>:
extern double w0;
extern double wa;

extern double vel_prefac_lam;
double Dplus_Interface(double z2, double z1, double omegamI, double omegaQI, double omegavI, double wQI, double wQpI, double hI);
// </JMK>.
