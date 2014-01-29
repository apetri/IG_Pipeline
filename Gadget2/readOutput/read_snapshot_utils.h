#ifndef __READ_SNAPSHOT_UTILS
#define __READ_SNAPSHOT_UTILS

struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
};

struct particle_data
{
  float Pos[3];
  float Vel[3];
  float Mass;
  int Type;

  float Rho, U, Temp, Ne;
};

/*Main interface!!!!!*/
int read_snapshot(char *path, char *basename,char *number, char *nfiles, struct io_header_1 *header1,int *NumPart,int *Ngas, struct particle_data **P,int **Id,double *Time,double *Redshift);

/*You can treat these as a black box!*/
int unit_conversion(int NumPart, struct particle_data **Q);
int load_snapshot(char *fname, int files, int *NumPart, int *Ngas, struct io_header_1 *header1, int **IdP, struct particle_data **Q,double *Time, double *Redshift);
int allocate_memory(struct particle_data **P,int **Id,int NumPart);
int reordering(int NumPart,int **IdP,struct particle_data **Q);

#endif