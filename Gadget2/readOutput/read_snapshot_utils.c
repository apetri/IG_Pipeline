#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "read_snapshot_utils.h"

/*Do not touch anything below here!*/


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int read_snapshot(char *path, char *basename,char *number, char *nfiles, struct io_header_1 *header1,int *NumPart,int *Ngas, struct particle_data **P,int **Id,double *Time,double *Redshift)
{
  char input_fname[1024];
  int snapshot_number, files;

  snapshot_number = atoi(number);		/* number of snapshot */
  files = atoi(nfiles);			/* number of files per snapshot */


  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  load_snapshot(input_fname, files,NumPart,Ngas,header1,Id,P,Time,Redshift);


  reordering(*NumPart,Id,P);			/* call this routine only if your ID's are set properly */

  unit_conversion(*NumPart,P);		/* optional stuff */

  return 0;
}

/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
int unit_conversion(int NumPart, struct particle_data **Q)
{
  double GRAVITY, BOLTZMANN, PROTONMASS;
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
  double G, Xh, HubbleParam;

  int i;
  double MeanWeight, u, gamma;

  struct particle_data *P = *Q;
  P--;

  /* physical constants in cgs units */
  GRAVITY = 6.672e-8;
  BOLTZMANN = 1.3806e-16;
  PROTONMASS = 1.6726e-24;

  /* internal unit system of the code */
  UnitLength_in_cm = 3.085678e21;	/*  code length unit in cm/h */
  UnitMass_in_g = 1.989e43;	/*  code mass unit in g/h */
  UnitVelocity_in_cm_per_s = 1.0e5;

  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);

  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);


  Xh = 0.76;			/* mass fraction of hydrogen */
  HubbleParam = 0.65;


  for(i = 1; i <= NumPart; i++)
    {
      if(P[i].Type == 0)	/* gas particle */
	{
	  MeanWeight = 4.0 / (3 * Xh + 1 + 4 * Xh * P[i].Ne) * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u = P[i].U * UnitEnergy_in_cgs / UnitMass_in_g;

	  gamma = 5.0 / 3;

	  /* get temperature in Kelvin */

	  P[i].Temp = MeanWeight / BOLTZMANN * (gamma - 1) * u;
	}
    }

    return 0;
}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files, int *NumPart, int *Ngas, struct io_header_1 *header1, int **IdP, struct particle_data **Q,double *Time, double *Redshift)
{
  FILE *fd;
  char buf[200];
  int i, k, dummy, ntot_withmasses;
  int n, pc, pc_new, pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
      if(files > 1)
	sprintf(buf, "%s.%d", fname, i);
      else
	sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  fprintf(stderr,"can't open file `%s`\n", buf);
	  exit(0);
	}

      fprintf(stderr,"reading `%s' ...\n", buf);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(header1, sizeof(struct io_header_1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files == 1)
	{
	  for(k = 0, *NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    *NumPart += header1->npart[k];
	    *Ngas = header1->npart[0];
	}
      else
	{
	  for(k = 0, *NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    *NumPart += header1->npartTotal[k];
	    *Ngas = header1->npartTotal[0];
	}

      for(k = 0, ntot_withmasses = 0; k < 6; k++)
	{
	  if(header1->mass[k] == 0)
	    ntot_withmasses += header1->npart[k];
	}

      if(i == 0)
	
  allocate_memory(Q,IdP,*NumPart);
  
  struct particle_data *P = *Q;
  int *Id = *IdP;

  P--;
  Id--;


      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1->npart[k]; n++)
	    {	
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	      
	    }

	}

      SKIP;

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1->npart[k]; n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;


      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1->npart[k]; n++)
	    {
	      fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses > 0)
	SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1->npart[k]; n++)
	    {
	      P[pc_new].Type = k;

	      if(header1->mass[k] == 0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass = header1->mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses > 0)
	SKIP;


      if(header1->npart[0] > 0)
	{
	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1->npart[0]; n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1->npart[0]; n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1->flag_cooling)
	    {
	      SKIP;
	      for(n = 0, pc_sph = pc; n < header1->npart[0]; n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n = 0, pc_sph = pc; n < header1->npart[0]; n++)
	      {
		P[pc_sph].Ne = 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }


  *Time = header1->time;
  *Redshift = header1->redshift;

  return 0;

}




/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(struct particle_data **P,int **Id,int NumPart)
{
  fprintf(stderr,"allocating memory...\n");

  if(!(*P = malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  //P--;				/* start with offset 1 */


  if(!(*Id = malloc(NumPart * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  //Id--;				/* start with offset 1 */

  fprintf(stderr,"allocating memory...done\n");

  return 0;
}




/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(int NumPart,int **IdP,struct particle_data **Q)
{
  int i;
  int idsource, idsave, dest;
  struct particle_data psave, psource;

  struct particle_data *P = *Q;
  int *Id = *IdP;

  P--;
  Id--;


  fprintf(stderr,"reordering....\n");

  for(i = 1; i <= NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource = P[i];
	  idsource = Id[i];
	  dest = Id[i];

	  do
	    {
	      psave = P[dest];
	      idsave = Id[dest];

	      P[dest] = psource;
	      Id[dest] = idsource;

	      if(dest == i)
		break;

	      psource = psave;
	      idsource = idsave;

	      dest = idsource;
	    }
	  while(1);
	}
    }

  fprintf(stderr,"done.\n");

  Id++;
  free(Id);

  fprintf(stderr,"space for particle ID freed\n");

  return 0;
}
