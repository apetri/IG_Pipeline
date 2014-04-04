/*                                                                                                                                                                                       
 *  endianness.c                                                                                                                                                                         
 *  Inspector Gadget                                                                                                                                                                       
 *                                                                                                                                                                                       
 *  Created by Jan Michael Kratochvil at Columbia University on 9/19/09.                                                                                                                                                
 *  Copyright 2009. All rights reserved.                                                                                                                               
 *                                                                                                                                                                                       
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <time.h>

#include "main.h"
#include "endianness.h"

#define LITTLE_END 0
#define BIG_END 1


struct endian_parameters endianness;


// Determines endianness of machine code is running on:                                                                                                                                  
int machineEndianness(){
  long int i = 1;
  const char *p = (const char *) &i;
  if (p[0] == 1)  // Lowest address contains the least significant byte                                                                                                                 
    {
      printf("Machine is little-endian.\n");
      fflush(stdout);
      return LITTLE_END;
    }
  else
    {
      printf("Machine is big-endian.\n");
      fflush(stdout);
      return BIG_END;
    }
}





// Byte reversal for standard types (short, int, float, double):                                                                                                                         

inline short ShortSwap( short s )
{
  unsigned char b1, b2;

  b1 = s & 255;
  b2 = (s >> 8) & 255;

  return (b1 << 8) + b2;
}


inline int IntSwap (int i)
{
  unsigned char b1, b2, b3, b4;

  b1 = i & 255;
  b2 = ( i >> 8 ) & 255;
  b3 = ( i>>16 ) & 255;
  b4 = ( i>>24 ) & 255;

  return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
}

inline long long LongSwap (long long i)
{
	printf("ERROR: LongSwap not implemented yet.\n");
	exit(1);
}

inline float FloatSwap( float f )
{
  union
  {
    float f;
    unsigned char b[4];
  } dat1, dat2;

  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.f;
}


inline double DoubleSwap( double f )
{
  union
  {
    double f;
    unsigned char b[8];
  } dat1, dat2;

  dat1.f = f;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.f;
}


// fread with automatic endianness correction for various types:                                                                                                                         
// (call single version only if one variable is read in (and even in that case the normal version can be called)).                                                                       

inline void singleint_fread(int *x, FILE *fd)
{
  fread(x, sizeof(int), 1, fd);
  if (endianness.byteswap!=0) *x=IntSwap(*x); // swap bytes in int if endianness of machine and snapshot file are different.                                              
}


inline void int_fread(int *x, size_t size, size_t number, FILE *fd)
{
  int i;
  fread(x, size, number, fd);
  if (endianness.byteswap!=0) for (i=0; i<number; i++) x[i]=IntSwap(x[i]); // swap bytes in int if endianness of machine and snapshot file are different.                 
}

inline void singlelong_fread(long long *x, FILE *fd)
{
  fread(x, sizeof(long long), 1, fd);
  if (endianness.byteswap!=0) *x=LongSwap(*x); // swap bytes in int if endianness of machine and snapshot file are different.                                              
}

inline void long_fread(long long *x, size_t size, size_t number, FILE *fd)
{
  int i;
  fread(x, size, number, fd);
  if (endianness.byteswap!=0) for (i=0; i<number; i++) x[i]=LongSwap(x[i]); // swap bytes in int if endianness of machine and snapshot file are different.                 
}

inline void singlefloat_fread(float *x, FILE *fd)
{
  fread(x, sizeof(float), 1, fd);
  if (endianness.byteswap!=0) *x=FloatSwap(*x); // swap bytes in float if endianness of machine and snapshot file are different.                                          
}


inline void float_fread(float *x, size_t size, size_t number, FILE *fd)
{
  int i;
  fread(x, size, number, fd);
  if (endianness.byteswap!=0) for (i=0; i<number; i++) x[i]=FloatSwap(x[i]); // swap bytes in float if endianness of machine and snapshot file are different.             
}


inline void singledouble_fread(double *x, FILE *fd)
{
  fread(x, sizeof(double), 1, fd);
  if (endianness.byteswap!=0) *x=DoubleSwap(*x); // swap bytes in double if endianness of machine and snapshot file are different.                                        
}


inline void double_fread(double *x, size_t size, size_t number, FILE *fd)
{
  int i;
  fread(x, size, number, fd);
  if (endianness.byteswap!=0) for (i=0; i<number; i++) x[i]=DoubleSwap(x[i]); // swap bytes in double if endianness of machine and snapshot file are different.           
}



inline void header_fread(struct io_header_1 *x, size_t size, size_t number, FILE *fd)
{
  if (size!=sizeof(struct io_header_1)) // extra precaution for function to be called correctly:                                                                                   
    {
      printf("ERROR: header_fread used with non-header1 data type.\n");
      exit(1002);
    }

  fread(x, size, number, fd);

  if (parameters.byteswap!=0) // swap bytes in header structure if endianness of machine and snapshot file are different:                                                 
    {

      int i;

      for (i=0; i<6; i++)
	{
	  x->npart[i]=IntSwap(x->npart[i]);
	  x->mass[i]=DoubleSwap(x->mass[i]);
	  x->npartTotal[i]=IntSwap(x->npartTotal[i]);
	}

      x->time=DoubleSwap(x->time);
      x->redshift=DoubleSwap(x->redshift);
      x->flag_sfr=IntSwap(x->flag_sfr);
      x->flag_feedback=IntSwap(x->flag_feedback);
      x->flag_cooling=IntSwap(x->flag_cooling);
      x->num_files=IntSwap(x->num_files);
      x->BoxSize=DoubleSwap(x->BoxSize);
      x->Omega0=DoubleSwap(x->Omega0);
      x->OmegaLambda=DoubleSwap(x->OmegaLambda);
      x->HubbleParam=DoubleSwap(x->HubbleParam);
	  x->w0=DoubleSwap(x->w0);
	  x->wa=DoubleSwap(x->wa);
	  x->comoving_distance=DoubleSwap(x->comoving_distance);

    }

}



// Determines endianness of Gadget-2 snapshot loaded relative to machine endianness (returns 0 if same endianness as machine, otherwise 1):
int snapshotEndianness(char *fname, int file_number)
{
  FILE *fd;
  char   buf[200];
  int    dummy;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);


  sprintf(buf,"%s.%d",fname,file_number);

  printf("Buffer:\n%s\n", buf);
  fflush(stdout);

  if(!(fd=fopen(buf,"r")))
    {
      printf("can't open file `%s`\n",buf);
      fflush(stdout);
#if defined(MPI_COMPILE)
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
      exit(1);
    }

  printf("Checking Endianness of Snapshot `%s' ...\n",buf); fflush(stdout);

  fread(&dummy, sizeof(dummy), 1, fd);
  printf("First integer in snapshot file (used for endianness determination): %d\n", dummy);
  fflush(stdout);
  fread(&header1, sizeof(header1), 1, fd);
  fclose(fd);
  if (dummy==256)
    {
      printf("Snapshot file has same endianness as machine. Determined based on first integer in file.\n");
      return 0;                                                                                                                                     
    }
  else if (dummy==65536)
    {
      printf("Snapshot file has opposite endianness than machine. Byte swap will be necessary upon read-in. Determined based on first integer in file.\n");
      return 1;                                                                                                                                           
    }
  else
    {
      printf("First endianness check ambiguous. Performing secondary check based on header content...\n");
      if (header1.time>=0 && header1.time<1.1 && header1.redshift>-0.1 && header1.redshift<10000.0 && header1.Omega0>=0.0 && header1.Omega0<=1.0) // relax conditions if studying some exotic cosmology or parameter range.                                                                                                                                            
	{
	  printf("Snapshot file has same endianness as machine;\n determined based on scale factor, redshift and Omega_matter in header: %e %e %e.\n", header1.time, header1.redshift, header1.Omega0);
	  return 0;                                                                                                                             
	}
      else
	{
	  printf("Snapshot file has opposite endianness than machine, byte swap will be necessary upon read-in;\n determined based on scale factor, redshift and Omega_matter in header: %e %e %e.\n", header1.time, header1.redshift, header1.Omega0);
	  return 1;                                                                                                                                   
	}

    }
}


