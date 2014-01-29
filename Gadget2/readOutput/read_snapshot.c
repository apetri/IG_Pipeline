#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "read_snapshot_utils.h"

/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */

int main(int argc, char **argv)
{
  
  char path[1024], basename[128], number[16], nfiles[8];

  struct io_header_1 header1;

  int NumPart, Ngas;

  struct particle_data *P;

  int *Id;

  double Time, Redshift;

  //Read path,basename and snapshot number from stdin

  if(fgets(path,sizeof(path),stdin)){
    *(strrchr(path,'\n')) = 0;
  } else{
    perror("couldn't read from stdin");
    exit(1);
  }

  if(fgets(basename,sizeof(path),stdin)){
    *(strrchr(basename,'\n')) = 0;
  } else{
    perror("couldn't read from stdin");
    exit(1);
  }

  if(fgets(number,sizeof(number),stdin)){
    *(strrchr(number,'\n')) = 0;
  } else{
    perror("couldn't read from stdin");
    exit(1);
  }

  if(fgets(nfiles,sizeof(nfiles),stdin)){
    *(strrchr(nfiles,'\n')) = 0;
  } else{
    perror("couldn't read from stdin");
    exit(1);
  }

  read_snapshot(path,basename,number,nfiles,&header1,&NumPart,&Ngas,&P,&Id,&Time,&Redshift);

  do_what_you_want(NumPart,P);

  //free memory
  free(P);

  return 0;
}