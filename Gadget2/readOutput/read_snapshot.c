#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "read_snapshot_utils.h"

 /* here the particle data is at your disposal: this is the only function you have to customize!
  */
int do_what_you_want(int NumPart,struct particle_data *P){
  
  int i;
  P--;

  for(i=1;i<=NumPart;i++){
    printf("%e %e %e\n",P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
  }

  return 0;
}

int main(int argc, char **argv)
{
  
  char path[1024], basename[128], number[16], nfiles[8];

  struct io_header_1 header1;

  int NumPart, Ngas;

  struct particle_data *P;

  int *Id;

  double Time, Redshift;

  //Read path,basename and snapshot number from stdin

  printf("Type the path of the directory with the snapshots:\n");

  if(fgets(path,sizeof(path),stdin)){
    *(strrchr(path,'\n')) = 0;
  } else{
    perror("couldn't read from stdin");
    exit(1);
  }

  printf("Type the basename of you snapshots (e.g. 'snapshot' if your snapshots are named snapshot_xxx):\n");

  if(fgets(basename,sizeof(path),stdin)){
    *(strrchr(basename,'\n')) = 0;
  } else{
    perror("couldn't read from stdin");
    exit(1);
  }

  printf("Type the number of the snapshot you want to read:\n");

  if(fgets(number,sizeof(number),stdin)){
    *(strrchr(number,'\n')) = 0;
  } else{
    perror("couldn't read from stdin");
    exit(1);
  }

  printf("Type the number of files per snapshot (usually it's 1 if you ran Gadget on your local computer):\n");

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