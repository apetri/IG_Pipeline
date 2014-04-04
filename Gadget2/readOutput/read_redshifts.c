#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "read_snapshot_utils.h"

/*This program reads the redshifts of all the snapshots and prints them in a file
A little inefficient, ok, but the best I could do without messing with the snapshots
binary format*/

int main(int argc,char **argv){

	if(argc<3){
		fprintf(stderr,"Usage %s <snapshots_path> <redshift_list_filename>\n",*argv);
		exit(1);
	}

	char *path = argv[1];
	char *basename = "snapshot";
	int i,numSnapshots,filesPerSnapshot;
	char number[16],nfiles[8];
	struct io_header_1 header1;
	int NumPart, Ngas;
	struct particle_data *P;
	int *Id;
	double Time, Redshift;

	//Prompts for user
	fprintf(stderr,"How many snapshots do you want to read in?\n");
	fscanf(stdin,"%d",&numSnapshots);
	fprintf(stderr,"How many files per snapshot?\n");
	fscanf(stdin,"%d",&filesPerSnapshot);

	//Convert int to char* to match read_snapshot signature
	sprintf(nfiles,"%d",filesPerSnapshot);
	fprintf(stderr,"You want to read %d snapshots, %s files per snapshot",numSnapshots,nfiles);

	//Open file for redshift list output
	FILE *out = fopen(argv[2],"w");
	if(out==NULL){
		perror(argv[2]);
		exit(1);
	}



	for(i=1;i<=numSnapshots;i++){

		sprintf(number,"%d",i);

		//Read in all the snapshot
		read_snapshot(path,basename,number,nfiles,&header1,&NumPart,&Ngas,&P,&Id,&Time,&Redshift);

		//Print to file
		fprintf(out,"Snapshot %d --> Time: %.5f Redshift: %.7f\n",i,Time,Redshift);

		//Free snapshot memory
		free(P);

	}

	fprintf(stderr,"\nDone!!\n\n");

	fclose(out);
	return 0;

}