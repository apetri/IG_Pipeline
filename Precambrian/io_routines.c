
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "main.h"
#include "io_routines.h"


#define MAXLINE 5000
#define MAXNAME 1000
#define MAXPATHNAME 5000

// Note: home_path, repository_path, mass_storage_path, and series are global variables. 


int fgetline(FILE *stream, char s[], int lim)
{
  int c, i, j;

  for (i=0;i<lim-1 && (c=getc(stream))!=EOF && c!='\n' && c!='#'; ++i) s[i]=c;
  if (c == '\n' || c == '#')
    {
      s[i] = '\n';
      i++;
    }
  if (c == '#') {for (j=0;j<lim-1-i && (c=getc(stream))!=EOF && c!='\n'; ++j);};
  s[i]='\0';
  return i;
}



void read_parameter_file(char filename[])
{

  char line[MAXLINE];
  int line_length;
  char parameter_name[MAXNAME];
  char path_name[MAXPATHNAME];
  float parameter_value;

  // int number_of_input_parameters, number_of_paths;
  int path_checksum=0;
  int feedback=0; // can turn on feedback here to print lines read in if have difficulties with parameter file.

  FILE *input_file;

  input_file=fopen(filename, "r");
 

      while ((line_length=fgetline(input_file, line, sizeof(line))) > 0)
	{
	  if (feedback >=2) printf("Reading a line from input file:\nLine length: %d \nLine Content: %s", line_length, line);
	  
	  //////////// Paths and filenames (strings read in): ////////////////////
	  
	  if (sscanf(line, "%s %s", parameter_name, path_name) == 2)
	    {
	      /////////////////////// Paths: ////////////////////////
	      
	      if (strcmp(parameter_name,"Home_path")==0)
		{
		  strncpy(home_path, path_name, MAXPATHNAME);
		  //parameter_assigned=1;
		  path_checksum++;
		}

	      if (strcmp(parameter_name,"Repository_path")==0)
		{
                  strncpy(repository_path, path_name, MAXPATHNAME);
                  //parameter_assigned=1;                                      
                  path_checksum++;
		}

	      if (strcmp(parameter_name,"Mass_storage_path")==0)
		{
                  strncpy(mass_storage_path, path_name, MAXPATHNAME);
                  //parameter_assigned=1;
                  path_checksum++;
		}

	      if (strcmp(parameter_name,"Series_name")==0)
		{
                  strncpy(series, path_name, MAXPATHNAME);
                  //parameter_assigned=1;
                  path_checksum++;
		}

	    } // end of large if loop.

	} // end of while loop.

   fclose(input_file);
}


#undef MAXLINE
#undef MAXNAME
#undef MAXPATHNAME

