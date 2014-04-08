/*
 *  io_routines.c
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 5/5/07.
 *  Copyright 2007. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "main.h"
#include "io_routines.h"


#define MAXLINE 2000
#define MAXRANDOMLINE 10000
////////////////////////////

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



void read_sampler(char filename[], int **readarray, int lines, int columns, int max_realizations, int subfield)
{
  // Note that columns must be two times the number of planes used (because always a pair of numbers is read in (ic, realization) for every plane.

  FILE *input_file;
  char line[MAXRANDOMLINE];
  int line_length;
  char *ptr = NULL;
  int i, j;

    int dummy;
    
  input_file=fopen(filename, "r");
  
    i=0;
    
    // Skip max_realizations*subfield many lines at the beginning of the file:
    for (i=0; i<max_realizations*subfield; i++)
    {
        line_length=fgetline(input_file, line, sizeof(line));
    }
 
  
    i=0;
  while ((line_length=fgetline(input_file, line, sizeof(line))) > 0 && i<lines)
    {
      // printf("Reading a line from input file:\nLine length: %d \nLine Content: %s", line_length, line);
      
        
            j=0;
            // extract numbers from line string:
            ptr = strtok(line, " "); // second argument contains delimiters by which tokens (numbers) are separated.
            while(ptr!=NULL)
            {
                // convert token to integer and insert into array
                if (j<columns) readarray[i][j]=atoi(ptr);
                else dummy=atoi(ptr);
                j++;

                ptr = strtok(NULL, " "); // strtok expects NULL pointer during subsequent calls.
            }
        
        i++;
    }
  
  fclose(input_file);
}



void read_analysis_parameter_file(FILE *input_file, int input_file_type)
{

	char line[MAXLINE];
	int line_length;
	char parameter_name[MAXNAME];
	char path_name[MAXPATHNAME];
	float parameter_value;
	int parameter_assigned, assignment_checksum, path_checksum;
	int snapshots, comoving_number;
	
	int number_of_input_parameters, number_of_paths;
	number_of_input_parameters = 17; // Adjust to the number of parameters above. Enter here number of input parameters to be read from input file. This value is not crucial for performance of the code, but is used with the checksum to confirm the correct number of assigned parameters. If set improperly, a false error message will be displayed, and assignment errors may not be caught.
	number_of_paths=22; // Adjust to the number of paths and filenames (sum of both types) to be read in from parameter file. Needs to be set properly or code will abort for safety reasons.
	// Numbers above need to be updated, otherwise warnings don't work correctly.
	
	path_checksum = 0;
	assignment_checksum = 0; // Used for file_type one to check that all parameters read.
	snapshots = 0;
	comoving_number = 0;     // Used to count through comoving distance read-ins in file_type 2.
	
	if (input_file_type==1)
	{
	while ((line_length=fgetline(input_file, line, sizeof(line))) > 0)
	{
		if (feedback >=2) printf("Reading a line from input file:\nLine length: %d \nLine Content: %s", line_length, line);
		
		//////////// Paths and filenames (strings read in): ////////////////////
		
		if (sscanf(line, "%s %s", parameter_name, path_name) == 2)
		{
			/////////////////////// Paths: ////////////////////////
			
			if (strcmp(parameter_name,"Snapshot_path")==0)
			{
				strncpy(parameters.snapshot_path, path_name, MAXPATHNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Parameter_path")==0)
			{
				strncpy(parameters.parameter_path, path_name, MAXPATHNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Plane_path")==0)
			{
				strncpy(parameters.plane_path, path_name, MAXPATHNAME);
				//parameter_assigned=1;
				path_checksum++;
			}

			if (strcmp(parameter_name,"Plane_output_path")==0)
			{
				strncpy(parameters.plane_output_path, path_name, MAXPATHNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Map_output_path")==0)
			{
				strncpy(parameters.map_output_path, path_name, MAXPATHNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
            
            if (strcmp(parameter_name,"Galaxy_catalogue_path")==0)
            {
                strncpy(parameters.galaxy_catalogue_path, path_name, MAXPATHNAME);
                //parameter_assigned=1;
                path_checksum++;
            }
            
            if (strcmp(parameter_name,"Galaxy_catalogue_output_path")==0)
            {
                strncpy(parameters.galaxy_catalogue_output_path, path_name, MAXPATHNAME);
                //parameter_assigned=1;
                path_checksum++;
            }
            
            if (strcmp(parameter_name,"Galaxy_catalogue_basename")==0)
            {
                strncpy(parameters.galaxy_catalogue_basename, path_name, MAXNAME);
                //parameter_assigned=1;
                path_checksum++;
            }
            
            if (strcmp(parameter_name,"Plane_randomizer_file_general")==0)
            {
                strncpy(parameters.plane_randomizer_file_general, path_name, MAXNAME);
                //parameter_assigned=1;
                path_checksum++;
            }
            
            if (strcmp(parameter_name,"Plane_randomizer_file_fiducial")==0)
            {
                strncpy(parameters.plane_randomizer_file_fiducial, path_name, MAXNAME);
                //parameter_assigned=1;
                path_checksum++;
            }
            
            if (strcmp(parameter_name,"Snapshot_rotation_randomizer_file")==0)
            {
                strncpy(parameters.snapshot_rotation_randomizer_file, path_name, MAXNAME);
                //parameter_assigned=1;
                path_checksum++;
            }
            

			
			//////////////////// Condor Paths: ////////////////////

                        if (strcmp(parameter_name,"Condor_script_path")==0)
                          {
                            strncpy(parameters.condor_script_path, path_name, MAXPATHNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Condor_script_name")==0)
                          {
                            strncpy(parameters.condor_script_name, path_name, MAXNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Storage_node_name")==0)
                          {
                            strncpy(parameters.storage_node_name, path_name, MAXNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Storage_Gadget_path")==0)
                          {
                            strncpy(parameters.storage_Gadget_path, path_name, MAXPATHNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Local_Gadget_path")==0)
                          {
                            strncpy(parameters.local_Gadget_path, path_name, MAXPATHNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Storage_IG_path")==0)
                          {
                            strncpy(parameters.storage_IG_path, path_name, MAXPATHNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Local_IG_path")==0)
                          {
                            strncpy(parameters.local_IG_path, path_name, MAXPATHNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Map_basename")==0)
                          {
                            strncpy(parameters.map_basename, path_name, MAXNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Comoving_file")==0)
                          {
                            strncpy(parameters.comoving_file, path_name, MAXNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Simfolder")==0)
                          {
                            strncpy(parameters.simfolder, path_name, MAXNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Planes_folder")==0)
                          {
                            strncpy(parameters.planes_folder, path_name, MAXNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }

                        if (strcmp(parameter_name,"Maps_folder")==0)
                          {
                            strncpy(parameters.maps_folder, path_name, MAXNAME);
                            //parameter_assigned=1;
                            path_checksum++;
                          }
            
			
			//////////////// File names: //////////////////////////
						
			if (strcmp(parameter_name,"Inspector_parameterfile")==0)
			{
				strncpy(parameters.parameterfile1, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Gadget_parameterfile")==0)
			{
				strncpy(parameters.parameterfile2, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Modelname")==0)
			{
				strncpy(parameters.modelname, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}

			/* In being passed as argument to executable now:
			if (strcmp(parameter_name,"Simulation_codename")==0)
			{
				strncpy(parameters.simulation_codename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			*/
			
			if (strcmp(parameter_name,"Snapshot_name")==0)
			{
				strncpy(parameters.snapshot_name, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Density_basename")==0)
			{
				strncpy(parameters.density_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
		
			if (strcmp(parameter_name,"Potential_basename")==0)
			{
				strncpy(parameters.potential_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			
			///////////// Weak lensing file name components: /////////////////
			if (strcmp(parameter_name,"Convergence_basename")==0)
			{
				strncpy(parameters.convergence_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Shear1_basename")==0)
			{
				strncpy(parameters.shear1_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Shear2_basename")==0)
			{
				strncpy(parameters.shear2_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
                        
			if (strcmp(parameter_name,"Omega_basename")==0)
			  {
			    strncpy(parameters.omega_basename, path_name, MAXNAME);
			    //parameter_assigned=1;
			    path_checksum++;
			  }

			if (strcmp(parameter_name,"Shear_modulus_basename")==0)
			{
				strncpy(parameters.shear_modulus_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			if (strcmp(parameter_name,"Shear_angle_basename")==0)
			{
				strncpy(parameters.shear_angle_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			if (strcmp(parameter_name,"Deflection1_basename")==0)
			{
				strncpy(parameters.deflection1_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Deflection2_basename")==0)
			{
				strncpy(parameters.deflection2_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Deflection_total_basename")==0)
			{
				strncpy(parameters.deflection_total_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Deflection_angle_basename")==0)
			{
				strncpy(parameters.deflection_winkel_basename, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"Plane_comment")==0)
			{
				strncpy(parameters.plane_comment, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
			
			if (strcmp(parameter_name,"WL_map_comment")==0)
			{
				strncpy(parameters.WL_map_comment, path_name, MAXNAME);
				//parameter_assigned=1;
				path_checksum++;
			}
		}
		
		//////////// Parameters: ///////////////////////
		
		if (sscanf(line, "%s %e", parameter_name, &parameter_value) == 2)
		{
			parameter_assigned = 0;
			
						
			if (strcmp(parameter_name,"Feedback")==0)
			{
				parameters.feedback=(int) parameter_value;
				feedback = parameters.feedback;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
			
			
			if (strcmp(parameter_name,"Species")==0)
			{
				parameters.species=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
			
			if (strcmp(parameter_name,"NbinsX")==0)
			{
				parameters.NbinsX=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
			
			
			if (strcmp(parameter_name,"NbinsY")==0)
			{
				parameters.NbinsY=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
			}
			
			if (strcmp(parameter_name,"nx")==0)
			{
				parameters.nx=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
			
			
			if (strcmp(parameter_name,"ny")==0)
			{
				parameters.ny=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
			}

			
			if (strcmp(parameter_name,"Survey_angle")==0)
			{
				parameters.survey_angle=parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
			}
			
			
			if (strcmp(parameter_name,"Plane_padding")==0)
			{
				parameters.plane_padding=parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
			}
			
			if (strcmp(parameter_name,"Plane_assignment_mode")==0)
			{
				parameters.plane_assignment_mode=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
			}
			
			if (strcmp(parameter_name,"Averaging_mode")==0)
			{
				parameters.averaging_mode=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
			}
			
			
			if (strcmp(parameter_name,"Scramble_mode")==0)
			{
				parameters.scramble_mode=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
			
			
			if (strcmp(parameter_name,"Ray_tracing")==0)
			{
				parameters.ray_tracing=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}			
			
			
			if (strcmp(parameter_name,"Generate_planes")==0)
			{
				parameters.generate_planes=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}		


			if (strcmp(parameter_name,"Generate_maps")==0)
			{
				parameters.generate_maps=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}		
			
						
			if (strcmp(parameter_name,"Plane_shift")==0)
			{
				parameters.plane_shift=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}	

			
			if (strcmp(parameter_name,"Seed")==0)
			{
				parameters.seed=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}	
			
			
			if (strcmp(parameter_name,"Cell_embedding")==0)
			{
				parameters.cell_embedding=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}	

			if (strcmp(parameter_name,"Raypoint_averaging")==0)
			{
				parameters.raypoint_averaging=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}	

                        if (strcmp(parameter_name,"Convergence_direct")==0)
			  {
			    parameters.convergence_direct=(int) parameter_value;
			    parameter_assigned=1;
			    assignment_checksum++;

			  }
            
            
			if (strcmp(parameter_name,"Source_redshift")==0)
			{
				parameters.source_redshift=(double) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
            
            // This parameter will need to be phased out.
			if (strcmp(parameter_name,"Plane_before_source")==0)
			{
				parameters.plane_before_source=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
            
            
			if (strcmp(parameter_name,"First_galaxy_subfield")==0)
			{
				parameters.first_galaxy_subfield=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
            
                        if (strcmp(parameter_name,"Last_galaxy_subfield")==0)
			{
			        parameters.last_galaxy_subfield=(int) parameter_value;
			        parameter_assigned=1;
			        assignment_checksum++;

			}


            if (strcmp(parameter_name,"Max_realizations")==0)
			{
				parameters.max_realizations=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
            
            
            if (strcmp(parameter_name,"Preload_planes")==0)
			{
				parameters.preload_planes=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
            
            
            if (strcmp(parameter_name,"Number_of_plane_realizations")==0)
			{
				parameters.number_of_plane_realizations=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
            
            
            if (strcmp(parameter_name,"Number_of_sim_ics")==0)
			{
				parameters.number_of_sim_ics=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
            
            
            if (strcmp(parameter_name,"First_sim_ic")==0)
			{
				parameters.first_sim_ic=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}
            

			//////////////////////////// MPI-version Parameters:  //////////

                        if (strcmp(parameter_name,"Mode")==0)
                          {
                            parameters.mode=(int) parameter_value;
                            parameter_assigned=1;
                            assignment_checksum++;

                          }

			if (strcmp(parameter_name,"Snapskip")==0)
                          {
                            parameters.snapskip=(int) parameter_value;
                            parameter_assigned=1;
                            assignment_checksum++;

                          }
			
			if (strcmp(parameter_name,"Fiducial")==0)
                          {
                            parameters.fiducial=(int) parameter_value;
			    parameter_assigned=1;
			    assignment_checksum++;
                          }


			////////////////////// Cosmological Parameters (to be set to the same values as in the N-body simulation analyzed): /////////////////

			if (strcmp(parameter_name,"H_0")==0)
			{
				parameters.H_0=(double) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}	
			
			if (strcmp(parameter_name,"Omega_m")==0)
			{
				parameters.Omega_m=(double) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}	


			/////////////////////// Additional variables for Multi mode: //////////////////////


			
			if (strcmp(parameter_name,"Number_of_seeds")==0)
			{
				parameters.number_of_seeds=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}	
			
			if (strcmp(parameter_name,"Seed_block")==0)
			{
				parameters.seed_block=(int) parameter_value;
				parameter_assigned=1;
				assignment_checksum++;
				
			}	


			/////////////////////////// Condor Parameters: /////////////////////////


                        if (strcmp(parameter_name,"First_snapshot")==0)
			  {
			    parameters.global_first_snapshot=(int) parameter_value;
			    parameter_assigned=1;
			    assignment_checksum++;

			  }

                        if (strcmp(parameter_name,"Last_snapshot")==0)
			  {
			    parameters.global_last_snapshot=(int) parameter_value;
			    parameter_assigned=1;
			    assignment_checksum++;

			  }

                        if (strcmp(parameter_name,"First_realization")==0)
			  {
			    parameters.first_realization=(int) parameter_value;
			    parameter_assigned=1;
			    assignment_checksum++;

			  }

                        if (strcmp(parameter_name,"Last_realization")==0)
			  {
			    parameters.last_realization=(int) parameter_value;
			    parameter_assigned=1;
			    assignment_checksum++;

			  }
            
                        if (strcmp(parameter_name,"Galaxy_catalogue_type")==0)
              {
			    parameters.galaxy_catalogue_type=(int) parameter_value;
			    parameter_assigned=1;
			    assignment_checksum++;
                
              }

			
			if (parameter_assigned ==1)
				{
					if (feedback>2) printf("Parameter %s read, has been assigned value: %e\n", parameter_name, parameter_value);
				}
				else
				{
					printf("WARNING: Parameter %s is an unknown parameter name. This parameter could not be asigned! Check parameter input file for spelling.\n", parameter_name);
				}
		}
	}
	if (assignment_checksum < number_of_input_parameters) printf("WARNING: Not all required input parameters assigned from input file. This error is probably fatal. Check input file.\n");
	if (assignment_checksum > number_of_input_parameters) printf("WARNING: At least one parameter multiply assigned. The code will probably run correctly, but with the last assigned values from parameter file. Check file for multiple definitions.\n");

	if (path_checksum < number_of_paths)
	{
		printf("Some path or file name was not read in. Code will not run properly. Aborting. Check parameter file for proper path and file name listings. Number of paths/names read in is %d.\n", path_checksum);
		exit(1);
	}
	
	}

	
	/*
	if (input_file_type==2)
	{
		while ((line_length=fgetline(input_file, line, sizeof(line))) > 0)
		{
			if (feedback >=2)
			  {
			    printf("Reading a line from input file:\nLine length: %d \nLine Content: %s", line_length, line);
			    fflush(stdout);
			  }

			if (sscanf(line, "%s %e", parameter_name, &parameter_value) == 2)
			{
				parameter_assigned = 0;
				
				
				if (strcmp(parameter_name,"Comoving_distance")==0)
				{
				        
					if (comoving_number<=parameters.global_last_snapshot+1)
					{
						boxcenter[comoving_number].comoving_distance=parameter_value;
						parameter_assigned=1;
						assignment_checksum++;
					}
					else
					{
						parameter_assigned=2;
					}
					comoving_number++;
					
				}
				
				if (parameter_assigned ==1)
				{
					if (feedback>2) printf("Parameter %s read, has been assigned value: %e\n", parameter_name, parameter_value);
				}
				else if (parameter_assigned ==2)
				{
					if (feedback>2) printf("Parameter %s read with value: %e, but has not been assigned (comoving distance for unused snapshot).\n", parameter_name, parameter_value);
				}
				else
				{
					printf("WARNING: Parameter %s is an unknown parameter name. This parameter could not be assigned! Check parameter input file for spelling.\n", parameter_name);
				}

			}
			
		}
	

			if(assignment_checksum<parameters.global_last_snapshot+1)
			{
			printf("ERROR: Insufficient number of comoving distances specified: more snapshots used than comoving distances exist. Code will not run properly. Aborting");
			fflush(stdout);
			exit(0);
			}
				
	
	}
	*/

	if (assignment_checksum==0) printf("No parameters assigned in this call of read_analysis_parameters,\nread_analysis_parameters probably called with an invalid input_file_type.\n");
	

}

#undef MAXLINE
#undef MAXRANDOMLINE

