#Python glue program to generate C options to be read from inifile

import StringIO

############################################
#########Type in your options here!!########
############################################

DEFAULT_OPTIONS_FILENAME = "default_parameters.ini"
HELP_FILENAME = "parameters_explanations.txt"
C_HANDLER_SOURCE = "analysis_parameters.c"
C_OPTIONS_HEADER = "analysis_parameters.h"

#Scalar options (name = value)

#Format is: (Category,Parameter_name,Parameter_type,Default_value,Help message)
mode=[
("mode","mode",int,2,"1: Potential Planes, 2: WL Maps"),
("mode","preload_planes",int,1,"0: read in potential planes as needed for each map (standard operation mode), 1: preload all potential planes at the beginning of the run (into an RMA window) - useful for Blue Gene/Q with NFS PB disk"),
("mode","galaxy_catalogue_type",int,1,"0: Maps with source redshift plane, 1: galaxy catalogue (with galaxies at various redshifts, and angular galaxy positions in radians from map center)")
]

survey_parameters=[
("survey_parameters","source_redshift",float,2.0,"For galaxy catalogue type 0: this is the redshift of the sources"),
("survey_parameters","plane_before_source",int,46,"Ray tracing needs to be performed up to this plane (0-->1-->....-->46)"),
("survey_parameters","survey_angle","float32",3.41,"If set to zero, automatically determined by code (typically set manually, so all cosmologies get maps of same angular size). Do not use this automatic feature")
]

paths=[
("paths","snapshot_path[2000]","char","/scratch/02918/apetri/Storage/sims/snapshots/mQ3-series","path from which Gadget-2 snapshots are read in"),
("paths","plane_output_path[2000]","char","/scratch/02918/apetri/Storage/wl/IG/mQ3-series","path where density and potential planes are written to after TSC insertion"),
("paths","plane_path[2000]","char","scratch/02918/apetri/Storage/wl/IG/mQ3-series","path from which potential planes are read for WL map generation (typically same as Plane_output_path)"),
("paths","planes_folder[200]","char","Planes","Name of folder where planes are stored"),
("paths","map_output_path[2000]","char","/scratch/02918/apetri/Storage/wl/IG/mQ3-series","path where WL maps are written to"),
("paths","maps_folder[200]","char","Maps","Name of folder in which maps are put"),
("paths","galaxy_catalogue_path[2000]","char","/home1/02918/apetri/Surveys/CFHT/13subfields/raytrace_subfields","where galaxy catalogues specifications are (externally supplied)"),
("paths","galaxy_catalogue_output_path[2000]","char","/scratch/02918/apetri/Storage/wl/IG/mQ3-series/galaxy_catalogue_subfields1-13","where simulated catalogues will be saved")
]

basenames=[
("basenames","snapshot_name[200]","char","snapshot","Gadget output prefix"),
("basenames","density_basename[200]","char","dens","this is prepended to the density planes filename"),
("basenames","potential_basename[200]","char","pot","this is prepended to the potential planes filename"),
("basenames","galaxy_catalogue_basename[2000]","char","raytrace_subfield","prefix of subfield specification file (excludes subfield number)"),
("basenames","map_basename[200]","char","WL-","this is prepended to all simulated map filenames"),
("basenames","convergence_basename[200]","char","conv","this is prepended to convergence map filenames"),
("basenames","shear1_basename[200]","char","shear1","this is prepended to shear 1 map filenames"),
("basenames","shear2_basename[200]","char","shear2","this is prepended to shear 2 map filenames"),
("basenames","omega_basename[200]","char","omega","this is prepended to omega map filenames"),
("basenames","shear_modulus_basename[200]","char","shear_abs","this is prepended to shear_abs map filenames"),
("basenames","shear_angle_basename[200]","char","shear_ang","this is prepended to shear_ang map filenames"),
("basenames","deflection1_basename[200]","char","defl1","this is prepended to deflection angle 1 filenames"),
("basenames","deflection2_basename[200]","char","defl2","this is prepended to deflection angle 2 filenames"),
("basenames","deflection_total_basename[200]","char","defl_total","this is prepended to total deflection angle filenames"),
("basenames","deflection_winkel_basename[200]","char","defl_ang","this is prepended to deflection winkel filenames")
]

io_size=[
("i/o_amount","global_first_snapshot",int,0,"First Gadget snapshot to be read"),
("i/o_amount","global_last_snapshot",int,58,"Last Gadget snapshot to be read"),
("i/o_amount","first_realization",int,1,"First map/galaxy/plane realization to be generated"),
("i/o_amount","last_realization",int,256,"In Mode 1: Set to number of lens planes per snapshot (9); in Mode 2, set to number of maps (or foreground realizations per galaxy catalogue subfield), i.e. typically 1000"),
("i/o_amount","first_galaxy_subfield",int,1,"count starts at 1; e.g. the 13 CFHT subfields for our map sizes from the CFHT survey"),
("i/o_amount","last_galaxy_subfield",int,2,"last galaxy subfield: i.e. 13 for CFHT"),
("i/o_amount","seed_block",int,9,"10 is a good value for 1GB memory / CPU on NYBlue/P (i.e. DUAL mode) and nx=2048, 3 for nx=4096 on NYBlue/P and DUAL mode (2 on NYBlue/L and CO mode). 8 for nx=4096 and SMP mode (NYBlue/P only). Use 9 for new minimal plane method (or as many planes as there are per snapshot"),
("i/o_amount","max_realizations",int,1000,"maximum number of realizations allowed during ray-tracing, before random numbers repeat between different subfields"),
("i/o_amount","number_of_plane_realizations",int,9,"This can be later automatized from random plane drawing file"),
("i/o_amount","number_of_sim_ics",int,5,"This can be later automatized from random plane drawing file; relevant in mode 2 only!!!!"),
("i/o_amount","first_sim_ic",int,1,"count starts at 1. First N-body simulation IC used; relevant in mode 2 only!!!!"),
("i/o_amount","snapskip",int,1,"Set to 1 if not snapshots are to be skipped, otherwise to correspondingly higher numbers (2 for j-series of simulations, 1 for the more modern i-series and m-series)"),
("i/o_amount","fiducial",int,1,"For Mode=2:  Submitted cosmologies beyond this number treated as nonfiducial models, cosmologies lower or equal to this number are fiducial (more ICs, fewer plane realizations per IC). Count starts at 1. Set to zero if all are nonfiducial. Leave this at zero. Control fiducial models by giving a a different random number file")
]

randomization=[
("randomization","snapshot_rotation_randomizer_file[2000]","char","/home1/02918/apetri/IG_Pipeline_0.1/Inspector_Gadget/IG_Snapshot_Rotations.txt","file that contains the random snapshot rotations"),
("randomization","plane_randomizer_file_general[2000]","char","/home1/02918/apetri/IG_Pipeline_0.1/Inspector_Gadget/IG_RandomPlanes_general_minimalplane.txt","file that contains the plane randomizations between differen ics; Note: this path is needed for mode=1 (lens plane generation)"),
("randomization","plane_randomizer_file_fiducial[2000]","char","/home1/02918/apetri/IG_Pipeline_0.1/Inspector_Gadget/IG_RandomPlanes_general_minimalplane.txt","set same as above; Note: the fiducial option may not work anymore")
]

comments=[
("comments","plane_comment[200]","char","(none)","a short comment you can add to planes FITS header"),
("comments","WL_map_comment[200]","char","(none)","a short comment you can add to maps FITS header")
]

rarely_changed=[
("rarely_changed","feedback",int,-1,"don't touch it!"),
("rarely_changed","nx",int,4096,"Number of grid points on 2D lens planes, 2048 is a good choice for high resolution maps"),
("rarely_changed","ny",int,4096,"Number of grid points on 2D lens planes, 2048 is a good choice for high resolution maps."),
("rarely_changed","NbinsX",int,2048,"Number of bins (pixels in 1 dimension) for weak lensing analysis"),
("rarely_changed","NbinsY",int,2048,"Number of bins (pixels in 1 dimension) for weak lensing analysis"),
]

not_changed=[
("not_changed","species",int,1,"0: gas (SPH), 1: dark matter (halo). So far only 1 used and tested"),
("not_changed","scramble_mode",int,2,"0: No box rotations or shifts, 1: simple random number (does not work well when code is run parallel), 2: random numbers from prefabricated list (up to 1000 realizations possible). Default is 2."),
("not_changed","plane_shift",int,1,"Shifts planes transversally (with periodic completion) before ray-tracing, causing additional randomization. Depreciated parameter, always set to 0"),
("not_changed","ray_tracing",int,1,"Ray-tracing is on if !=0. Default is 1 (on)"),
("not_changed","cell_embedding",int,2,"Particles on Grid: 1: Cloud-in-cell, 2: TSC. Best is TSC"),
("not_changed","raypoint_averaging",int,3,"1: for linear averaging, 2: for bicubic averaging, 3 (and everything else): for TSC averaging over derivative values on grid"),
("not_changed","convergence_direct",int,0,"If !=0, calculate convergence directly from density planes (for testing mostly, default is 0 for this parameter)"),
("not_changed","plane_padding",float,0,"0 is no padding, i.e. exact size visible by survey, 1.0 is size of survey (extra half on each side), etc. Alway leave this at zero")
]

not_used=[
("nonused","parameter_path[2000]","char",".","DON'T TOUCH IT!"),
("nonused","seed",int,1,"DON'T TOUCH IT!")
]

################################################################################################################
################################################################################################################

#These are NOT read from the ini parameter file, but are filled in automatically by the code, so the 
#analysis_parameter struct needs to have them declared

padding=[
("nxny",int),
("plane_assignment_mode",int),
("plane_averaging_mode",int),
("averaging_mode",int),
("survey_angle_in_rad","float32"),
("snapshots",int),
("number_of_planes",int),
("boxsize","float32"),
("generate_planes",int),
("generate_maps",int),
("first_plane",int),
("last_plane",int),
("parallel",int),
("process_number",int),
("number_of_seeds",int),
("realization",int),
("first_snapshot",int),
("last_snapshot",int),
("number_of_snapshots",int),
("number_of_realizations",int),
("global_first_realization",int),
("global_last_realization",int),
("H_0",float),
("h",float),
("Omega_m",float),
("Omega_Lambda",float),
("w0",float),
("wa",float),
("ns",float),
("sigma_8",float),
("initial_condition",float),
("source_scale_factor",float),
("source_comoving_distance",float),
("ASPP_resolution",float),
("PPAM_resolution",float),
("NumPartTotal[6]",float),
("mass[6]",float),
("plane_urpath[2000]","char"),
("parameterfile1[200]","char"),
("parameterfile2[200]","char"),
("modelname[200]","char"),
("simulation_codename[200]","char"),
("start_time","time"),
("runtime","float32"),
("condor_script_path[1000]","char"),
("condor_script_name[200]","char"),
("condor_script_mode",int),
("storage_node_name[200]","char"),
("storage_Gadget_path[2000]","char"),
("local_Gadget_path[2000]","char"),
("storage_IG_path[2000]","char"),
("local_IG_path[2000]","char"),
("comoving_file[200]","char"),
("simfolder[200]","char"),
("number_of_boxcenters",int),
("number_of_source_planes",int),
("ThisTask",int),
("NTasks",int),
("byteswap",int),
("endianness_set",int),
("diagnostic",int),
("number_of_cosmologies",int),
("cosmology_number",int),
("number_of_processes",int),
("extension[20]","char"),
("darkenergy_initialized",int),
("chi_initialized",int),
("snapshot_allocated",int),
("pass",int),
("galaxy_subfield",int)
]

################################################################################################################
################################################################################################################

scalar_options=[
("IG mode of operation",mode),
("Survey parameters",survey_parameters),
("Input/output paths",paths),
("Base names for files",basenames),
("FITS header comments",comments),
("Input/output amounts",io_size),
("Randomization of planes",randomization),
("Useful parameters, but rarely changed",rarely_changed),
("Unusual parameters, typically not changed",not_changed),
("These parameters are NOT used, don't even touch them!!!",not_used)
]

#Vector options (num_par = Size_of_the_vector; arr_name = Name_of_the_vector)

#Format is (Category,Number_of_parameters_name,array_parameter_name,num_elements,Parameter_type,Default_values) (only int,float supported)

vector_options=[]

###############################################################
#########Do not modify anything below here!!!!#################
###############################################################

def ini_string(value,ptype):
	if ptype==str or ptype=="char":
		return value
	elif ptype in(float,int,"float32"):
		return str(value)
	else:
		raise ValueError("Unknown parameter type")

def declaration_string(ptype):
	if ptype==str:
		return 'char*'
	elif ptype=="char":
		return 'char'
	elif ptype==int:
		return 'int'
	elif ptype==float:
		return 'double'
	elif ptype=="float32":
		return 'float'
	elif ptype=="time":
		return 'time_t'
	else:
		raise ValueError("Unknown parameter type")

def print_format(ptype):
	if ptype==str or ptype=="char":
		return '%s'
	elif ptype==int:
		return '%d'
	elif ptype==float or ptype=="float32":
		return '%f'
	else:
		raise ValueError("Unknown parameter type")

def declaration_match(ptype):
	if ptype==str:
		return 'strdup(value)'
	elif(ptype=="char"):
		return 'value'
	elif ptype==int:
		return 'atoi(value)'
	elif ptype==float or ptype=="float32":
		return 'atof(value)'
	else:
		raise ValueError("Unknown parameter type")

def generate_handler(scalar_options,vector_options):
	S=StringIO.StringIO()
	T=StringIO.StringIO()
	
	S.write("""/* DO NOT EDIT THIS FILE - IT IS GENERATED FROM build_options.py. Edit that. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "analysis_parameters.h"

int handler(void *user,const char *section,const char *name, const char *value){
		
		struct analysis_parameters *options = (struct analysis_parameters *) user;
		
		#define MATCH(s,n) strcmp(section,s)==0 && strcmp(name,n)==0
""")

	T.write("""\n\nvoid print_options(FILE* stream,struct analysis_parameters *options){
""")
	
	S.write("""\n 
		if(0){
		}""")
	
	#Scalar options
	for set_name,options_set in scalar_options:
		
		T.write("""    fprintf(stream,"\\n[%s]\\n\\n");\n"""%options_set[0][0])

		for set,name,ptype,default,help in options_set:
			
			if(ptype=="char"):
				name = name.split("[")[0]

				S.write(""" else if(MATCH("%s","%s")){
					strcpy(options->%s,%s);
					}"""%(set,name,name,declaration_match(ptype)))
			else:
				S.write(""" else if(MATCH("%s","%s")){
					options->%s = %s;
					}"""%(set,name,name,declaration_match(ptype)))


			T.write("""    fprintf(stream,"%s = %s\\n",options->%s);\n"""%(name,print_format(ptype),name))

	#Vector options
	for set_name,options_set in vector_options:

		T.write("""    fprintf(stream,"\\n[%s]\\n\\n");\n"""%options_set[0][0])

		for set,num_name,arr_name,num_elements,ptype,default in options_set:
			
			S.write(""" else if(MATCH("%s","%s")){
				options->%s = %s;
				}"""%(set,num_name,num_name,declaration_match(ptype)))
			T.write("""    fprintf(stream,"%s = %s\\n",options->%s);\n"""%(num_name,print_format(int),num_name))

			for i in range(num_elements):
				S.write(""" else if(MATCH("%s","%s[%d]")){
				options->%s[%d] = %s;
				}"""%(set,arr_name,i,arr_name,i,declaration_match(ptype)))
				T.write("""    fprintf(stream,"%s[%d] = %s\\n",options->%s[%d]);\n"""%(arr_name,i,print_format(ptype),arr_name,i))
			T.write("""    fprintf(stream,"\\n");\n""")


	S.write(""" else{
		return 0;
		}
		return 1;
}""")

	T.write("""}""")

	S.seek(0)
	T.seek(0)
	return S.read(),T.read()

def generate_options_type(scalar_options,vector_options):
	S=StringIO.StringIO()
	S.write("""/* DO NOT EDIT THIS FILE - IT IS GENERATED FROM build_options.py. Edit that. */

#ifndef __ANALYSIS_PARAMETERS_H
#define __ANALYSIS_PARAMETERS_H

#include <time.h>

struct analysis_parameters
{
		
""")

	for set_name,options_set in scalar_options:
		for set,name,ptype,default,help in options_set:
			S.write("""%s %s;\n"""%(declaration_string(ptype),name))

		S.write("""\n""")

	for name,ptype in padding:
		S.write("""%s %s;\n"""%(declaration_string(ptype),name))

	for set_name,options_set in vector_options:
		for set,num_name,arr_name,num_elements,ptype,default in options_set:
			S.write("""int %s;\n"""%num_name)
			S.write("""%s %s[%d];\n"""%(declaration_string(ptype),arr_name,num_elements))

	S.write("""
};

int handler(void *,const char*,const char*,const char*);
void print_options(FILE *stream,struct analysis_parameters *options);

#endif
		
""")

	S.seek(0)
	return S.read()

def generate_default_parameter_file(scalar_options,vector_options):
	S=StringIO.StringIO()
	T=StringIO.StringIO()

	S.write("""# Default options ini file generated from build_options.py #\n""")
	T.write("""IG Parameters guide\n\n\n""")

#scalar options
	for set_name,options_set in scalar_options:
		S.write("""
##########################
##########%s##############
##########################
[%s]

"""%(set_name,options_set[0][0]))
		
		T.write("""
##########################
##########%s##############
##########################

"""%(set_name))

		for set,name,ptype,default,help in options_set:

			if(ptype=="char"):
				name = name.split("[")[0]

			S.write("""%s = %s\n""" %(name,ini_string(default,ptype)))
			T.write("""%s 	: 	%s\n""" %(name,help))

		S.write("""\n""")
		T.write("""\n\n""")

#Vector options
	for set_name,options_set in vector_options:
		S.write("""
##########################
##########%s##############
##########################
[%s]

"""%(set_name,options_set[0][0]))
		for set,num_name,arr_name,num_elements,ptype,default in options_set:
			S.write("""%s = %s\n""" %(num_name,ini_string(num_elements,int)))
			for i in range(num_elements):
				S.write("""%s[%d] = %s\n"""%(arr_name,i,ini_string(default[i],ptype)))
			S.write("""\n""")

		S.write("""\n""")

	
	S.seek(0)
	T.seek(0)
	
	return S.read(),T.read()
	
if(__name__=="__main__"):
	
	to_print = generate_handler(scalar_options,vector_options)
	
	file(C_HANDLER_SOURCE,"w").write(to_print[0])
	file(C_HANDLER_SOURCE,"a").write(to_print[1])
	
	file(C_OPTIONS_HEADER,"w").write(generate_options_type(scalar_options,vector_options))
	
	to_print = generate_default_parameter_file(scalar_options,vector_options)

	file(DEFAULT_OPTIONS_FILENAME,"w").write(to_print[0])
	file(HELP_FILENAME,"w").write(to_print[1])
