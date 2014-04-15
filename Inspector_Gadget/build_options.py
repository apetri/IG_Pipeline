#Python glue program to generate C options to be read from inifile

import StringIO

############################################
#########Type in your options here!!########
############################################

DEFAULT_OPTIONS_FILENAME = "default_options.ini"
HELP_FILENAME = "parameters_explanations.txt"
C_HANDLER_SOURCE = "analysis_parameters.c"
C_OPTIONS_HEADER = "analysis_parameters.h"

#Scalar options (name = value)

#Format is: (Category,Parameter_name,Parameter_type,Default_value,Help message)

paths=[
("paths","snapshot_path",str,"/scratch/02918/apetri/Storage/sims/snapshots/mQ3-series","path from which Gadget-2 snapshots are read in"),
("paths","plane_output_path",str,"/scratch/02918/apetri/Storage/wl/IG/mQ3-series","path where density and potential planes are written to after TSC insertion"),
("paths","plane_path",str,"scratch/02918/apetri/Storage/wl/IG/mQ3-series","path from which potential planes are read for WL map generation (typically same as Plane_output_path)"),
("paths","map_output_path",str,"/scratch/02918/apetri/Storage/wl/IG/mQ3-series","path where WL maps are written to"),
("paths","galaxy_catalogue_path",str,"/home1/02918/apetri/Surveys/CFHT/13subfields/raytrace_subfields","where galaxy catalogues specifications are (externally supplied)")
]

basenames=[
("basenames","snapshot_name",str,"snapshot","Gadget output prefix"),
("basenames","density_basename",str,"dens","this is prepended to the density planes filename"),
("basenames","potential_basename",str,"pot","this is prepended to the potential planes filename"),
("basenames","galaxy_catalogue_basename",str,"raytrace_subfield","prefix of subfield specification file (excludes subfield number)"),
("basenames","map_basename",str,"WL-","this is prepended to all simulated map filenames"),
("basenames","convergence_basename",str,"conv","this is prepended to convergence map filenames"),
("basenames","shear1_basename",str,"shear1","this is prepended to shear 1 map filenames"),
("basenames","shear2_basename",str,"shear2","this is prepended to shear 2 map filenames"),
("basenames","omega_basename",str,"omega","this is prepended to omega map filenames"),
("basenames","shear_modulus_basename",str,"shear_abs","this is prepended to shear_abs map filenames"),
("basenames","shear_angle_basename",str,"shear_ang","this is prepended to shear_ang map filenames"),
("basenames","deflection1_basename",str,"defl1","this is prepended to deflection angle 1 filenames"),
("basenames","deflection2_basename",str,"defl2","this is prepended to deflection angle 2 filenames"),
("basenames","deflection_total_basename",str,"defl_total","this is prepended to total deflection angle filenames"),
("basenames","deflection_winkel_basename",str,"defl_ang","this is prepended to deflection winkel filenames")
]

comments=[
("comments","plane_comment",str,"(none)","a short comment you can add to planes FITS header"),
("comments","WL_map_comment",str,"(none)","a short comment you can add to maps FITS header")
]

rarely_changed=[
("rarely_changed","feedback",int,-1,"don't touch it!"),
("rarely_changed","nx",int,4096,"Number of grid points on 2D lens planes, 2048 is a good choice for high resolution maps"),
("rarely_changed","nx",int,4096,"Number of grid points on 2D lens planes, 2048 is a good choice for high resolution maps."),
("rarely_changed","NbinsX",int,2048,"Number of bins (pixels in 1 dimension) for weak lensing analysis"),
("rarely_changed","NbinsY",int,2048,"Number of bins (pixels in 1 dimension) for weak lensing analysis"),
]

not_changed=[
("not_changed","species",int,1,"0: gas (SPH), 1: dark matter (halo). So far only 1 used and tested")
]

not_used=[
("nonused","parameter_path",str,".","DON'T TOUCH IT!"),
("nonused","parameterfile1",str,"(none)","DON'T TOUCH IT!"),
("nonused","parameterfile2",str,"(none)","DON'T TOUCH IT!"),
("nonused","modelname",str,"(none)","DON'T TOUCH IT!"),
("nonused","seed",int,1,"DON'T TOUCH IT!")
]

scalar_options=[
("Input/output paths",paths),
("Base names for files",basenames),
("FITS header comments",comments),
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
	if ptype==str:
		return value
	elif ptype in(float,int):
		return str(value)
	else:
		raise ValueError("Unknown parameter type")

def declaration_string(ptype):
	if ptype==str:
		return 'char*'
	elif ptype==int:
		return 'int'
	elif ptype==float:
		return 'float'
	else:
		raise ValueError("Unknown parameter type")

def print_format(ptype):
	if ptype==str:
		return '%s'
	elif ptype==int:
		return '%d'
	elif ptype==float:
		return '%f'
	else:
		raise ValueError("Unknown parameter type")

def declaration_match(ptype):
	if ptype==str:
		return 'strdup(value)'
	elif ptype==int:
		return 'atoi(value)'
	elif ptype==float:
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

	T.write("""\n\nvoid print_options(sys_options *options){
""")
	
	S.write("""\n 
		if(0){
		}""")
	
	#Scalar options
	for set_name,options_set in scalar_options:
		
		T.write("""    printf("\\n[%s]\\n\\n");\n"""%options_set[0][0])

		for set,name,ptype,default,help in options_set:
			S.write(""" else if(MATCH("%s","%s")){
				options->%s = %s;
				}"""%(set,name,name,declaration_match(ptype)))
			T.write("""    printf("%s = %s\\n",options->%s);\n"""%(name,print_format(ptype),name))

	#Vector options
	for set_name,options_set in vector_options:

		T.write("""    printf("\\n[%s]\\n\\n");\n"""%options_set[0][0])

		for set,num_name,arr_name,num_elements,ptype,default in options_set:
			
			S.write(""" else if(MATCH("%s","%s")){
				options->%s = %s;
				}"""%(set,num_name,num_name,declaration_match(ptype)))
			T.write("""    printf("%s = %s\\n",options->%s);\n"""%(num_name,print_format(int),num_name))

			for i in range(num_elements):
				S.write(""" else if(MATCH("%s","%s[%d]")){
				options->%s[%d] = %s;
				}"""%(set,arr_name,i,arr_name,i,declaration_match(ptype)))
				T.write("""    printf("%s[%d] = %s\\n",options->%s[%d]);\n"""%(arr_name,i,print_format(ptype),arr_name,i))
			T.write("""    printf("\\n");\n""")


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

struct analysis_parameters
{
		
""")

	for set_name,options_set in scalar_options:
		for set,name,ptype,default,help in options_set:
			S.write("""%s %s;\n"""%(declaration_string(ptype),name))

	for set_name,options_set in vector_options:
		for set,num_name,arr_name,num_elements,ptype,default in options_set:
			S.write("""int %s;\n"""%num_name)
			S.write("""%s %s[%d];\n"""%(declaration_string(ptype),arr_name,num_elements))

	S.write("""
};

int handler(void *,const char*,const char*,const char*);

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
	file(C_OPTIONS_HEADER,"w").write(generate_options_type(scalar_options,vector_options))
	
	to_print = generate_default_parameter_file(scalar_options,vector_options)

	file(DEFAULT_OPTIONS_FILENAME,"w").write(to_print[0])
	file(HELP_FILENAME,"w").write(to_print[1])