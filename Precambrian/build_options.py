#Python glue program to generate C options to be read from inifile

import StringIO

############################################
#########Type in your options here!!########
############################################

DEFAULT_OPTIONS_FILENAME = "default_options.ini"
C_HANDLER_SOURCE = "options.c"
C_OPTIONS_HEADER = "options.h"

#Scalar options (name = value)

#Format is: (Category,Parameter_name,Parameter_type,Default_value)

submission=[
("submission","submission_style",int,2),
("submission","machine_endianness",str,"little")
]

paths=[
("paths","home_path",str,"/bgsys/home2/jank"),
("paths","repository_relative_path",str,"Documents/GIT/IG_Pipeline_0.1"),
("paths","mass_storage_path",str,"/bgusr/data01/jank")
]

series=[
("series","series_name",str,"mQ3"),
("series","num_particles_side",int,512),
("series","models_file",str,"cosmologies.txt")
]

settings=[
("settings","power_spectrum_at_zini",int,0),
("settings","flat_universe",int,1),
("settings","remove_old",int,0)
]

power3D=[
("power3D","num_files_snapshot",int,16),
("power3D","num_snapshots",int,60),
("power3D","FFT_grid_size",int,256),
("power3D","number_of_bins",int,256),
("power3D","box_size_snapshot_mpc",int,240),
("power3D","particle_buffer_length",int,10000),
("power3D","model_basename",str,"Om0.260_Ol0.740_w-1.000_ns0.960_si0.798"),
("power3D","realization_number",int,1),
("power3D","power_spectrum_savepath",str,"."),
("power3D","power_spectrum_filebase",str,"power_spectrum_3D"),
("power3D","projection_savepath",str,"."),
("power3D","projection_filebase",str,"projection_2D")
]

scalar_options=[
("Submission type",submission),
("Path names (user specific)",paths),
("Simulation series",series),
("Various settings",settings),
("3D output power spectrum",power3D)
]

#Vector options (num_par = Size_of_the_vector; arr_name = Name_of_the_vector)

#Format is (Category,Number_of_parameters_name,array_parameter_name,num_elements,Parameter_type,Default_values) (only int,float supported)

boxes=[
("box_number","Nboxsize","boxsize",1,float,[240.0]),
("box_number","Nz","z",1,float,[100.0]),
#("box_number","Nseed","seed",10,int,[168757,580133,311652,325145,222701,194340,705031,674951,495306,105884])
("box_number","Nseed","seed",1,int,[168757])
]

vector_options=[
("Number of box sizes,initial redshifts,random seeds to be evaluated",boxes)
]

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
		return 'double'
	else:
		raise ValueError("Unknown parameter type")

def print_format(ptype):
	if ptype==str:
		return '%s'
	elif ptype==int:
		return '%d'
	elif ptype==float:
		return '%lf'
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
#include "options.h"

int handler(void *user,const char *section,const char *name, const char *value){
		
		sys_options *options = (sys_options*) user;
		
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

		for set,name,ptype,default in options_set:
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

#ifndef __OPTIONS_H
#define __OPTIONS_H

typedef struct{
		
""")

	for set_name,options_set in scalar_options:
		for set,name,ptype,default in options_set:
			S.write("""%s %s;\n"""%(declaration_string(ptype),name))

	for set_name,options_set in vector_options:
		for set,num_name,arr_name,num_elements,ptype,default in options_set:
			S.write("""int %s;\n"""%num_name)
			S.write("""%s %s[%d];\n"""%(declaration_string(ptype),arr_name,num_elements))

	S.write("""
} sys_options;

int handler(void *,const char*,const char*,const char*);
void print_options(sys_options *options);

#endif
		
""")

	S.seek(0)
	return S.read()

def generate_default_parameter_file(scalar_options,vector_options):
	S=StringIO.StringIO()
	S.write("""# Default options ini file generated from build_options.py #\n""")

#scalar options
	for set_name,options_set in scalar_options:
		S.write("""
##########################
##########%s##############
##########################
[%s]

"""%(set_name,options_set[0][0]))
		for set,name,ptype,default in options_set:
			S.write("""%s = %s\n""" %(name,ini_string(default,ptype)))

		S.write("""\n""")

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
	return S.read()
	
if(__name__=="__main__"):
	to_print = generate_handler(scalar_options,vector_options)
	file(C_HANDLER_SOURCE,"w").write(to_print[0])
	file(C_HANDLER_SOURCE,"a").write(to_print[1])
	
	file(C_OPTIONS_HEADER,"w").write(generate_options_type(scalar_options,vector_options))
	
	file(DEFAULT_OPTIONS_FILENAME,"w").write(generate_default_parameter_file(scalar_options,vector_options))
