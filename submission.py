import sys,os,shutil,stat
import ConfigParser
import StringIO

################################################################
#This function returns the BGQ corner ID given the corner index#
################################################################
def corner(index):
	if(index==1):
		return "R00-M0-N00-J00"
	elif(index==2):
		return "R00-M0-N00-J29"
	elif(index==3):
		return "R00-M0-N01-J01"
	elif(index==4):
		return "R00-M0-N01-J28"
	elif(index==5):
		return "R00-M0-N02-J12"
	elif(index==6):
		return "R00-M0-N02-J17"
	elif(index==7):
		return "R00-M0-N03-J13"
	elif(index==8):
		return "R00-M0-N03-J16"
	else:
		raise ValueError("The corner index must be less than or equal to 8!")

####################################################################
#This function returns the number of needed sub blocks given the####
#number of cosmological models, the number of simulations per model#
#and the maximum number of simulations a block can run##############
####################################################################
def needed(num_simulations,max_in_block):
	return num_simulations/max_in_block + cmp(num_simulations%max_in_block,0)

##########################################################################################################
#This function checks how many simulations should be run in total, and prints according prompts to screen#
##########################################################################################################
def num_simulations_check(options):
	
	cosmologies = options.options("cosmologies_gadget")
	sims_per_model = options.getint("series","simulations_per_model")
	max_sims_sub_block = options.getint("topology","max_sims_sub_block")

	print "You want to run %d models, %d simulations per model, for a total of %d simulations"%(len(cosmologies),sims_per_model,len(cosmologies)*sims_per_model)
	
	needed_blocks = needed(len(cosmologies)*sims_per_model,max_sims_sub_block)

	num_sub_blocks = options.getint("topology","num_sub_blocks")
	if(needed_blocks>num_sub_blocks):
		
		print "You need more than %d sub-blocks! You will need to split up your work!"%num_sub_blocks
		return needed_blocks
	
	else:
		
		Nside = options.getint("series","num_particles_side")
		print "Each sub-block can run %d simulations with %dx%dx%d particles"%(max_sims_sub_block,Nside,Nside,Nside)
		print "You will need %d sub-blocks for this batch"%needed_blocks
		print "You will need to select the sub-blocks you want to use among these (make sure they are free)"
		for i in range(1,num_sub_blocks+1):
			print "id %d --> %s"%(i,corner(i))
		return needed_blocks  

############################################################
#This function generates the CAMB submission script for BGQ#
############################################################
def generate_BGQ_camb_submission(options):
	
	S = StringIO.StringIO()
	
	S.write("""#!/bin/sh

#Do not edit!! This script is generated automatically by submission.py

""")

	#Write blockid to use
	S.write("""BLOCKID=%s\n\n"""%options.get("topology","blockid"))

	#Write relevant paths
	S.write("""HOMEPATH=%s
REPOSITORY_DIR=%s
WORKING_DIR=$HOMEPATH/$REPOSITORY_DIR/localStorage/ics/%s-series/data_CAMB/Output_Data

EXECUTABLE=$HOMEPATH/$REPOSITORY_DIR/camb/camb
LOGFILE_ROOT=log_camb_%s-series_OMP16
"""%(options.get("paths","home_path"),options.get("paths","repository_path"),options.get("series","series_name"),options.get("series","series_name")))

	#Check for how many cosmological models to calculate power spectra
	cosmologies = options.options("cosmologies_camb")
	#Adjust nodes allocation accordingly 
	S.write("""\nCORES_PER_NODE=%d
NUM_MPI_TASKS=%d
"""%(options.getint("computing_resources","cores_per_node_camb"),len(cosmologies)))

	#Select corner
	S.write("""
CORNER=%s
SHAPE=%s
"""%(corner(options.getint("topology","camb_corner")),options.get("topology","corner_shape")))

	#Set Parameters and Logs directories
	S.write("""
LOGSDIR=$HOMEPATH/$REPOSITORY_DIR/localStorage/ics/mQ3-series/data_CAMB/Logs
PARAMETER_DIR=$HOMEPATH/$REPOSITORY_DIR/localStorage/ics/mQ3-series/data_CAMB/Parameters

""")

	#Set parameter filenames
	i=1
	all_parameters = ""
	
	for cosmology_id in cosmologies:
		S.write("""FILE_%d=params_%s-%s.ini\n"""%(i,options.get("series","series_name"),options.get("cosmologies_camb",cosmology_id)))
		all_parameters = all_parameters + "$PARAMETER_DIR/$FILE_%d "%i
		i = i+1

	S.write("""
ARGS="%s"
"""%all_parameters)

	#Write the execution part of the script
	S.write("""
#execution

echo "You will be executing this command:"
echo ""
echo "runjob --block $BLOCKID --corner $CORNER --shape $SHAPE --exe $EXECUTABLE -p $CORES_PER_NODE --np $NUM_MPI_TASKS --args $ARGS --cwd $WORKING_DIR --envs OMP_NUM_THREADS=16 > $LOGSDIR/$LOGFILE_ROOT.out 2> $LOGSDIR/$LOGFILE_ROOT.err"
echo ""
echo "Do you wish to proceed? (y/n)"

read ANSWER
# ANSWER="y"

if [ $ANSWER == "y" ]; then
	runjob --block $BLOCKID --corner $CORNER --shape $SHAPE --exe $EXECUTABLE -p $CORES_PER_NODE --np $NUM_MPI_TASKS --args $ARGS --cwd $WORKING_DIR --envs OMP_NUM_THREADS=16 > $LOGSDIR/$LOGFILE_ROOT.out 2> $LOGSDIR/$LOGFILE_ROOT.err &
else
	echo "Aborting"
	exit 1
fi
""")

	S.seek(0)
	return S.read()


###############################################################
#This function generates the N-GenIC submission script for BGQ#
###############################################################
def generate_BGQ_ngenic_submission(options):

	S = StringIO.StringIO()

	#Set blockid
	S.write("""#!/bin/sh

#Do not edit!! This script is generated automatically by submission.py		

BLOCKID=%s
CORNER=%s 
SHAPE=%s
"""%(options.get("topology","blockid"),corner(options.getint("topology","camb_corner")),options.get("topology","corner_shape")))

	#Set paths,executable name and log file basenames and directories
	S.write("""
HOMEPATH=%s
REPOSITORY_DIR=%s

EXECUTABLE=$HOMEPATH/$REPOSITORY_DIR/N-GenIC/N-GenICq
LOGFILE_ROOT=log_N-GenIC_%s-series
LOGSDIR=$HOMEPATH/$REPOSITORY_DIR/localStorage/ics/mQ3-series/data_N-GenIC/Logs
"""%(options.get("paths","home_path"),options.get("paths","repository_path"),options.get("series","series_name")))

	#Set path for N-GenIC parameter files
	S.write("""
IC_ROOT=$HOMEPATH/$REPOSITORY_DIR/localStorage/ics/%s-series/data_N-GenIC/Parameters
"""%options.get("series","series_name"))

	#Set computing resources usage
	S.write("""

CORES_PER_NODE=%d
NUM_MPI_TASKS=%d
"""%(options.getint("computing_resources","cores_per_node_ngenic"),options.getint("computing_resources","tasks_per_simulation_ngenic")))

	#Create a list with all the Gadget parameter file names (one for each simulation to run)
	parameter_filenames = []
	series_name = options.get("series","series_name")
	cosmologies = options.options("cosmologies_gadget")
	num_particles_side = options.getint("series","num_particles_side")
	box_size_kpc = options.getint("series","box_size_kpc")
	simulations_per_model = options.getint("series","simulations_per_model")
	tasks_per_sim = options.getint("computing_resources","tasks_per_simulation_ngenic")
	
	for cosmology_id in cosmologies:
		for i in range(1,simulations_per_model+1):
			parameter_filenames.append("ics_%s-%db%d_%s_ic%d.param"%(series_name,num_particles_side,box_size_kpc,options.get("cosmologies_gadget",cosmology_id),i))

	#parameter_filenames now contains all the names of the N-GenIC parameter files, one for each simulation

	#Write the execution part of the script
	ngenic_arguments = ""
	for filename in parameter_filenames:
		ngenic_arguments = ngenic_arguments + "%s "%filename

	runjob_command = "runjob --block $BLOCKID --corner $CORNER --shape $SHAPE --exe $EXECUTABLE -p $CORES_PER_NODE --np $NUM_MPI_TASKS"
	runjob_command = runjob_command + " --args %s --cwd $IC_ROOT > $LOGSDIR/$LOGFILE_ROOT.out 2> $LOGSDIR/$LOGFILE_ROOT.err"%ngenic_arguments

	S.write("""

echo "You will be executing this command:"
echo ""
echo "%s"
echo ""
echo "Do you wish to proceed? (y/n)"

read ANSWER

if [ $ANSWER == "y" ]; then

        %s &

else
	echo "Aborting"
	exit 1
fi
"""%(runjob_command,runjob_command))

	S.seek(0)
	return S.read()


##############################################################
#This function generates the Gadget submission script for BGQ#
##############################################################
def generate_BGQ_Gadget_submission(options,used_blocks,first,last):

	S = StringIO.StringIO()

	#Set blockid
	S.write("""#!/bin/sh

#Do not edit!! This script is generated automatically by submission.py		

BLOCKID=%s
"""%options.get("topology","blockid"))

	#Set paths,executable name and log file basenames and directories
	S.write("""
HOMEPATH=%s
REPOSITORY_DIR=%s

EXECUTABLE=$HOMEPATH/$REPOSITORY_DIR/Gadget2/Gadget2q_OMP2_G800_TOPNODE16
LOGFILE_ROOT=log_Gadget_%s-series
LOGSDIR=$HOMEPATH/$REPOSITORY_DIR/localStorage/ics/mQ3-series/data_Gadget/Logs
"""%(options.get("paths","home_path"),options.get("paths","repository_path"),options.get("series","series_name")))

	#Set path for Gadget parameter files
	S.write("""
G_ROOT=$HOMEPATH/$REPOSITORY_DIR/localStorage/ics/%s-series/data_Gadget/Parameters
"""%options.get("series","series_name"))

	#Set computing resources usage
	S.write("""

CORES_PER_NODE=%d
"""%(options.getint("computing_resources","cores_per_node_gadget")))

	#Create a list with all the Gadget parameter file names (one for each simulation to run)
	parameter_filenames = []
	series_name = options.get("series","series_name")
	mass_storage_path = "%s/Storage/sims/snapshots/%s-series"%(options.get("paths","mass_storage_path"),series_name)
	cosmologies = options.options("cosmologies_gadget")
	num_particles_side = options.getint("series","num_particles_side")
	box_size_kpc = options.getint("series","box_size_kpc")
	simulations_per_model = options.getint("series","simulations_per_model")
	tasks_per_sim = options.getint("computing_resources","tasks_per_simulation_gadget")
	
	for cosmology_id in cosmologies:
		for i in range(1,simulations_per_model+1):
			filename_root = "%s-%db%d_%s_ic%d"%(series_name,num_particles_side,box_size_kpc,options.get("cosmologies_gadget",cosmology_id),i)
			parameter_filenames.append(filename_root+".param")
			#Create corresponding snapshot directory on mass storage disk
			try:
				os.mkdir("%s/%s"%(mass_storage_path,filename_root))
				print "Created %s/%s"%(mass_storage_path,filename_root)
			except OSError:
				print "%s/%s already exists!!"%(mass_storage_path,filename_root)

	#parameter_filenames now contains all the names of the Gadget parameter files, one for each simulation
	parameter_filenames = parameter_filenames[first-1:last]

	#print list of simulations that will run
	print "\nIn this sub-batch you will run the following simulations:"
	for filename in parameter_filenames:
		print filename
	print ""

	#Copy the outputs_xxx-series.txt into the Gadget parameters directory
	try:
		outputs_filename = "%s/%s/outputs_%s-series.txt"%(options.get("paths","home_path"),options.get("paths","repository_path"),series_name)
		target_filename = "%s/%s/localStorage/ics/%s-series/data_Gadget/Parameters/outputs_%s-series.txt"%(options.get("paths","home_path"),options.get("paths","repository_path"),series_name,series_name)
		shutil.copy(outputs_filename,target_filename)
		print "Copied Gadget outputs file in %s"%target_filename
	except IOError:
		print "%s doesn not exist! Create it and rerun submission.py!!"%outputs_filename
		exit(1)
	
	#Check for correct number of blocks
	max_sims_sub_block = options.getint("topology","max_sims_sub_block")
	if(len(used_blocks)!=needed(len(parameter_filenames),max_sims_sub_block)):
		raise ValueError("Provided number of blocks doesn't match the one required by the number of simulations to run!")
	
	#Now divide work between sub-blocks, each one runs an instance of runjob, with the appropriate number of simulations
	for sub_block in used_blocks:
		gadget_arguments = ""
		for num_sims in range(1,max_sims_sub_block+1):
			try:
				filename = parameter_filenames.pop(0)
				#Add filename as argument for Gadget
				gadget_arguments = gadget_arguments+"%s "%filename
			except IndexError:
				#Extremely bad programming practice, sorry :(
				num_sims = num_sims - 1 
				break

		gadget_arguments = "%d %d "%(num_sims,tasks_per_sim) + gadget_arguments
		total_mpi_tasks = num_sims*tasks_per_sim

		runjob_command = "runjob --block $BLOCKID --corner %s --shape %s --exe $EXECUTABLE -p $CORES_PER_NODE --np %d"%(corner(sub_block),options.get("topology","corner_shape"),total_mpi_tasks) 
		runjob_command = runjob_command + " --args %s --cwd $G_ROOT"%(gadget_arguments)
		runjob_command = runjob_command + " > $LOGSDIR/${LOGFILE_ROOT}_Bl%d.out 2> $LOGSDIR/${LOGFILE_ROOT}_Bl%d.err"%(sub_block,sub_block)

		#Write the execution part of the script
		S.write("""


echo "You will be executing this command:"
echo ""
echo "%s"
echo ""
echo "Do you wish to proceed? (y/n)"

read ANSWER

if [ $ANSWER == "y" ]; then
	%s &
else
	echo "Aborting"
	exit 1 
fi
"""%(runjob_command,runjob_command)) 

	S.seek(0)
	return S.read()

#################################################################
#Here we write the submission scripts to the appropriate folders#
#################################################################
if(__name__=="__main__"):

	#Check if options file is provided
	if(len(sys.argv)<3):
		print "\nUsage: python %s <ini_options_file> <mode>"%sys.argv[0]
		print "Operation modes:\n"
		print "1: Generate CAMB submission script for Blue Gene Q"
		print "2: Generate N-GenIC submission script for Blue Gene Q"
		print "3: Generate Gadget2 submission script for Blue Gene Q"
		print "\n"
		exit(1)

	#Parse options from ini file
	options = ConfigParser.RawConfigParser()
	options.readfp(file(sys.argv[1],"r"))
	
	if(sys.argv[2]=="1"):
		#CAMB submission script
		print "Generating CAMB submission script..."

		camb_script_directory = "%s/%s/localStorage/ics/%s-series/data_CAMB/Jobs/"%(options.get("paths","home_path"),options.get("paths","repository_path"),options.get("series","series_name"))
		camb_script_filename = "jobsubmitQ_CAMB_%s-series.sh"%options.get("series","series_name")
		file(camb_script_directory+camb_script_filename,"w").write(generate_BGQ_camb_submission(options))
		print "\n%s job script written!\n"%(camb_script_directory+camb_script_filename)
		os.chmod(camb_script_directory+camb_script_filename,stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)

	elif(sys.argv[2]=="2"):
		#N-GenIC submission script
		print "Generating N-GenIC submission script..."

		ngenic_script_directory = "%s/%s/localStorage/ics/%s-series/data_N-GenIC/Jobs/"%(options.get("paths","home_path"),options.get("paths","repository_path"),options.get("series","series_name"))
		ngenic_script_filename = "jobsubmitQ_N-GenIC_%s-series.sh"%options.get("series","series_name")
		file(ngenic_script_directory+ngenic_script_filename,"w").write(generate_BGQ_ngenic_submission(options))
		print "\n%s job script written!\n"%(ngenic_script_directory+ngenic_script_filename)
		os.chmod(ngenic_script_directory+ngenic_script_filename,stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)

	elif(sys.argv[2]=="3"):
		#Gadget2 submission script
		print "Generating Gadget2 submission script...\n"
		
		#Check that the partition can handle your work
		total_simulations = len(options.options("cosmologies_gadget"))*options.getint("series","simulations_per_model")
		needed_blocks = num_simulations_check(options)
		num_sub_blocks = options.getint("topology","num_sub_blocks")
		
		if(needed_blocks>num_sub_blocks):
			
			print "Your job is too large: I'm selecting for you all the blocks"
			print "You will still need to split your job"
			
			used_blocks = range(1,num_sub_blocks+1)

		else:
			#If there are enough sub blocks, proceed to select the blocks you want to run on
			print "\n"
			used_blocks = []
			
			for i in range(needed_blocks):
				
				print "Please select the id of block %d"%(i+1)
				selected_block = int(raw_input("-->"))
				while(True):
					if(selected_block in used_blocks):
						print "sub-block with id %d already selected! Choose another one!"%selected_block
						print "Selected blocks so far"
						print used_blocks
						print "Please select the id of block %d"%(i+1)
						selected_block = int(raw_input("-->"))
					else:
						used_blocks.append(selected_block)
						break

		print "\nThese are the sub-blocks you selected:"
		for i in used_blocks:
			print "%d --> %s"%(i,corner(i))
		print ""

		if(needed_blocks>num_sub_blocks):
			
			#Need to split the job in this case!!
			answer = "y"

		else:

			#Prompt user if he wants to further split the job
			print "Do you want to further split this batch in several sub-batches?(y/n)"
			answer = raw_input("-->")

		if(answer=="n"):

			#Now the sub-blocks are selected, we can proceed in writing the submission script
			gadget_script_directory = "%s/%s/localStorage/ics/%s-series/data_Gadget/Jobs/"%(options.get("paths","home_path"),options.get("paths","repository_path"),options.get("series","series_name"))
			gadget_script_filename = "jobsubmitQ_Gadget_%s-series.sh"%options.get("series","series_name")
			file(gadget_script_directory+gadget_script_filename,"w").write(generate_BGQ_Gadget_submission(options,used_blocks,1,total_simulations))
			print "\n%s job script written!\n"%(gadget_script_directory+gadget_script_filename)
			os.chmod(gadget_script_directory+gadget_script_filename,stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)

		elif(answer=="y"):

			#Prompt in how many parts user want to split the batch
			print "In how many parts do you want to split this batch?"
			nparts = int(raw_input("-->"))

			#Generate a submission script for each part
			for i in range(nparts):
					
				print "\nPart %d: select simulations to run (remember there are) %d in total"%(i+1,total_simulations)
				print "First:"
				first = int(raw_input("-->"))
				print "Last:"
				last = int(raw_input("-->"))

				subjob_nsim = last-first+1
				subjob_needed_blocks = needed(subjob_nsim,options.getint("topology","max_sims_sub_block"))
				print "There are %d simulations in this sub-batch, you will need %d sub-blocks, please select them among:"%(subjob_nsim,subjob_needed_blocks)
				for j in used_blocks:
					print "%d --> %s"%(j,corner(j))
				print ""

				subjob_used_blocks = []

				for j in range(subjob_needed_blocks):
				
					print "Please select the id of block %d"%(j+1)
					selected_block = int(raw_input("-->"))
					while(True):
						if(selected_block in subjob_used_blocks):
							print "sub-block with id %d already selected! Choose another one!"%selected_block
							print "Selected blocks so far"
							print subjob_used_blocks
							print "Please select the id of block %d"%(j+1)
							selected_block = int(raw_input("-->"))
						else:
							subjob_used_blocks.append(selected_block)
							break

				print "\nThese are the sub-blocks you selected for part %d:"%(i+1)
				for j in subjob_used_blocks:
					print "%d --> %s"%(j,corner(j))
				print ""

				#Now the sub-blocks are selected, we can proceed in writing the submission script
				gadget_script_directory = "%s/%s/localStorage/ics/%s-series/data_Gadget/Jobs/"%(options.get("paths","home_path"),options.get("paths","repository_path"),options.get("series","series_name"))
				gadget_script_filename = "jobsubmitQ_Gadget_%s-series_%d.sh"%(options.get("series","series_name"),i+1)
				file(gadget_script_directory+gadget_script_filename,"w").write(generate_BGQ_Gadget_submission(options,subjob_used_blocks,first,last))
				print "\n%s job script written!\n"%(gadget_script_directory+gadget_script_filename)
				os.chmod(gadget_script_directory+gadget_script_filename,stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)


		else:
			print "Please select y or n!"
			exit(1)



	else:
		print "%s is not a valid option!"%sys.argv[2]
		exit(1)

	print "\nDone!"