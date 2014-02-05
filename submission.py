import sys,os,stat,ConfigParser
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

##########################################################################################################
#This function checks how many simulations should be run in total, and prints according prompts to screen#
##########################################################################################################
def num_simulations_check(options):
	
	cosmologies = options.options("cosmologies_gadget")
	sims_per_model = options.getint("series","simulations_per_model")
	print "You want to run %d models, %d simulations per model"%(len(cosmologies),sims_per_model)
	needed_blocks = (len(cosmologies)*sims_per_model + 1)/2
	num_sub_blocks = options.getint("topology","num_sub_blocks")
	if(needed_blocks>num_sub_blocks):
		
		print "You need more than %d sub-blocks! You will need to split up your work!"%num_sub_blocks
		return needed_blocks
	
	else:
		
		print "Each sub-block can run 2 simulations with 512x512x512 particles"
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
fi
""")

	S.seek(0)
	return S.read()

##############################################################
#This function generates the Gadget submission script for BGQ#
##############################################################
def generate_BGQ_Gadget_submission(options,used_blocks):

	S = StringIO.StringIO()

	#Set blockid
	S.write("""#!/bin/sh

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

	#Set computing resources usage (need to fix this, so far this job submission file runs two simulations per sub-block job)
	S.write("""
NUM_SIMS=2
TASKS_PER_SIM=%d

CORES_PER_NODE=%d
NUM_MPI_TASKS=%d
"""%(options.getint("computing_resources","tasks_per_simulation_gadget"),options.getint("computing_resources","cores_per_node_gadget"),options.getint("computing_resources","block_mpi_tasks_gadget")))

	#Set path for Gadget parameter files
	S.write("""
G_ROOT=$HOMEPATH/$REPOSITORY_DIR/localStorage/ics/%s-series/data_Gadget/Parameters
"""%options.get("series","series_name"))

	#Create a list with all the initial conditions file names (one for each simulation to run)
	parameter_filenames = []
	series_name = options.get("series","series_name")
	cosmologies = options.options("cosmologies_gadget")
	num_particles_side = options.getint("series","num_particles_side")
	box_size_kpc = options.getint("series","box_size_kpc")
	simulations_per_model = options.getint("series","simulations_per_model")
	for i in range(1,simulations_per_model+1):
		for cosmology_id in cosmologies:
			parameter_filenames.append("%s-%db%d_%s_ic%d.param"%(series_name,num_particles_side,box_size_kpc,options.get("cosmologies_gadget",cosmology_id),i))

	#parameter_filenames contains all the names of the Gadget parameter files; each block takes care of two of these
	for filename in parameter_filenames:
		S.write("""%s\n"""%filename)

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
		os.chmod(camb_script_directory+camb_script_filename,stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)

	elif(sys.argv[2]=="3"):
		#Gadget2 submission script
		print "Generating Gadget2 submission script...\n"
		
		#Check that the partition can handle your work
		needed_blocks = num_simulations_check(options)
		num_sub_blocks = options.getint("topology","num_sub_blocks")
		
		if(needed_blocks>num_sub_blocks):
			exit(1)
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

			#Now the sub-blocks are selected, we can proceed in writing the submission script
			gadget_script_directory = "%s/%s/localStorage/ics/%s-series/data_Gadget/Jobs/"%(options.get("paths","home_path"),options.get("paths","repository_path"),options.get("series","series_name"))
			gadget_script_filename = "jobsubmitQ_Gadget_%s-series.sh"%options.get("series","series_name")
			file(gadget_script_directory+gadget_script_filename,"w").write(generate_BGQ_Gadget_submission(options,used_blocks))
			os.chmod(gadget_script_directory+gadget_script_filename,stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)



	else:
		print "%s is not a valid option!"%sys.argv[2]
		exit(1)

	print "Done!"