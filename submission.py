import sys,os,stat,ConfigParser
import StringIO

#This function returns the BGQ corner ID given the corner index
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

#This function checks how many simulations should be run in total, and prints according prompts to screen
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
		
		print "You will need %d sub-blocks for this batch"%needed_blocks
		print "You will need to select the sub-blocks you want to use among these (make sure they are free)"
		for i in range(1,num_sub_blocks+1):
			print "%d --> %s"%(i,corner(i))
		return needed_blocks  

#This function generates CAMB submission script for BGQ
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




#Here we write the submission scripts to the appropriate folders

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
		print "Generating Gadget2 submission script..."
		num_simulations_check(options)


	else:
		print "%s is not a valid option!"%sys.argv[2]
		exit(1)

	print "Done!"