import sys,os,stat
import ConfigParser
import StringIO

import math

qsys = "SBATCH"
starter = "ibrun"

##################################################################
##########Parser for options contained in INI file################
##################################################################
def parseOptions(filename):
	
	optfile = file(filename,"r")
	
	options = ConfigParser.RawConfigParser()
	options.readfp(optfile) 

	optfile.close()

	return options

##################################################################
##########Check number of simulations in input####################
##################################################################
def numSimulationsCheck(models):

	numSims = 0
	cosmo_id = []
	first_ic = []
	last_ic = []
	
	for model in models:

		model = model.strip("\n")
		model = model.split(" ")
		cosmo_id.append(model[0])
		first_ic.append(int(model[1]))
		last_ic.append(int(model[2]))

		numSims += int(model[2]) - int(model[1]) + 1

	return numSims,cosmo_id,first_ic,last_ic



####################################################
############CAMB submission script##################
####################################################
def generateCAMBSubmission(options):

	S = StringIO.StringIO()
	S.write("""#!/bin/bash\n\n""")

	repositoryPath = options.get("user","IG_repository")

	########Submission script directives########

	#Job name#
	S.write("""
##########################################
#############Directives###################
##########################################

#%s -J %sCAMB\n
"""%(qsys,options.get("user","username")))

	#Output and error logs#
	logPath = "%s/%s/localStorage/ics/%s-series/data_CAMB/Logs/"%(options.get("user","home"),repositoryPath,options.get("series","series"))

	S.write("""#%s -o %s%sCAMB.o%%j\n"""%(qsys,logPath,options.get("user","username")))
	S.write("""#%s -e %s%sCAMB.e%%j\n\n"""%(qsys,logPath,options.get("user","username")))

	#Computing resources#
	cosmofile = file(options.get("camb","models_file"),"r")
	cosmologies = cosmofile.readlines()
	cosmofile.close()

	nCores = len(cosmologies)

	S.write("""#%s -n %d\n"""%(qsys,nCores))
	S.write("""#%s -p %s\n"""%(qsys,options.get("camb","queue")))
	S.write("""#%s -t %s\n\n"""%(qsys,options.get("camb","wall_time")))

	#Email notifications
	S.write("""#%s --mail-user=%s\n"""%(qsys,options.get("user","email")))
	S.write("""#%s --mail-type=all\n"""%qsys)

	########Execution directives##########
	S.write("""

###################################################
#################Execution#########################
###################################################

""")

	S.write("""cd %s/%s/localStorage/ics/%s-series/data_CAMB/Output_Data\n\n"""%(options.get("user","home"),repositoryPath,options.get("series","series")))

	for i in range(nCores):

		cosmology_id = cosmologies[i].strip("\n")
		parameterFile = "params_%s-%s.ini"%(options.get("series","series"),cosmology_id)
		homePath = options.get("user","home")

		S.write("""%s -n 1 -o %d %s/%s/camb/%s %s/%s/localStorage/ics/%s-series/data_CAMB/Parameters/%s &\n"""%(starter,i,homePath,repositoryPath,options.get("camb","executable"),homePath,repositoryPath,options.get("series","series"),parameterFile))

	S.write("""
wait
""")

	S.seek(0)
	return S.read()

#######################################################
############N-GenIC submission script##################
#######################################################
def generateNgenICSubmission(options):

	S = StringIO.StringIO()
	S.write("""#!/bin/bash\n\n""")

	repositoryPath = options.get("user","IG_repository")

	########Submission script directives########

	#Job name#
	S.write("""
##########################################
#############Directives###################
##########################################

#%s -J %sN-GenIC\n
"""%(qsys,options.get("user","username")))

	#Output and error logs#
	logPath = "%s/%s/localStorage/ics/%s-series/data_N-GenIC/Logs/"%(options.get("user","home"),repositoryPath,options.get("series","series"))

	S.write("""#%s -o %s%sN-GenIC.o%%j\n"""%(qsys,logPath,options.get("user","username")))
	S.write("""#%s -e %s%sN-GenIC.e%%j\n\n"""%(qsys,logPath,options.get("user","username")))

	#Computing resources#
	cosmofile = file(options.get("ngenic","models_file"),"r")
	cosmologies = cosmofile.readlines()
	cosmofile.close()

	S.write("""#%s -n %d\n"""%(qsys,options.getint("ngenic","cores")))
	S.write("""#%s -p %s\n"""%(qsys,options.get("ngenic","queue")))
	S.write("""#%s -t %s\n\n"""%(qsys,options.get("ngenic","wall_time")))

	#Email notifications
	S.write("""#%s --mail-user=%s\n"""%(qsys,options.get("user","email")))
	S.write("""#%s --mail-type=all\n"""%qsys)

	########Execution directives##########
	S.write("""

###################################################
#################Execution#########################
###################################################

""")

	parameterDir = "%s/%s/localStorage/ics/%s-series/data_N-GenIC/Parameters"%(options.get("user","home"),repositoryPath,options.get("series","series"))
	executable = "%s/%s/N-GenIC/%s"%(options.get("user","home"),repositoryPath,options.get("ngenic","executable"))

	S.write("""
cd %s

"""%parameterDir)
	
	for cosmology in cosmologies:

		ngenic_args = ""
		
		cosmology = cosmology.strip("\n")
		cosmology = cosmology.split(" ")
		cosmology_id,first_ic,last_ic = cosmology[0],int(cosmology[1]),int(cosmology[2])

		for i in range(first_ic,last_ic+1):

			ngenic_args += "ics_%s-%db%d_%s_ic%d.param "%(options.get("series","series"),options.getint("series","particles_side"),options.getint("series","box_size_mpc"),cosmology_id,i)

		S.write("""%s -n %d -o 0 %s %s
"""%(starter,options.getint("ngenic","cores"),executable,ngenic_args))

	S.seek(0)
	return S.read()


#######################################################
############Gadget2 submission script##################
#######################################################
def generateGadgetSubmission(options,models,breakdown_parts):

	S = StringIO.StringIO()
	S.write("""#!/bin/bash\n\n""")

	repositoryPath = options.get("user","IG_repository")

	########Submission script directives########

	#Job name#
	S.write("""
##########################################
#############Directives###################
##########################################

#%s -J %sGadget\n
"""%(qsys,options.get("user","username")))

	#Output and error logs#
	logPath = "%s/%s/localStorage/ics/%s-series/data_Gadget/Logs/"%(options.get("user","home"),repositoryPath,options.get("series","series"))

	S.write("""#%s -o %s%sGadget.o%%j\n"""%(qsys,logPath,options.get("user","username")))
	S.write("""#%s -e %s%sGadget.e%%j\n\n"""%(qsys,logPath,options.get("user","username")))

	#Computing resources#

	#First look at how many simulations are there to run#
	numSims,cosmo_id,first_ic,last_ic = numSimulationsCheck(models)
	if(breakdown_parts>numSims):
		raise ValueError("breakdown_parts must always be less or equal to the total number of simulations!!")

	#Construct list with all the directory root names that refer to the simulations to run
	dirRoots = []
	for i in range(len(cosmo_id)):
		for j in range(first_ic[i],last_ic[i]+1):
			dirRoots.append("%s-%db%d_%s_ic%d"%(options.get("series","series"),options.getint("series","particles_side"),options.getint("series","box_size_mpc"),cosmo_id[i],j))

	#Directory root names are available here#

	#Create directories in which the Gadget snapshots will be saved#
	for root in dirRoots:
		dirToMake = "%s/Storage/sims/snapshots/%s-series/%s"%(options.get("user","scratch"),options.get("series","series"),root)
		try:
			os.mkdir(dirToMake)
		except OSError:
			print "%s already exists, or you don't have write privileges on %s"%(dirToMake,options.get("user","scratch"))

	#Directories created, now we can proceed with the execution part of the script
	simsPerPart = numSims/breakdown_parts + cmp(numSims%breakdown_parts,0)
	nCores = options.getint("gadget","cores_per_sim") * simsPerPart

	#Request appropriate number of cores
	S.write("""#%s -n %d\n"""%(qsys,nCores))
	S.write("""#%s -p %s\n"""%(qsys,options.get("ngenic","queue")))
	S.write("""#%s -t %s\n\n"""%(qsys,options.get("ngenic","wall_time")))

	#Email notifications
	S.write("""#%s --mail-user=%s\n"""%(qsys,options.get("user","email")))
	S.write("""#%s --mail-type=all\n"""%qsys)

	########Execution directives##########
	S.write("""

###################################################
#################Execution#########################
###################################################

""")

	parameterDir = "%s/%s/localStorage/ics/%s-series/data_Gadget/Parameters"%(options.get("user","home"),repositoryPath,options.get("series","series"))
	executable = "%s/%s/Gadget2/%s"%(options.get("user","home"),repositoryPath,options.get("gadget","executable"))

	S.write("""
cd %s

"""%parameterDir)

	#Proceed to write the execution commands (which will be executed in series)
	for part in range(breakdown_parts):

		gadget_args = []
		arg_string = ""

		for i in range(simsPerPart):

			try:
				gadget_args.append(dirRoots.pop(0)+".param ")
			except IndexError:
				break

		for arg in gadget_args:
			arg_string += "%s "%arg

		if(arg_string!=""):
			S.write("""%s -n %d -o 0 %s %d %d %s\n"""%(starter,options.getint("gadget","cores_per_sim")*len(gadget_args),executable,len(gadget_args),options.getint("gadget","cores_per_sim"),arg_string))

	#Done writing the script, go ahead and return
	S.seek(0)
	return S.read()



###############################################################
############Inspector Gadget submission script: planes#########
###############################################################
def generatePlanesSubmission(options,models,blockSize):

	S = StringIO.StringIO()
	S.write("""#!/bin/bash\n\n""")

	repositoryPath = options.get("user","IG_repository")

	########Submission script directives########

	#Job name#
	S.write("""
##########################################
#############Directives###################
##########################################

#%s -J %sIGPlanes\n
"""%(qsys,options.get("user","username")))

	#Output and error logs#
	logPath = "%s/%s/localStorage/ics/%s-series/data_Inspector_Gadget/Logs/"%(options.get("user","home"),repositoryPath,options.get("series","series"))

	S.write("""#%s -o %s%sIGPlanes.o%%j\n"""%(qsys,logPath,options.get("user","username")))
	S.write("""#%s -e %s%sIGPlanes.e%%j\n\n"""%(qsys,logPath,options.get("user","username")))

	#Parse Inspector Gadget options file#
	IG_options = parseOptions(options.get("raytracing","IG_parameter_file"))

	#Throw error if IG parameter file is configured for the wrong mode
	if(IG_options.getint("mode","mode") != 1):
		raise ValueError("The parameter file you supplied is not configured in mode 1 (Plane generation)! Quitting...")

	#Decide the number of snapshots to read and the seed block
	nSnapshots = IG_options.getint("i/o_amount","global_last_snapshot") - IG_options.getint("i/o_amount","global_first_snapshot") + 1

	#Check models file for cosmological models to generate the planes of
	numSims,cosmo_id,first_ic,last_ic = numSimulationsCheck(models) 

	#Construct list with all the directory root names that refer to the simulations to run
	dirRoots = []
	for i in range(len(cosmo_id)):
		for j in range(first_ic[i],last_ic[i]+1):
			dirRoots.append("%s-%db%d_%s_ic%d"%(options.get("series","series"),options.getint("series","particles_side"),options.getint("series","box_size_mpc"),cosmo_id[i],j))

	#Directory root names are available here#

	#Create directories in which the lensing Planes will be saved#
	for root in dirRoots:
		
		dirToMake = "%s/Storage/wl/IG/%s-series/%s"%(options.get("user","scratch"),options.get("series","series"),root)
		try:
			os.mkdir(dirToMake)
		except OSError:
			print "%s already exists, or you don't have write privileges on %s"%(dirToMake,options.get("user","scratch"))

		dirToMake = "%s/Storage/wl/IG/%s-series/%s/Planes"%(options.get("user","scratch"),options.get("series","series"),root)
		try:
			os.mkdir(dirToMake)
		except OSError:
			print "%s already exists, or you don't have write privileges on %s"%(dirToMake,options.get("user","scratch"))

	#Directories created, now we can proceed with the execution part of the script

	#Decide number of cores and nodes to request, based on I/O
	nCores = nSnapshots * blockSize

	#Give IG a little more memory than what it needs, for safety (this line may be changed as needed)
	nNodes = int(nCores/16.0 + 2*(nCores/236.0))

	#Request resources accordingly
	S.write("""#%s -n %d\n"""%(qsys,nCores))
	S.write("""#%s -N %d\n"""%(qsys,nNodes))
	S.write("""#%s -p %s\n"""%(qsys,options.get("raytracing","queue")))
	S.write("""#%s -t %s\n\n"""%(qsys,options.get("raytracing","wall_time")))

	#Email notifications
	S.write("""#%s --mail-user=%s\n"""%(qsys,options.get("user","email")))
	S.write("""#%s --mail-type=all\n"""%qsys)

	########Execution directives##########
	S.write("""

###################################################
#################Execution#########################
###################################################

""")

	parameterDir = "%s/%s/localStorage/ics/%s-series/Inspector_Gadget/Parameters"%(options.get("user","home"),repositoryPath,options.get("series","series"))
	executable = "%s/%s/Inspector_Gadget/%s"%(options.get("user","home"),repositoryPath,options.get("raytracing","executable"))

	S.write("""
cd %s

"""%parameterDir)

	#Run commands
	while(True):

		IG_args = ""
		j = 0
		broken = False

		for i in range(blockSize):

			try:
				root = dirRoots.pop(0)
				IG_args += "%s "%root
				j += 1
			except IndexError:
				broken = True
				break

		if j>0:
			S.write("""%s -n %d -o 0 %s %d %d %s %s\n"""%(starter,j*nSnapshots,executable,j,nSnapshots,options.get("raytracing","IG_parameter_file"),IG_args))

		if broken:
			break

	#Done generating script, return
	
	S.seek(0)
	return S.read()


####################################################################
############Inspector Gadget submission script: ray tracing#########
############Note: this works for one cosmology at a time############
####################################################################
def generateRaySubmission(options,models):

	S = StringIO.StringIO()
	S.write("""#!/bin/bash\n\n""")

	repositoryPath = options.get("user","IG_repository")

	########Submission script directives########

	#Job name#
	S.write("""
##########################################
#############Directives###################
##########################################

#%s -J %sIGRay\n
"""%(qsys,options.get("user","username")))

	#Output and error logs#
	logPath = "%s/%s/localStorage/ics/%s-series/data_Inspector_Gadget/Logs/"%(options.get("user","home"),repositoryPath,options.get("series","series"))

	S.write("""#%s -o %s%sIGRay.o%%j\n"""%(qsys,logPath,options.get("user","username")))
	S.write("""#%s -e %s%sIGRay.e%%j\n\n"""%(qsys,logPath,options.get("user","username")))

	#Parse Inspector Gadget options file#
	IG_options = parseOptions(options.get("raytracing","IG_parameter_file"))

	#Throw error if IG parameter file is configured for the wrong mode
	if(IG_options.getint("mode","mode") != 2):
		raise ValueError("The parameter file you supplied is not configured in mode 2 (Ray tracing)! Quitting...")

	#Check how many realizations in total we want to produce: this will be the number of cores requested
	nRealizations = IG_options.getint("i/o_amount","last_realization") - IG_options.getint("i/o_amount","first_realization") + 1

	#Run on every other processor to ensure enough memory (i.e request twice as many nodes as nRealizations/16)
	nNodes = math.ceil(2*nRealizations/16.0)

	#Request resources accordingly
	S.write("""#%s -n %d\n"""%(qsys,nRealizations))
	S.write("""#%s -N %d\n"""%(qsys,nNodes))
	S.write("""#%s -p %s\n"""%(qsys,options.get("raytracing","queue")))
	S.write("""#%s -t %s\n\n"""%(qsys,options.get("raytracing","wall_time")))

	#Email notifications
	S.write("""#%s --mail-user=%s\n"""%(qsys,options.get("user","email")))
	S.write("""#%s --mail-type=all\n"""%qsys)

	########Execution directives##########
	S.write("""

###################################################
#################Execution#########################
###################################################

""")

	parameterDir = "%s/%s/localStorage/ics/%s-series/Inspector_Gadget/Parameters"%(options.get("user","home"),repositoryPath,options.get("series","series"))
	executable = "%s/%s/Inspector_Gadget/%s"%(options.get("user","home"),repositoryPath,options.get("raytracing","executable"))

	S.write("""
cd %s

"""%parameterDir)

	##########Look at models to run###########
	numSims,cosmo_id,first_ic,last_ic = numSimulationsCheck(models)

	##########Execution commands, might need to be modified in case of memory overflow problems#############
	for i in range(len(cosmo_id)):
		
		ig_arg = "%s-%db%d_%s"%(options.get("series","series"),options.getint("series","particles_side"),options.getint("series","box_size_mpc"),cosmo_id[i])

		#First create the directory which will contain the outputs
		if(IG_options.getint("mode","galaxy_catalogue_type")==0):
			
			dirToMake = "%s/Storage/wl/IG/%s-series/%s"%(options.get("user","scratch"),options.get("series","series"),ig_arg)
			try:
				os.mkdir(dirToMake)
			except OSError:
				print "%s already exists, or you don't have write privileges on %s"%(dirToMake,options.get("user","scratch"))

			dirToMake = "%s/Storage/wl/IG/%s-series/%s/Maps"%(options.get("user","scratch"),options.get("series","series"),ig_arg)
			try:
				os.mkdir(dirToMake)
			except OSError:
				print "%s already exists, or you don't have write privileges on %s"%(dirToMake,options.get("user","scratch"))

		elif(IG_options.getint("mode","galaxy_catalogue_type")==1):

			dirToMake = IG_options.get("paths","galaxy_catalogue_output_path")
			try:
				os.mkdir(dirToMake)
			except OSError:
				print "%s already exists, or you don't have write privileges on %s"%(dirToMake,options.get("user","scratch"))

			dirToMake = IG_options.get("paths","galaxy_catalogue_output_path") + "/%s"%ig_arg
			try:
				os.mkdir(dirToMake)
			except OSError:
				print "%s already exists, or you don't have write privileges on %s"%(dirToMake,options.get("user","scratch"))


		#Now write the appropriate execution command
		S.write("""%s -n %d -o 0 %s 1 %d %s %s\n"""%(starter,nRealizations,executable,nRealizations,options.get("raytracing","IG_parameter_file"),ig_arg))

	#Done generating script, return

	S.seek(0)
	return S.read()


####################################################
###############Main execution#######################
####################################################

if(__name__=="__main__"):

	#Prompt for correct number of arguments
	if(len(sys.argv)<2):
		print "Usage: python %s <ini_options_file>"%sys.argv[0]
		exit(1)

	#Parse ini options file
	options = parseOptions(sys.argv[1])
	repositoryPath = options.get("user","IG_repository")

	#Display usage instructions
	print "\nWelcome to the SLURM submission generator! Please select an operation mode"
	print ""
	print "1: Generate CAMB sumbission script"
	print "2: Generate N-GenIC submission script"
	print "3: Generate Gadget submission script"
	print "4: Generate 3D Power Spectrum Measurer submission script (not functional yet)"
	print "5: Generate Inspector Gadget submission script (Planes)"
	print "6: Generate Inspector Gadget submission script (Ray tracing)"

	mode = int(raw_input("-->"))

	#Directives for the different operation modes
	if(mode==1):

		#CAMB
		scriptFileName = "%s/%s/localStorage/ics/%s-series/data_CAMB/Jobs/%s_camb_SLURM.sh"%(options.get("user","home"),repositoryPath,options.get("series","series"),options.get("user","username"))
		
		scriptFile = file(scriptFileName,"w")
		scriptFile.write(generateCAMBSubmission(options))
		scriptFile.close()

		print "CAMB submission script generated and saved in %s!!"%scriptFileName
		print ""
		print "Do you want to sbatch-it now? (y/n)"

		answer = raw_input("-->")

		#Maybe submit the script directly?
		if(answer=="y"):
			os.execl("/usr/bin/sbatch","sbatch",scriptFileName)
		else:
			print "Goodbye! sumbission.py exited normally\n"

	elif(mode==2):

		#N-GenIC
		scriptFileName = "%s/%s/localStorage/ics/%s-series/data_N-GenIC/Jobs/%s_ngenic_SLURM.sh"%(options.get("user","home"),repositoryPath,options.get("series","series"),options.get("user","username"))

		scriptFile = file(scriptFileName,"w")
		scriptFile.write(generateNgenICSubmission(options))
		scriptFile.close()

		print "N-GenIC submission script generated and saved in %s!!"%scriptFileName
		print ""
		print "Do you want to sbatch-it now? (y/n)"

		answer = raw_input("-->")

		#Maybe submit the script directly?
		if(answer=="y"):
			os.execl("/usr/bin/sbatch","sbatch",scriptFileName)
		else:
			print "Goodbye! sumbission.py exited normally\n"

	elif(mode==3):

		#Gadget
		
		#First parse cosmological models to run from txt file
		modelFilename = options.get("gadget","models_file")
		modelFile = file(modelFilename,"r")
		models = modelFile.readlines()
		modelFile.close()

		#Prompt user for splitting job in multiple scripts
		print "There are %d models to run as indicated in %s"%(len(models),modelFilename)
		print "In how many scripts do you want to split this submission?(1-%d)"%len(models)

		nSplit = int(raw_input("-->"))
		if(nSplit>len(models)):
			"You had to type a number between 1 and %d! Quitting..."%len(models)
			exit(1)

		modelsPerPart = len(models)/nSplit + cmp(len(models)%nSplit,0)

		#Generate a submission script for each part 
		for i in range(nSplit):

			modelsSub = []
			for j in range(modelsPerPart):

				try:
					modelsSub.append(models.pop(0))
				except IndexError:
					break

			numSims,cosmo_id,first_ic,last_ic = numSimulationsCheck(modelsSub)

			#Show user which simulations will be run in this part
			print "These simulations will be run in part %d:\n"%(i+1)
			for j in range(len(cosmo_id)):
				print "%s ics:%d to %d included"%(cosmo_id[j],first_ic[j],last_ic[j])

			print "Part %d: there are %d simulations to run, do you want to further split the submission script?(y/n)"%(i+1,numSims)
			answer = raw_input("-->")

			if(answer=="y"):
				print "In how many parts? (1-%d)"%numSims
				breakdown_parts = int(raw_input("-->"))
			else:
				breakdown_parts = 1

			#Generate the script
			scriptFileName = "%s/%s/localStorage/ics/%s-series/data_Gadget/Jobs/%s_gadget_SLURM_%d.sh"%(options.get("user","home"),repositoryPath,options.get("series","series"),options.get("user","username"),i+1)

			scriptFile = file(scriptFileName,"w")
			scriptFile.write(generateGadgetSubmission(options,modelsSub,breakdown_parts))
			scriptFile.close()

			print "Gadget submission script generated and saved in %s\n!!"%scriptFileName

		print "Goodbye! sumbission.py exited normally\n"

	elif(mode==4):

		#3D power spectrum measurer
		print "Coming soon"

	elif(mode==5):

		#IG lens planes
		scriptFileName = "%s/%s/localStorage/ics/%s-series/data_Inspector_Gadget/Jobs/%s_IGPlanes_SLURM.sh"%(options.get("user","home"),repositoryPath,options.get("series","series"),options.get("user","username"))

		#Check models file for cosmological models to run
		modelFilename = options.get("raytracing","models_file")
		modelFile = file(modelFilename,"r")
		models = modelFile.readlines()
		modelFile.close()

		numSims,cosmo_id,first_ic,last_ic = numSimulationsCheck(models)

		#Inform user about which planes will be generated in this submission
		print "This submission will process %d independent simulations"%numSims
		print "This submission will generate the planes for these models:\n"
		for i in range(len(cosmo_id)):
			print "%s ICs %d to %d included"%(cosmo_id[i],first_ic[i],last_ic[i])

		print ""

		print "How many simulations do you want to process in parallel?"
		print "Note on memory scaling: to process 59 snapshots, 512x512x512 particles,5 simulations in parallel you need 20 nodes with 32GB of RAM each"
		blockSize = int(raw_input("-->"))

		#Generate script
		scriptFile = file(scriptFileName,"w")
		scriptFile.write(generatePlanesSubmission(options,models,blockSize))
		scriptFile.close()

		print "\nIG Planes submission script generated and saved in %s\n"%scriptFileName
		print ""
		print "Do you want to sbatch-it now? (y/n)"

		answer = raw_input("-->")

		#Maybe submit the script directly?
		if(answer=="y"):
			os.execl("/usr/bin/sbatch","sbatch",scriptFileName)
		else:
			print "Goodbye! sumbission.py exited normally\n"

	elif(mode==6):

		#IG ray tracing
		scriptFileName = "%s/%s/localStorage/ics/%s-series/data_Inspector_Gadget/Jobs/%s_IGRay_SLURM.sh"%(options.get("user","home"),repositoryPath,options.get("series","series"),options.get("user","username"))

		#Check models file for cosmological models to run
		modelFilename = options.get("raytracing","models_file")
		modelFile = file(modelFilename,"r")
		models = modelFile.readlines()
		modelFile.close()

		numSims,cosmo_id,first_ic,last_ic = numSimulationsCheck(models)

		#Inform user about which planes will be generated in this submission
		print "This submission will generate the maps/catalogs for these models:\n"
		for i in range(len(cosmo_id)):
			print "%s mixing ICs %d to %d included"%(cosmo_id[i],first_ic[i],last_ic[i])

		print ""

		#Generate script
		scriptFile = file(scriptFileName,"w")
		scriptFile.write(generateRaySubmission(options,models))
		scriptFile.close()

		print "\nIG Ray tracing submission script generated and saved in %s\n"%scriptFileName
		print ""
		print "Do you want to sbatch-it now? (y/n)"

		answer = raw_input("-->")

		#Maybe submit the script directly?
		if(answer=="y"):
			os.execl("/usr/bin/sbatch","sbatch",scriptFileName)
		else:
			print "Goodbye! sumbission.py exited normally\n"


	else:
		print "Mode has to be between 1 and 6! Quitting...\n"
		exit(1)
