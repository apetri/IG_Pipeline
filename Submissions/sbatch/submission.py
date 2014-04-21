import sys,os
import ConfigParser
import StringIO

qsys = "SBATCH"

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

	########Submission script directives########

	#Job name#
	S.write("""
##########################################
#############Directives###################
##########################################

#%s -J %sCAMB\n
"""%(qsys,options.get("user","username")))

	#Output and error logs#
	logPath = "%s/IG_Pipeline_0.1/localStorage/ics/%s-series/data_CAMB/Logs/"%(options.get("user","home"),options.get("series","series"))

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

	S.write("""cd %s/IG_Pipeline_0.1/localStorage/ics/%s-series/data_CAMB/Output_Data\n\n"""%(options.get("user","home"),options.get("series","series")))

	for i in range(nCores):

		cosmology_id = cosmologies[i].strip("\n")
		parameterFile = "params_%s-%s.ini"%(options.get("series","series"),cosmology_id)
		homePath = options.get("user","home")

		S.write("""ibrun -n 1 -o %d %s/camb/%s %s/IG_Pipeline_0.1/localStorage/ics/%s-series/data_CAMB/Parameters/%s &\n"""%(i,homePath,options.get("camb","executable"),homePath,options.get("series","series"),parameterFile))

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

	########Submission script directives########

	#Job name#
	S.write("""
##########################################
#############Directives###################
##########################################

#%s -J %sN-GenIC\n
"""%(qsys,options.get("user","username")))

	#Output and error logs#
	logPath = "%s/IG_Pipeline_0.1/localStorage/ics/%s-series/data_N-GenIC/Logs/"%(options.get("user","home"),options.get("series","series"))

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

	parameterDir = "%s/IG_Pipeline_0.1/localStorage/ics/%s-series/data_N-GenIC/Parameters"%(options.get("user","home"),options.get("series","series"))
	executable = "%s/IG_Pipeline_0.1/N-GenIC/%s"%(options.get("user","home"),options.get("ngenic","executable"))

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

		S.write("""ibrun -n %d -o 0 %s %s
"""%(options.getint("ngenic","cores"),executable,ngenic_args))

	S.seek(0)
	return S.read()


#######################################################
############Gadget2 submission script##################
#######################################################
def generateGadgetSubmission(options,models,breakdown_parts):

	S = StringIO.StringIO()
	S.write("""#!/bin/bash\n\n""")

	########Submission script directives########

	#Job name#
	S.write("""
##########################################
#############Directives###################
##########################################

#%s -J %sGadget\n
"""%(qsys,options.get("user","username")))

	#Output and error logs#
	logPath = "%s/IG_Pipeline_0.1/localStorage/ics/%s-series/data_Gadget/Logs/"%(options.get("user","home"),options.get("series","series"))

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

	parameterDir = "%s/IG_Pipeline_0.1/localStorage/ics/%s-series/data_Gadget/Parameters"%(options.get("user","home"),options.get("series","series"))
	executable = "%s/IG_Pipeline_0.1/Gadget2/%s"%(options.get("user","home"),options.get("gadget","executable"))

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
			S.write("""ibrun -n %d -o 0 %s %d %d %s\n"""%(options.getint("gadget","cores_per_sim")*len(gadget_args),executable,len(gadget_args),options.getint("gadget","cores_per_sim"),arg_string))

	#Done writing the script, go ahead and return
	S.seek(0)
	return S.read()



#######################################################
############Inspector Gadget submission script#########
#######################################################

####################################################
###############Main execution#######################
####################################################