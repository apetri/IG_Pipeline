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

#######################################################
############Gadget2 submission script##################
#######################################################

#######################################################
############Inspector Gadget submission script#########
#######################################################

####################################################
###############Main execution#######################
####################################################