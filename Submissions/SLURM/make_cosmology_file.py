import sys,glob,ConfigParser
import numpy as np

if(len(sys.argv)<3):
	print "Usage: python %s <precambrian_cosmologies_file> <ini_options_file>"%sys.argv[0]
	exit(1)

#Parse options
options = ConfigParser.RawConfigParser()
options.readfp(file(sys.argv[2],"r"))

cosmopar = np.loadtxt(sys.argv[1],comments="#")

if(len(cosmopar.shape)==1):
	Ncosmo = 1
else:
	Ncosmo = cosmopar.shape[0]

cambFileName = options.get("camb","models_file")
gadgetFileName = options.get("gadget","models_file")

print "What's the first initial condition you want to generate? Enter an integer:"
firstIC = int(raw_input("-->"))

print "What's the last initial condition you want to generate? Enter an integer:"
lastIC = int(raw_input("-->"))

cambIds = []
gadgetIds = []

#Construct cosmology IDs for CAMB/N-GenIC/Gadget files

if Ncosmo==1:
	cambIds.append("Om%.3f_Ol%.3f_w%.3f_ns%.3f"%(cosmopar[1],cosmopar[2],cosmopar[3],cosmopar[5]))
	gadgetIds.append("Om%.3f_Ol%.3f_w%.3f_ns%.3f_si%.3f"%(cosmopar[1],cosmopar[2],cosmopar[3],cosmopar[5],cosmopar[7]))
else:
	for i in range(Ncosmo):
		cambIds.append("Om%.3f_Ol%.3f_w%.3f_ns%.3f"%(cosmopar[i,1],cosmopar[i,2],cosmopar[i,3],cosmopar[i,5]))
		gadgetIds.append("Om%.3f_Ol%.3f_w%.3f_ns%.3f_si%.3f"%(cosmopar[i,1],cosmopar[i,2],cosmopar[i,3],cosmopar[i,5],cosmopar[i,7]))

cambIds = list(set(cambIds))

#Write the files

cambFile = file(cambFileName,"w")
gadgetFile = file(gadgetFileName,"w")

for Id in cambIds:
	#Throw user a warning if the Precambrian didn't generate this parameter file for some reason
	cambParamFile = "%s/%s/localStorage/ics/%s-series/data_CAMB/Parameters/params_%s-%s.ini"%(options.get("user","home"),options.get("user","IG_Repository"),options.get("series","series"),options.get("series","series"),Id)
	if(len(glob.glob(cambParamFile))==0): print "WARNING! The CAMB parameter file associated with %s does not exist!\n"%Id
	
	cambFile.write("%s\n"%Id)

for Id in gadgetIds:
	#Throw user a warning if the Precambrian didn't generate this parameter file for some reason
	gadgetParamFile = "%s/%s/localStorage/ics/%s-series/data_Gadget/Parameters/%s-%db%d_%s"%(options.get("user","home"),options.get("user","IG_Repository"),options.get("series","series"),options.get("series","series"),options.getint("series","particles_side"),options.getint("series","box_size_mpc"),Id)
	if(len(glob.glob(gadgetParamFile+"*"))==0): print "WARNING! The Gadget parameter file associated with %s does not exist!\n"%Id

	gadgetFile.write("%s %d %d\n"%(Id,firstIC,lastIC))

cambFile.close()
gadgetFile.close()

print "\nDONE!\n"
