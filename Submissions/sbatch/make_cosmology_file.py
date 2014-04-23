import sys
import numpy as np

if(len(sys.argv)<2):
	print "Usage: python %s <file_with_cosmology_parameters>"%sys.argv[0]
	exit(1)

cosmopar = np.loadtxt(sys.argv[1],comments="#")

if(len(cosmopar.shape)==1):
	Ncosmo = 1
else:
	Ncosmo = cosmopar.shape[0]

print "What's the filename you want for the CAMB cosmology file?"
cambFileName = raw_input("-->")

print "What's the filename you want for the N-GenIC/Gadget cosmology file?"
gadgetFileName = raw_input("-->")

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
	cambFile.write("%s\n"%Id)

for Id in gadgetIds:
	gadgetFile.write("%s %d %d\n"%(Id,firstIC,lastIC))

cambFile.close()
gadgetFile.close()

print "\nDONE!\n"