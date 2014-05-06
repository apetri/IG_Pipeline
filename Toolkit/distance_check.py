######################################################################
########This utility uses astropy to check the########################
########outputs spacings between different cosmologies################
######################################################################

import sys
import numpy as np

try:
	from astropy import cosmology
except ImportError:
	print "This utility requires the astropy package!! Install it and rerun"
	exit(1)

def planeSpacing(z,cosmoParams):

	#Set cosmological model
	obh2,Om,Ol,w0,wa,ns,As,si8,h = cosmoParams
	cosmo = cosmology.FlatwCDM(H0=h*100,Om0=Om,w0=w0,Neff=0)

	#Compute comoving distances
	distances = cosmo.comoving_distance(z)*h

	#Compute spacings between redshifts
	return distances[:len(z)-1] - distances[1:]

##################################################################################
##########################Main execution##########################################
##################################################################################

if __name__=="__main__":

	if(len(sys.argv)<3):
		print "Usage: python %s <outputs_file> <cosmological_models_file>"%sys.argv[0]
		exit(1)

	#Load scale factors from argv[1] and convert into redshifts
	a = np.loadtxt(sys.argv[1],comments="#")
	z = 1/a - 1

	#Load cosmological parameters from argv[2]
	cosmoParams = np.loadtxt(sys.argv[2],comments="#")

	#Construct output filename
	outFilename = sys.argv[2].split(".")[0] + "_distances.txt"

	#Open output filename and write first row
	outFile = file(outFilename,"w")
	outFile.write("#Obh2   Om    Ol     w0     wa   ns     As     si8    h    Mean plane spacing (h^-1 Mpc)\n\n")

	#Write mean plane spacing
	for i in range(cosmoParams.shape[0]):

		print "Processing model %d..."%(i+1)

		meanSpacing = (planeSpacing(z,cosmoParams[i]).mean()).value
		formattedOutput = tuple(cosmoParams[i]) + tuple([meanSpacing])

		outFile.write("%.4f %.3f %.3f %.3f %.3f %.3f %.2e %.3f %.3f %.3f\n"%formattedOutput)

	#Close file and exit
	outFile.close()

	print "Read %s and written plane spacings to %s\n"%(sys.argv[2],outFilename)
