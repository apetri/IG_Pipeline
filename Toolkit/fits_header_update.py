#################################################################
#This tool searches for a cosmological parameter header keyword##
#and updates it with the value that it finds in the filename#####
#################################################################

#Needs pyfits or astropy.io#

import sys

try:
	from astropy.io import fits
except ImportError:
	import pyfits as fits

#This is the function that performs the magic#

def headerUpdate(filename,keyword,filenameKey,hdu=0,separator="_"):

	#Open FITS file
	hdulist = fits.open(filename)

	#Find the filenameKey string in the FITS file
	spl = filename.split(separator)
	for s in spl:
		if s.find(filenameKey) != -1:
			break
	
	if s.find(filenameKey) == -1:
		raise ValueError("%s does not contain the %s keyword!"%(filename,filenameKey))

	#Compute value to update the header with
	value = float(s.strip(filenameKey))

	#Update the header and close hdu
	hdulist[hdu].header[keyword] = value
	hdulist.close()