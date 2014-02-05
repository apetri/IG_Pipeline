import sys,os,ConfigParser

#Check if options file is provided
if(len(sys.argv)<2):
	print "Usage: python %s <ini_options_file>"%sys.argv[0]
	exit(1)

#Parse options from ini file
options = ConfigParser.RawConfigParser()
options.readfp(file(sys.argv[1],"r"))

#This function generates CAMB submission script for BGQ
print options.options("cosmologies_camb")