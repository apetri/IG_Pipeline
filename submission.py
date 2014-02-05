import sys,os,stat,ConfigParser
import StringIO

#Check if options file is provided
if(len(sys.argv)<2):
	print "Usage: python %s <ini_options_file>"%sys.argv[0]
	exit(1)

#Parse options from ini file
options = ConfigParser.RawConfigParser()
options.readfp(file(sys.argv[1],"r"))

#This function generates CAMB submission script for BGQ
def generate_BGQ_camb_submission(options):
	
	S = StringIO.StringIO()
	
	S.write("""#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <BLOCKID>" >&2
    exit 1
fi""")

	S.seek(0)
	return S.read()


#Here we write the submission scripts to the appropriate folders

#CAMB submission script
camb_script_directory = "%s/%s/localStorage/ics/%s-series/data_CAMB/Jobs/"%(options.get("paths","home_path"),options.get("paths","repository_path"),options.get("series","series_name"))
camb_script_filename = "jobsubmitQ_CAMB_%s-series.sh"%options.get("series","series_name")
file(camb_script_directory+camb_script_filename,"w").write(generate_BGQ_camb_submission(options))
os.chmod(camb_script_directory+camb_script_filename,stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)