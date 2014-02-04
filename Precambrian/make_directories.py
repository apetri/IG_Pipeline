import sys,os
import ConfigParser

#Check if ini options file is given 
if(len(sys.argv) < 2):
	print 'Usage: python %s <ini_options_file>' %sys.argv[0]
	exit(1)

#Parse options using config parser
options = ConfigParser.RawConfigParser()
options.readfp(file(sys.argv[1],"r"))

############Useful names#############
data_folders = ["data_CAMB","data_N-GenIC","data_Gadget"]
type_folders = ["Jobs","Logs","Parameters"]
camb_output = "Output_Data"
power_spectra = "Power_Spectra"

#############Create localStorage directory structure#############
repository_location = options.get("paths","home_path") + options.get("paths","repository_relative_path")
localStorage = repository_location + "/localStorage"
try: 
	print "Creating %s"%localStorage
	os.mkdir(localStorage)
except OSError:
	print "%s already exists!"%localStorage

try:
	print "Creating %s"%localStorage+"/ics"
	os.mkdir(localStorage+"/ics")
except OSError:
	print "%s already exists!"%localStorage+"/ics"

series_folder = localStorage+"/ics/%s-series"%options.get("series","series_name")
try:
	print "Creating %s"%series_folder
	os.mkdir(series_folder)
except OSError:
	print "%s already exists!"%series_folder

for folder in data_folders:
	try:
		data_folder = series_folder+"/"+folder
		print "Creating %s"%data_folder
		os.mkdir(data_folder)
	except OSError:
		print "%s already exists!"%data_folder

for data_name in data_folders:
	for  type_name in type_folders:
		try:
			folder = series_folder+"/"+data_name+"/"+type_name
			print "Creating %s"%folder
			os.mkdir(folder)
		except OSError:
			print "%s already exists!"%folder

try:
	folder = series_folder + "/%s"%data_folders[0]+"/%s"%camb_output
	print "Creating %s"%folder
	os.mkdir(folder)
except OSError:
	print "%s already exists!"%folder

try:
	folder = series_folder + "/%s"%data_folders[1] + "/%s"%power_spectra
	print "Creating %s"%folder
	os.mkdir(folder)
except OSError:
	print "%s already exists!"%folder

print 'Done!'
