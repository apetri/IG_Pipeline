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
mass_storage_folders = ["ics","snapshots"]
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

###########Create mass Storage directory structure####################
massStorage = options.get("paths","mass_storage_path") + "/Storage"
try:
	print "Creating %s"%massStorage
	os.mkdir(massStorage)
except OSError:
	print "%s already exists!"%massStorage

try:
	folder = massStorage+"/sims"
	print "Creating %s"%folder
	os.mkdir(folder)
except OSError:
	print "%s already exists!"%folder

for type_name in mass_storage_folders:
	
	folder = massStorage+"/sims"+"/%s"%type_name
	try:
		print "Creating %s"%folder
		os.mkdir(folder)
	except OSError:
		print "%s already exists!"%folder

	folder = folder + "/%s-series"%options.get("series","series_name")
	try:
		print "Creating %s"%folder
		os.mkdir(folder)
	except OSError:
		print "%s already exists!"%folder

folder = massStorage+"/sims/%s/%s-series/data_Gadget"%(mass_storage_folders[0],options.get("series","series_name"))
try:
	print "Creating %s"%folder
	os.mkdir(folder)
except OSError:
	print "%s already exists!"%folder


folder = massStorage+"/sims/%s/%s-series/data_Gadget/IC_Files"%(mass_storage_folders[0],options.get("series","series_name"))
try:
	print "Creating %s"%folder
	os.mkdir(folder)
except OSError:
	print "%s already exists!"%folder
	

print 'Done!'
