import os,sys,glob
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

#Options

snapshotPathDefault = '/Users/andreapetri/Documents/Cosmology_software/IG_Pipeline_0.1/Storage/sims/snapshots/mQ2-series/m-32b15_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic1'
snapshotBaseDefault = 'snapshot'
frameDirectory = sys.argv[1]
frameBase = sys.argv[2]

#Execution
snapshotPath = raw_input('Enter full path to the directory with the snapshots (default: %s) -->'%snapshotPathDefault)
if(snapshotPath == ''):
	snapshotPath = snapshotPathDefault

snapshotBase = raw_input('Enter base name for snapshot files (defalut: %s) -->'%snapshotBaseDefault)
if(snapshotBase == ''):
	snapshotBase = snapshotBaseDefault

snapshots = glob.glob('%s/%s_*'%(snapshotPath,snapshotBase))

try:
	os.mkdir(frameDirectory)
except OSError:
	print '%s already exists!'%frameDirectory

pipeName = 'mypipe'

for snapshot in snapshots:
	
	os.mkfifo(pipeName)
	pid = os.fork()

	snapshotID = snapshot.rsplit('_',1)[1]
	snapshotNumber = str(int(snapshotID))
	
	if(pid==0):
	
		os.execl('./read_snapshot.sh','read_snapshot.sh',snapshotPath,snapshotBase,snapshotNumber,pipeName)

	else:
	
		x,y,z = np.loadtxt(pipeName).transpose()/1.0e3

		os.waitpid(pid,0)
		os.remove(pipeName)

		fig = plt.figure()
		ax = fig.add_subplot(111,projection='3d')
		ax.scatter(x,y,z,s=1)
		ax.set_xlabel(r'$x(\mathrm{kpc})$')
		ax.set_ylabel(r'$y(\mathrm{kpc})$')
		ax.set_zlabel(r'$z(\mathrm{kpc})$')

		plt.savefig('%s/%s%s.png'%(frameDirectory,frameBase,snapshotID))
		plt.clf()