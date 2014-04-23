import sys

#This version works with numpy only, maybe fix in the future
import numpy as np

cosmo = np.loadtxt(sys.argv[1],comments="#")

if(len(cosmo.shape)==1):
	N = 1
else:
	N = cosmo.shape[0]

cosmoFile = file(sys.argv[2],"a")

cosmoFile.write("%d\n"%N)
np.savetxt(cosmoFile,cosmo,fmt="%.3e")

cosmoFile.close()