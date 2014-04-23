import sys

#This version works with numpy only, maybe fix in the future
import numpy as np

cosmo = np.loadtxt(sys.argv[1],comments="#")

if(len(cosmo.shape)==1):
	N = 1
else:
	N = cosmo.shape[0]

file(sys.argv[2],"w").write("%d\n"%N)
np.savetxt(sys.argv[2],cosmo,fmt="%.3e")