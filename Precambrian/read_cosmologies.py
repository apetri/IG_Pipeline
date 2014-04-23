import sys

#This version works with numpy only, maybe fix in the future
import numpy as np

cosmo = np.loadtxt(sys.argv[1],comments="#")

file(sys.argv[2],"w").write("%d\n"%len(cosmo))
np.savetxt(sys.argv[2],cosmo,fmt="%.3e")