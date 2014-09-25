#######################################################################
################Cut lens planes out of a Gadget2 snapshot##############
#######################################################################

import os
from lenstools.simulations import Gadget2Snapshot
from lenstools.utils import MPIWhirlPool
import numpy as np

########################################################################

#TODO options are hardcoded for now
snapshot_path = "/scratch/02918/apetri/Storage/sims/snapshots/cfhtcov-series/cfhtcov-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.800_ic1"
snapshot_file = "snapshot_001"
save_path = "/work/02918/apetri/Planes"

#########################################################################

#Initialize MPI environment
try:
	pool = MPIWhirlPool()
except:
	pool = None

#Open the snapshot
snap = Gadget2Snapshot.open(os.path.join(snapshot_path,snapshot_file),pool=pool)

if pool.is_master():
	print("Reading snapshot from {0}".format(snap.header["files"][0]))

#Get the positions of the particles
snap.getPositions()

#Decide where to cut the lenses
cut_points = np.array([70.0,140.0,210.0]) * snap.Mpc_over_h
thickness = 10.0*snap.Mpc_over_h

#Cut the lenses
for cut,pos in enumerate(cut_points):
	for normal in range(3):

		if pool.is_master():
			print("Cutting plane at {0} with normal {1}, of size {2} x {2}".format(pos,normal,snap.lensMaxSize()))

		#Do the cutting
		plane,resolution = snap.cutLens(normal=normal,center=pos,thickness=thickness,left_corner=np.zeros(3)*snap.Mpc_over_h,plane_size=snap.lensMaxSize(),plane_resolution=4096,thickness_resolution=1,smooth=2,kind="potential")

		if pool.is_master():
			
			#Save the result
			plane_file = os.path.join(save_path,"plane{0}_normal{1}.npy".format(cut,normal))
			print("Saving plane to {0}".format(plane_file))
			np.save(plane_file,plane)


#Close the snapshot
snap.close()
