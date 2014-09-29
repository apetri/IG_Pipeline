#######################################################################
################Cut lens planes out of a Gadget2 snapshot##############
#######################################################################
from __future__ import division

import os
from lenstools.simulations import Gadget2Snapshot,PotentialPlane
from lenstools.utils import MPIWhirlPool
import numpy as np
from astropy.units import deg

from mpi4py import MPI

########################################################################

#TODO options are hardcoded for now
snapshot_path = "/scratch/02918/apetri/Storage/sims/snapshots/cfhtcov-series/cfhtcov-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.800_ic1"
snapshot_file = "snapshot_"
save_path = "/work/02918/apetri/Planes"

#########################################################################

#Split the communicator to handle multiple snapshots
comm = MPI.COMM_WORLD
newcomm = comm.Split(color=comm.rank//16,key=comm.rank)

#Initialize MPI environment
try:
	pool = MPIWhirlPool(comm=newcomm)
except:
	pool = None

#Open the snapshot
snap = Gadget2Snapshot.open(os.path.join(snapshot_path,snapshot_file+"{0:03d}".format(comm.rank//16)),pool=pool)

if pool is not None:
	print("Rank {0} reading snapshot from {1}".format(comm.rank,snap.header["files"][0]))

#Get the positions of the particles
snap.getPositions()

#Decide where to cut the lenses
cut_points = np.array([70.0,140.0,210.0]) * snap.Mpc_over_h
thickness = 10.0*snap.Mpc_over_h

#Cut the lenses
for cut,pos in enumerate(cut_points):
	for normal in range(3):

		if pool is not None and pool.is_master():
			print("Cutting plane at {0} with normal {1}, of size {2} x {2}".format(pos,normal,2.9*deg))

		#Do the cutting
		plane,resolution = snap.cutLens(normal=normal,center=pos,thickness=thickness,left_corner=np.zeros(3)*snap.Mpc_over_h,plane_size=2.9*deg,plane_resolution=512,thickness_resolution=1,smooth=2,kind="potential")

		if pool is None or pool.is_master():
			
			#Wrap the plane in a PotentialPlane object
			potential_plane = PotentialPlane(plane,angle=2.9*deg,redshift=snap.header["redshift"],cosmology=snap.cosmology)

			#Save the result
			plane_file = os.path.join(save_path,"snap{0}_potentialPlane{1}_normal{2}.fits".format(comm.rank//16,cut,normal))
			print("Saving plane to {0}".format(plane_file))
			potential_plane.save(plane_file)
			
			
		if pool is not None:
			
			#Safety barrier sync
			pool.comm.Barrier()


#Close the snapshot
snap.close()
