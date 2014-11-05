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
snapshot_path = "/scratch/02918/apetri/Storage/sims/snapshots/dev-series/dev-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.800_ic1"
snapshot_file = "snapshot_"
save_path = "/scratch/02918/apetri/Planes512"
plane_resolution = 512
plane_size = 3.5*deg
first_snapshot = 11
last_snapshot = 58

#########################################################################

#Split the communicator to handle multiple snapshots
comm = MPI.COMM_WORLD

#Initialize MPI environment
try:
	pool = MPIWhirlPool(comm=comm)
except:
	pool = None

for n in range(first_snapshot,last_snapshot+1):

	#Open the snapshot
	snap = Gadget2Snapshot.open(os.path.join(snapshot_path,snapshot_file+"{0:03d}".format(n)),pool=pool)

	if pool is not None:
		print("Rank {0} reading snapshot from {1}".format(comm.rank,snap.header["files"][0]))

	#Get the positions of the particles
	snap.getPositions()

	#Decide where to cut the lenses
	cut_points = np.array([40.0,120.0,200.0]) * snap.Mpc_over_h
	thickness = 80.0*snap.Mpc_over_h

	#Cut the lenses
	for cut,pos in enumerate(cut_points):
		for normal in range(3):

			if pool is not None and pool.is_master():
				print("Cutting plane at {0} with normal {1},thickness {2}, of size {3} x {3}".format(pos,normal,thickness,plane_size))

			#Do the cutting
			plane,resolution,NumPart = snap.cutLens(normal=normal,center=pos,thickness=thickness,left_corner=np.zeros(3)*snap.Mpc_over_h,plane_size=plane_size,plane_resolution=plane_resolution,thickness_resolution=1,smooth=2,kind="potential")

			if pool is None or pool.is_master():
			
				#Wrap the plane in a PotentialPlane object
				potential_plane = PotentialPlane(plane,angle=plane_size,redshift=snap.header["redshift"],cosmology=snap.cosmology,num_particles=NumPart)

				#Save the result
				plane_file = os.path.join(save_path,"snap{0}_potentialPlane{1}_normal{2}.fits".format(n,cut,normal))
				print("Saving plane to {0}".format(plane_file))
				potential_plane.save(plane_file)
			
			
			if pool is not None:
			
				#Safety barrier sync
				pool.comm.Barrier()


	#Close the snapshot
	snap.close()


if pool is not None and pool.is_master():
	print("DONE!!")
