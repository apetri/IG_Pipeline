#######################################################################
################Cut lens planes out of a Gadget2 snapshot##############
#######################################################################
from __future__ import division

import os
from lenstools.simulations import Gadget2Snapshot,PotentialPlane
from lenstools.utils import MPIWhirlPool
import numpy as np
from astropy.units import deg,rad,Mpc

from mpi4py import MPI

########################################################################

#TODO options are hardcoded for now
snapshot_path = "/scratch/02918/apetri/Storage/sims/snapshots/cmb512-series/cmb512-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.800_ic1"
snapshot_file = "snapshot_"
save_path = "/scratch/02918/apetri/PlanesAdaptive"
plane_resolution = 2048
neighbors = 64
first_snapshot = 0
last_snapshot = 59

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

	#Cut the lenses
	for normal in range(3):

		if pool is not None and pool.is_master():
			print("Cutting plane with normal {0}, of size {1} x {1}".format(normal,snap.header["box_size"]))

		#Do the cutting
		plane,resolution,numPart = snap.cutPlaneAdaptive(normal=normal,left_corner=np.zeros(3)*snap.Mpc_over_h,plane_resolution=plane_resolution,neighbors=neighbors,projectAll=True,kind="potential")

		if pool is None or pool.is_master():
			
			#Wrap the plane in a PotentialPlane object
			potential_plane = PotentialPlane(plane.value,angle=snap.header["box_size"],redshift=snap.header["redshift"],cosmology=snap.cosmology,num_particles=numPart,unit=plane.unit)

			#Save the result
			plane_file = os.path.join(save_path,"snap{0}_potentialPlane{1}_normal{2}_neighbors{3}.fits".format(n,0,normal,neighbors))
			print("Saving plane to {0}".format(plane_file))
			potential_plane.save(plane_file)
			
			
		if pool is not None:
			
			#Safety barrier sync
			pool.comm.Barrier()

	
	#Close the snapshot
	snap.close()


if pool is not None and pool.is_master():
	print("DONE!!")
