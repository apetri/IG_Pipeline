from lenstools.simulations.raytracing import RayTracer,PotentialPlane,DeflectionPlane
import numpy as np
from astropy.units import deg,rad

import os
import logging
import time

logging.basicConfig(level=logging.DEBUG)

#TODO These are hardcoded, parse from options file in the future
plane_path = "/scratch/02918/apetri/Planes512"
save_path = "/work/02918/apetri/Maps"
np.random.seed(0)

#Instantiate the RayTracer
tracer = RayTracer(lens_mesh_size=512)

start = time.time()
last_timestamp = start

#Add the lenses to the system (and perform FFT)
for i in range(11,57):
	
	plane_name = os.path.join(plane_path,"snap{0}_potentialPlane{1}_normal{2}.fits".format(i,np.random.randint(0,3),np.random.randint(0,3)))
	logging.info("Reading plane from {0}...".format(plane_name))
	tracer.addLens(PotentialPlane.load(plane_name))

now = time.time()
logging.info("Plane loading and FFT completed in {0:.3f}s".format(now-start))
last_timestamp = now

for i in range(tracer.Nlenses):
	logging.debug("Lens {0} pixels on a side {1} space {2}".format(i,tracer.lens[i].data.shape[0],tracer.lens[i].space))

#Rearrange the lenses according to redshift and roll them randomly along the axes
tracer.reorderLenses()

now = time.time()
logging.info("Reordering completed in {0:.3f}s".format(now-last_timestamp))
last_timestamp = now

tracer.randomRoll()

now = time.time()
logging.info("Rolling completed in {0:.3f}s".format(now-last_timestamp))
last_timestamp = now

#Start a bucket of light rays from these positions
b = np.linspace(0.0,tracer.lens[0].side_angle.value,512)
xx,yy = np.meshgrid(b,b)
pos = np.array([xx,yy]) * deg

#Trace the ray deflections
fin = tracer.shoot(pos,z=2.0)

now = time.time()
logging.info("Ray tracing completed in {0:.3f}s".format(now-last_timestamp))
last_timestamp = now

#Build the deflection plane
dfl = DeflectionPlane(fin.value-pos.value,angle=tracer.lens[0].side_angle,redshift=tracer.redshift[-1],cosmology=tracer.lens[0].cosmology,unit=pos.unit)

#Compute shear and convergence
conv = dfl.convergence()
shear = dfl.shear()
omega = dfl.omega()

now = time.time()
logging.info("Weak lensing calculations completed in {0:.3f}s".format(now-last_timestamp))
logging.info("Total runtime {0:.3f}s".format(now-start))

#Save the result
np.save(os.path.join(save_path,"conv.npy"),conv.data)
np.save(os.path.join(save_path,"shear.npy"),shear.data)
np.save(os.path.join(save_path,"omega.npy"),omega.data)
