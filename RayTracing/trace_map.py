from lenstools.simulations.raytracing import RayTracer,PotentialPlane,DeflectionPlane
import numpy as np
from astropy.units import deg,rad

import logging
import time

logging.basicConfig(level=logging.DEBUG)

#TODO These are hardcoded, parse from options file in the future
plane_path = "/work/02918/apetri/Planes4096"
save_path = "/work/02918/apetri/Maps"

#Instantiate the RayTracer
tracer = RayTracer()

start = time.time()
last_timestamp = start

#Add the lenses to the system
for i in range(11,57):
	tracer.addLens(PotentialPlane.load(os.path.join(plane_path,"snap{0}_potentialPlane0_normal0.fits".format(i))))

now = time.time()
logging.info("Plane loading completed in {0:.3f}s".format(now-last_timestamp))
last_timestamp = now

#Rearrange the lenses according to redshift and roll them randomly along the axes
tracer.reorderLenses()
tracer.randomRoll()

now = time.time()
logging.info("Reordering and rolling completed in {0:.3f}s".format(now-last_timestamp))
last_timestamp = now

#Start a bucket of light rays from these positions
b = np.linspace(0.0,2.9,512)
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
