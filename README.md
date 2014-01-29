Inspector Gadget Pipeline (Version 0.1)
===============

This is a pipeline for Weak Gravitational Lensing simulations: given a set of cosmological parameters, it produces multiple realizations of convergence and shear maps that can then be used for further statistical analysis. The IG pipeline consists of three steps:
 - 3D cosmological box generation using Gadget 2
 - Projection of the 3D simulation boxes on 2D lensing planes
 - Ray tracing and production of 2D convergence and shear maps

**0) Precambrian**

**1) 3D box generation**

This step is broken down in three smaller sub-steps: running the CAMB code to generate the matter power spectra, running an initial condition generator, and evolve the initial conditions using Gadget. 

_1.1) Running CAMB_

In the camb directory, run make to compile the code; to run it, assuming you have openmpi installed, you can type:

_make_

_mpiexec -N <numTasks> ./camb params1.ini ... paramsN.ini_

to generate, in parallel, N power spectra, one for each parameter file. If you don't specify a parameter file for each task, camb quits and throws an error message. 

_1.2) Generate the initial conditions_

_1.3) Run Gadget for gravitational evolution_

_1.4) Read in a snapshot_

The directory Gadget2/readOutput contains a slight modification of the read_snapshot.c code provided with Gadget: it consists in a library of functions (coded in read_snapshot_utilities.c) that read in a single Gadget snapshot, doing the right thing (skipping headers, padding, etc...); the information about the particles is stored in a heap allocated array of type struct particle_data, which has to be freed at the end of usage. (The allocation of the array is done automatically by the call of the read_snapshot function). A test driver (read_snapshot.c) is provided for testing, and it can be compiled and linked with read_snapshot_utilities.c just typing 'make'; you can do a test run just running the read_snapshot binary that is produced. In addition, if you have ffpmeg installed, you can run 'make movie' to make a movie of the simulation snapshots to see the particle evolution in real time (at this stage, it works efficiently only for 32x32x32 particles or smaller only, though).    

**2) Generation of lens planes**

**3) Ray tracing**

