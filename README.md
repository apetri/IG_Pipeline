Inspector Gadget Pipeline (Version 0.1)
===============

This is a pipeline for Weak Gravitational Lensing simulations: given a set of cosmological parameters, it produces multiple realizations of convergence and shear maps that can then be used for further statistical analysis. The IG pipeline consists of three steps:
 - 3D cosmological box generation using Gadget 2
 - Projection of the 3D simulation boxes on 2D lensing planes
 - Ray tracing and production of 2D convergence and shear maps

**1) 3D box generation**

This step is broken down in three smaller sub-steps: running the CAMB code to generate the matter power spectra, running an initial condition generator, and evolve the initial conditions using Gadget. 

_1.1) Running CAMB_

In the camb directory, run make to compile the code; to run it, assuming you have openmpi installed, you can type:

_mpiexec -N <numTasks> ./camb params1.ini ... paramsN.ini_

to generate, in parallel, N power spectra, one for each parameter file. If you don't specify a parameter file for each task, camb quits and throws an error message. 

_1.2) Generate the initial conditions_

_1.3) Run Gadget for gravitational evolution_

**2) Generation of lens planes**

**3) Ray tracing**

