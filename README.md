Inspector Gadget Pipeline (Version 0.1)
===============

This is a pipeline for Weak Gravitational Lensing simulations: given a set of cosmological parameters, it produces multiple realizations of convergence and shear maps that can then be used for further statistical analysis. The IG pipeline consists of three steps:
 - 3D cosmological box generation using Gadget 2
 - Projection of the 3D simulation boxes on 2D lensing planes
 - Ray tracing and production of 2D convergence and shear maps

**1) 3D box generation: workflow**

This step is broken down in three four sub-steps: setting the environment and directory structure, running the CAMB code to generate the matter power spectra, running an initial condition generator, and evolve the initial conditions using Gadget. These steps are glued together by Precambrian, an application that takes care of generating parameter files and reformatting the output of a particular step making it a suitable input for the next step.  

_1.0) Set the environment_

A batch of N-body simulations is defined by a set of cosmological models to run, along with other parameters, such as the box size, the number of independent simulations to run, etc... All these specifications are tunable in the options section of the script build_options.py, which resides in the Precambrian directory. Once you decided a setting, and have given a name to the simulation batch (for example "mQ2"), also in the Precambrian directory, run 

make

This will compile and link a Precambrian application suited for your simulation batch; the next step is setting up the directory structure that will host the outputs of the box generation: typically the small files, like the power spectra and the submission scripts will be hosted on a subdirectory of your home directory, called "localStorage", while the big files, such as the 3D snapshots written by Gadget, will be hosted on some other mass storage. You have to set all this path options in an ini options file, which for reference we will call "example\_options.ini"; you will notice that the "make" you ran before generated a file called "default\_options.ini", which will serve as a blueprint for your ini options file. Once you set your "example_options.ini", you will run 

python make\_directories.py example\_options.ini

that will generate the directory structure for your simulation batch. Now you are ready to start! The next step will be producing the CAMB executable.  

_1.1) Compile and run CAMB_

In the camb directory, run

make

this will compile the code and produce a "camb" executable (you may have to adjust the paths to compilers or libraries in the Makefile). Before running CAMB, we have to specify the parameter files that will contain the appropriate directives; these files are generated for you by Precambrian. You have to go in the Precambrian directory and run

./Precambrian example\_options.ini 1

This will generate the CAMB parameter files in localStorage/ics/xxx-series/data_CAMB/Parameters; on a multi-core machine, usually one would run camb as 

mpiexec -N <numTasks> ./camb params1.ini ... paramsN.ini

to generate, in parallel, N power spectra, one for each parameter file (if you don't specify a parameter file for each task, camb quits and throws an error message). In the case where you are on a computer cluster (such as Blue Gene Q in this case), the CAMB runs have to be submitted to the cluster via a submission shell script. This will be taken care of for you. Notice that in the top level directory there is a python script, "submission.py", as long as a blueprint ini options file "submission\_sample\_options.ini", that will serve as a blueprint for an optons file to be passed to "submission.py", let's call it "submission\_options.ini". In this ini file you will adjust your paths as directed, you will select the block on which to run CAMB, the block corner and other options, such as the cosmological models for which to generate the power spectra. Once you are done run

python submission.py submission_options.ini 1

This will generate a submission script in localStorage/ics/xxx-series/data_CAMB/Jobs. Go in this directory and submit your job running the generated submission script

./jobsubmitQ\_CAMB\_xxx-series.sh

Once CAMB finishes running, the spectra are saved and stored in localStorage/ics/xxx-series/data\_CAMB/Output\_Data, and you are ready for the next step!

_1.2) Run N-GenIC: generate the initial conditions_

This step will take care of generating the initial conditions for the simulations using the power spactra computed with CAMB at the previous step. First we start with building the N-GenIC executable, by going in the N-GenIC directory and typing

make -f Makefile\_q

which will compile and link the code in an executable named N-GenICq. This executable will need its own parameter file and power spectrum format in order to run, which can be quite a pain to write by hand; luckily Precambrian will take care of this step for us. Just run, in the Precambrian directory

./Precambrian example\_options.ini 2

to convert the CAMB power spectra in a N-GenIC suitable format (these will be saved in data\_N-GenIC/Power\_Spectra), and 

./Precambrian example\_options.ini 3

to generate the appropriate N-GenIC parameter files (which will be written in data_N-GenIC/Parameters). Normally one would run the initial condition generator with mpiexec as usual

mpiexec -np number\_of\_tasks ./N-GenICq   parameters1.param   parameters2.param   ...   parametersN.param

but on a computer cluster such as Blue Gene Q we have to submit our runs via a job submission script. The generation of this script will be taken care of by submission.py, once you tune the appropriate knobs in submission_options.ini. You just have to run, in the top level repository

python submission.py submission\_options.ini 2

and this will generate a job submission script in data\_N-GenIC/Jobs, called jobsubmitQ\_N-GenIC\_xxx-series.sh. Go in this directory and run it

./jobsubmitQ\_N-GenIC\_xxx-series.sh

If you did everything right your job is on its way to the BGQ compute nodes! Wait till it is over and the files with the initial conditions will have been written to the mass storage disk. You are now ready to run Gadget2, the actual N-body code!!

_1.3) Run Gadget for gravitational evolution_

_1.4) Read in a snapshot_

The directory Gadget2/readOutput contains a slight modification of the read\_snapshot.c code provided with Gadget: it consists in a library of functions (coded in read\_snapshot\_utilities.c) that read in a single Gadget snapshot, doing the right thing (skipping headers, padding, etc...); the information about the particles is stored in a heap allocated array of type struct particle\_data, which has to be freed at the end of usage. (The allocation of the array is done automatically by the call of the read\_snapshot function). A test driver (read\_snapshot.c) is provided for testing, and it can be compiled and linked with read\_snapshot\_utilities.c just typing 'make'; you can do a test run just running the read_snapshot binary that is produced. In addition, if you have ffpmeg installed, you can run 'make movie' to make a movie of the simulation snapshots to see the particle evolution in real time (at this stage, it works efficiently only for 32x32x32 particles or smaller only, though).    

**2) Generation of lens planes**

**3) Ray tracing**

