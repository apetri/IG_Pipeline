Inspector Gadget Pipeline (Version 0.2)
===============

This is a pipeline for Weak Gravitational Lensing simulations: given a set of cosmological parameters, it produces multiple realizations of convergence and shear maps that can then be used for further statistical analysis. The IG pipeline consists of three steps:
 - 3D cosmological box generation using Gadget 2
 - Projection of the 3D simulation boxes on 2D lensing planes
 - Ray tracing and production of 2D convergence and shear maps (or catalogs)

1) 3D box generation: workflow
------------------------------

This step is broken down in three four sub-steps: setting the environment and directory structure, running the CAMB code to generate the matter power spectra, running an initial condition generator, and evolve the initial conditions using Gadget. These steps are glued together by Precambrian, an application that takes care of generating parameter files and reformatting the output of a particular step making it a suitable input for the next step.

**Before you do anything**

In order to make this pipeline portable between different computer clusters, we separated the main part of the Makefiles (which is always called Makefile.main) from the system specific one (which is called Makefile.xxx). If not existent, for each step in the pipeline you must create a Makefile.xxx file that points to the compilers and libraries specific to your system. After you do that, you need to set the environmeny variable "THIS" to "xxx"

    export THIS=xxx

This will select the correct system specific Makefile once you run "make". 

**1.0) Set the environment**

A batch of N-body simulations is defined by a set of cosmological models to run, along with other parameters, such as the box size, the number of independent simulations to run, etc... All these specifications are tunable in the options section of the script build_options.py, which resides in the Precambrian directory. Once you decided a setting, and have given a name to the simulation batch (for example "mQ2"), also in the Precambrian directory, run 

    make

This will compile and link a Precambrian application suited for your simulation batch; the next step is setting up the directory structure that will host the outputs of the box generation: typically the small files, like the power spectra and the submission scripts will be hosted on a subdirectory of your home directory, called "localStorage", while the big files, such as the 3D snapshots written by Gadget, will be hosted on some other mass storage. You have to set all this path options in an ini options file, which for reference we will call "example\_options.ini"; you will notice that the "make" you ran before generated a file called "default\_options.ini", which will serve as a blueprint for your ini options file. Once you set your "example_options.ini", you will run 

    python make_directories.py example_options.ini

that will generate the directory structure for your simulation batch. Now you are ready to start! The next step will be producing the CAMB executable. 

__IMPORTANT__: The cosmological parameters to run (which are specified by different combinations of nine different parameters) MUST be written in a file (see cosmologies.txt for an example) which name must be indicated in your example_options.ini file. Moreover parameters like the number of random seeds, the number of box sizes and the number of initial redshifts must be indicated in build_options.py, and the Precambrian executable must be re-made once you change one of these. If you simply change it in your ini options file, without re-running "make", it is likely that you'll bump in a segmentation fault at runtime!

**1.1) Compile and run CAMB**

In the camb directory, run

    make

this will compile the code and produce a "camb" executable. Before running CAMB, we have to specify the parameter files that will contain the appropriate directives; these files are generated for you by Precambrian. You have to go in the Precambrian directory and run

    ./Precambrian example_options.ini 1

This will generate the CAMB parameter files in localStorage/ics/xxx-series/data_CAMB/Parameters; on a multi-core machine, usually one would run camb as 

    mpiexec -np <numTasks> ./camb params1.ini ... paramsN.ini

to generate, in parallel, N power spectra, one for each parameter file (if you don't specify a parameter file for each task, camb quits and throws an error message). In the case where you are on a computer cluster, follow the instructions in the Submissions directory in order to generate the cluster specific submission scripts. Once CAMB finishes running, the spectra are saved and stored in localStorage/ics/xxx-series/data\_CAMB/Output\_Data, and you are ready for the next step!

**1.2) Run N-GenIC: generate the initial conditions**

This step will take care of generating the initial conditions for the simulations using the power spactra computed with CAMB at the previous step. First we start with building the N-GenIC executable, by going in the N-GenIC directory and typing

    make

which will compile and link the code in an executable named N-GenIC (or the name you provided in the system specific Makefile). This executable will need its own parameter file and power spectrum format in order to run, which can be quite a pain to write by hand; luckily Precambrian will take care of this step for us. Just run, in the Precambrian directory

    ./Precambrian example_options.ini 2

to convert the CAMB power spectra in a N-GenIC suitable format (these will be saved in data\_N-GenIC/Power\_Spectra), and 

    ./Precambrian example_options.ini 3

to generate the appropriate N-GenIC parameter files (which will be written in data_N-GenIC/Parameters). Normally one would run the initial condition generator with mpiexec as usual

    mpiexec -np <numTasks> ./N-GenICq   parameters1.param   parameters2.param   ...   parametersN.param

but on a computer cluster we have to submit our runs via a job submission script. Again, look in the Submissions directory for the machine specific instructions on how to do this. Once you submit, wait till the computation finishes and the files with the initial conditions will have been written to the mass storage disk. You are now ready to run Gadget2, the actual N-body code!!

**1.3) Run Gadget for gravitational evolution**

If you got to this step, it means now you have generated the initial conditions, and you are ready to evolve them with Gadget and generate a series of 3D snapshots which will be taken during the nonlinear evolution of the dark matter particles. Before even building the Gadget executable, you will need to create a text file in the repository top directory, called "outputs\_xxx-series.txt", with a list of numbers that will represent the time instants at which the snapshots will be taken (an example called "outputs\_mQ3-series.txt" is already provided). Now it's time to build the Gadget executable, just go in the Gadget2 directory and run

    make 

This will compile and link the code to an executable called Gadget2 (or other, the name of the executable is set by you in the system specific Makefile); of course this too will need its own formatted parameter files, which are a pain to write by hand. Luckily Precambrian already did it for us when we ran it last time: the Gadget parameter files are saved in data\_Gadget/Parameters. For how to generate the submission scripts for the specific machine you are running on, again, look in the appropriate subdirectory of Submissions. Once you have submitted the job, you have to wait till it completes. When done, you are ready for step 2, the plane generation. 

**1.4) Read in a snapshot**

The directory Gadget2/readOutput contains a slight modification of the read\_snapshot.c code provided with Gadget: it consists in a library of functions (coded in read\_snapshot\_utilities.c) that read in a single Gadget snapshot, doing the right thing (skipping headers, padding, etc...); the information about the particles is stored in a heap allocated array of type struct particle\_data, which has to be freed at the end of usage. (The allocation of the array is done automatically by the call of the read\_snapshot function). A test driver (read\_snapshot.c) is provided for testing, and it can be compiled and linked with read\_snapshot\_utilities.c just typing 'make'; you can do a test run just running the read_snapshot binary that is produced. In addition, if you have ffpmeg installed, you can run 'make movie' to make a movie of the simulation snapshots to see the particle evolution in real time (at this stage, it works efficiently only for 32x32x32 particles or smaller only, though).

**1.5) 3D shapshot power spectrum measurement**

Go in the Gadget2/Power\_Spectrum\_3D directory and run

    make

this will compile and link an executable called "3D\_Power\_Spectrum\_Calculator": this is the tool you will use to measure the 3D power spectrum of a Gadget snapshot. The "make" command, will also produce a "default\_options.ini" file, that serves as blueprint for the options file you will have to pass to "3D\_Power\_Spectrum\_Calculator"; let's call "power\_spectrum.ini" this options file. Adjust it to the settings corresponding to your simulation batch (the options should be self explanatory) and run it

    mpiexec -n <number_of_tasks> ./3D_Power_Spectrum_Calculator   power_spectrum.ini 

2) Generation of lens planes
----------------------------

Before we can generate shear and convergence maps by ray tracing through the Gadget boxes we generated in the previous steps, we need to make 2D projections of the 3D boxes. This is because the ray tracing will procede in discrete steps, in which a light ray traverses subsequent density planes, from Earth to the particular galaxy and gets deflected at each step. The deflection of the light ray once it traverses a plane is proportional to the gravitational potential at that particular plane; the computation of the gravitational potential planes (lens planes) is the main goal of step 2: the code will slice each Gadget snapshot in 9 slices (3 for each direction), and compute the gravitational potential solving the Poission equation with FFTs.
What you have to do is go in the Inspector_Gadget directory and build the executable running

    make

**Important**

In this part, you should build Inspector_Gadget __without__ the openMP flags! The behaviour is not tested yet with openMP on!! Go to the Inspector_Gadget directory for further instructions on how to tune the parameter files for the execution. In order to generate the submission scripts, look in the appropriate subdirectory of Submissions: like before this step will be taken care for you. For reference, the way you run Inspector_Gadget in planes generation mode is the following

    mpirun -n <number_of_processors> ./Inspector_Gadget <arguments>

where the arguments refer to the following

    argv[1]: Number of simulations to process (simulations with different initial conditions are considered different)
    argv[2]: Number of processors per simulation (the argv[1]*argv[2]=<number_of _processors> is enforced)
    argv[3]: IG ini parameter file (see the Inspector_Gadget directory for a complete explanation of parameters)
    argv[4...4+argv[1]-1]: Cosmology identifiers (example Om0.260_Ol0.740_w-1.000_ns0.960_si0.800_ic1)

The code will run in this mode if the "mode" options in argv[3] is set to 1.    

3) Ray tracing: creating the shear maps
---------------------------------------

Once we generated and saved the lens planes, we can proceed to the final step of the pipeline, the actual ray tracing: this operation will be performed by the same Inspector_Gadget executable as before (which will operate in a different mode); you need to remake it

    make

with the openMP flags on this time! For reference, the way you run Inspector_Gadget in planes generation mode is the following

    mpirun -n <number_of_processors> ./Inspector_Gadget <arguments>

where the arguments refer to the following

    argv[1]: Number of different cosmological models to process (Warning: tested with 1 only, more cosmologies at once coming soon...)
    argv[2]: Number of processors per simulation (the argv[1]*argv[2]=<number_of _processors> is enforced)
    argv[3]: IG ini parameter file (see the Inspector_Gadget directory for a complete explanation of parameters)
    argv[4...4+argv[1]-1]: Cosmology identifiers (example Om0.260_Ol0.740_w-1.000_ns0.960_si0.800)

The code will run in this mode if the "mode" options in argv[3] is set to 2. Once you run this executable, and it completes, this will ray trace through the potential planes produced in step 2 and produce the simulated 2D maps or galaxy catalogs that you need for your project.  

