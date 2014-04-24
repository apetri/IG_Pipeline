Blue Gene Q job submissions
===============

Before doing anything, remember that you should create a txt file named "<blockid>\_sub\_blocks.txt" that contains the names of all the sub-block corners you want to use (and which are specific to the block-ID which have been assigned to you). A sample file for the blockid R00-M0-N00-128, called "R00-M0-N00-128\_sub\_blocks.txt", has already been created for you. 

CAMB
---------------
The CAMB runs have to be submitted to the Blue Gene Cluster via a submission shell script. This will be taken care of for you. Notice that in this directory there is a python script, "submission.py", as long as a blueprint ini options file "submission\_sample\_options.ini", that will serve as a blueprint for an optons file to be passed to "submission.py", let's call it "submission\_options.ini". In this ini file you will adjust your paths as directed, you will select the block on which to run CAMB, the block corner and other options, such as the cosmological models for which to generate the power spectra. 


    python submission.py submission_options.ini 1

This will generate a submission script in localStorage/ics/xxx-series/data_CAMB/Jobs. Go in this directory and submit your job running the generated submission script

    ./jobsubmitQ_CAMB_xxx-series.sh

This will submit your job to the Blue Gene Q compute nodes. 

N-GenIC
---------------
The generation of the N-GenIC submission script will be taken care of by submission.py, once you tune the appropriate knobs in submission_options.ini. You just have to run, in this directory

    python submission.py submission_options.ini 2

and this will generate a job submission script in data\_N-GenIC/Jobs, called jobsubmitQ\_N-GenIC\_xxx-series.sh (or jobsubmitQ\_N-GenIC\_xxx-series\_n.sh if you prompted submission.py to split the job) . Go in this directory and run it

    ./jobsubmitQ_N-GenIC_xxx-series_n.sh (BGQ)

If you did everything right your job(s) is(are) on its(their) way to the BGQ compute nodes!

Gadget
---------------
You could in principle run the Gadget executable manually with mpiexec as usual, but Blue Gene Q requires that you submit your runs through a job submission script. This script will be generated automatically for you by submission.py. You have to be particularly careful in tuning the knobs in submission\_options.ini: the Blue Gene Q cluster connections have a complicated topology, and we need to specify the shapes and corners of the sub-blocks that make up our computing partition. In particular, in your submission\_options.ini file you need to specify the maximum number of simulations a sub block can handle; this of course depends on the size of the simulations (mainly the number of particles), so you need to choose this parameter carefully (for 512x512x512 particles this number is 2). After you do this just run

    python submission.py submission_options.ini 3

This will tell you how many sub-blocks you need for your job, and you will need to specify which ones you want to use (you have to make sure no one is using those, this script unfortunately does not check for that!); if you wish, this python script will allow you to split the simulation batch in multiple sub-batches that can be run independently, and will generate a submission script for each of these sub-batches. If submission.py completes succesfully, you will have your submission script ready in data\_Gadget/Jobs, and it will be called jobsubmitQ\_Gadget\_xxx-series\_n.sh (n is the sub-batch number, that may be absent if you run everything at once). Run it

    ./jobsubmitQ_Gadget_xxx-series_n.sh (BGQ)

and your Gadget jobs will be on their way to the Blue Gene Q compute nodes!