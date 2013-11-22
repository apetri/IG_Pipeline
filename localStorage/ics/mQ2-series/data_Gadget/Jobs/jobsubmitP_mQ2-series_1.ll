# @ job_type = bluegene
# @ class = normal
#
# The executable that will run your parallel application should always be specified as per the next line.
# @ executable = /usr/bin/mpirun
#
# Run job on 512 compute nodes. LoadLeveler will dynamically allocate the partition.
# @ bg_size = 512
#
# initialdir will be the initial directory. LoadLeveler changes to this
# directory before running the job. If you don't specify this, it defaults to your current working directory
# at the time the job was submitted.
# File names mentioned in the batch script that do not begin with a slash ( / ) are relative to the initial
# directory.
# The initial directory must exist on both the fen and the compute nodes.
# @ initialdir = /gpfs/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_Gadget/Logs
#
# If for example your jobid is 82, your output and error will be written in
# directory /home/johndoe/app1/runs, to files 82.out and 82.err respectively.
# @ input = /dev/null
# @ output = $(jobid).out
# @ error = $(jobid).err
#
# Maximum wall clock time for job will be 48 hours.
# @ wall_clock_limit = 48:00:00
#
# Send email to me when job has completed.
# @ notification = complete
# @ notify_user = apetri@phys.columbia.edu
#
# Specify executable for your parallel application, and arguments to that executable.
# Note that the arguments to specify for the executable will vary depending upon the executable.
# Specify any special environment variables for your application that need to be in the environment 
# presented to the job on the compute nodes, they will vary
# depending upon the application -  some applications will not require any - so delete or modify the
# -env specification below.
# @ arguments = -np 2048 -exe /gpfs/home3/a/apetri/IG_Pipeline_0.1/Gadget2/Gadget2p \ 
-cwd /gpfs/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_Gadget/Parameters \ 
-mode VN \ 
-args "8 256 m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic1.param m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic2.param m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic3.param m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic4.param m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic5.param m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic6.param m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic7.param m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic8.param" 
#
# The next statement marks the end of the job step. This example is a one-job-step batch job,
# so this is equivalent to saying that the next statement marks the end of the batch job.
# @ queue
