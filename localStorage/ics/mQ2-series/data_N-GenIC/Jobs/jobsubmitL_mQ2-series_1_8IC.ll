# @ job_type = bluegene
# @ class = long
#
# Classes above are currently: short (1024, 2048, 3072, or 4096 nodes, 24h) , normaldyn: (512 nodes, 48h), long (32 or 128 nodes, 72h).
#
# The executable that will run your parallel application should always be specified as per the next line.
# @ executable = /bgl/BlueLight/ppcfloor/bglsys/bin/mpirun
#
# Run on 512 nodes using a dynamic partition.
# Specify partition size using the following statement. This statement is the only way a partition size 
# should ever be specified in a LoadLeveler job control file, i.e. use of bg_partition
# has been eliminated.
# It is possible to run on fewer processors than those afforded by the partition size, see the 
# description of -np in the "Notes for Sample" below. BUT NOTE THAT you must use # @ bg_size 
# and not -np to specify the partition size.  In other words, you must use # @ bg_size
# to allocate the partition. 
# Then you can optionally use -np to run on fewer processors if this is necessary 
# for your run -- but you will charged for the entire partition that you allocated.
#
# @ bg_size = 128
#
# initialdir (see the next active line) will be used as the working directory for this batch job. 
# @ initialdir = /gpfs/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_N-GenIC/Logs 
#
# If for example your jobid is 82, your output and error will be written in
# directory /home/johndoe/app1/runs, to files 82.out and 82.err respectively.
# @ input = /dev/null
# @ output = $(jobid).out
# @ error = $(jobid).err
#
# Maximum wall clock time for job will be 20 minutes.
# @ wall_clock_limit = 12:00:00
#
# Send email to johndoe@bnl.gov when job has completed.
# @ notification = complete
# @ notify_user = apetri@phys.columbia.edu
#
# Specify executable for your parallel application, and arguments to that executable.
# Note that the arguments to specify for the executable will vary depending upon the executable.
# Specify any special environment variables for your application that need to be in the environment 
# presented to the job on the compute nodes, they will vary
# depending upon the application -  some applications will not require any - so delete or modify the
# -env specification below.
# @ arguments = -np 256 -exe /gpfs/home3/a/apetri/IG_Pipeline_0.1/N-GenIC/N-GenICl \ 
-cwd /gpfs/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_N-GenIC/Parameters \ 
-mode VN \ 
-args "ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic1.param ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic2.param ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic3.param ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic4.param ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic5.param ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic6.param ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic7.param ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic8.param" 
#
# The next statement marks the end of the job step. This example is a one-job-step batch job,
# so this is equivalent to saying that the next statement marks the end of the batch job.
# @ queue
