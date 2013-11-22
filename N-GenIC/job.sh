#!/bin/tcsh
#$ -j n
#$ -cwd
#$ -pe mvapich2 4
#$ -m be
#$ -M volker@mpa-garching.mpg.de
#$ -N ics
#

module load mvapich2-1.2-sdr-gnu/4.1.2

mpiexec -np $NSLOTS ./N-GenIC ics.param


 
