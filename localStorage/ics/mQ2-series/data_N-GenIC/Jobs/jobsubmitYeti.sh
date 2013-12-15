#!/bin/sh

#Directives

#PBS -N N-GenIC
#PBS -W group_list=yetiastro
#PBS -l nodes=2:ppn=16,walltime=12:00:00,mem=131072mb
#PBS -M apetri@phys.columbia.edu
#PBS -m abe
#PBS -V

#Output and error directories

#PBS -o localhost:/u/4/a/ap3020/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_N-GenIC/Logs/
#PBS -e localhost:/u/4/a/ap3020/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_N-GenIC/Logs/

mpirun -np 256 /u/4/a/ap3020/IG_Pipeline_0.1/N-GenIC/N-GenIC /u/4/a/ap3020/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_N-GenIC/Parameters/ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic1.param

#End of script
