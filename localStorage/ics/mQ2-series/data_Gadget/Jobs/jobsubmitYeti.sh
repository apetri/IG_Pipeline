#!/bin/sh

#Directives

#PBS -N Gadget
#PBS -W group_list=yetiastro
#PBS -l nodes=2:ppn=16:ib,walltime=60:00:00,mem=131072mb
#PBS -M apetri@phys.columbia.edu
#PBS -m abe
#PBS -V

#Output and error directories

#PBS -o localhost:/u/4/a/ap3020/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_Gadget/Logs/
#PBS -e localhost:/u/4/a/ap3020/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_Gadget/Logs/

mpirun -np 256 /u/4/a/ap3020/IG_Pipeline_0.1/Gadget2/Gadget2 1 256 /u/4/a/ap3020/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_Gadget/Parameters/m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic1.param

#End of script
