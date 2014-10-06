#!/bin/bash

#SBATCH -A TG-AST140041

#SBATCH -J RayTraceMap
#SBATCH -o raymap.out
#SBATCH -e raymap.err

#SBATCH -n 16
#SBATCH -p development
#SBATCH -t 01:00:00

#SBATCH --mail-user=apetri@phys.columbia.edu
#SBATCH --mail-type=all

#Enable MIC offload
#export MKL_MIC_ENABLE=1
#export OFFLOAD_REPORT=2
#export OMP_NUM_THREADS=16
#export MIC_OMP_NUM_THREADS=240

ibrun -n 1 -o 0 python trace_map.py
