#!/bin/bash

#SBATCH -A TG-AST140041

#SBATCH -J PlanesAdaptive
#SBATCH -o planesAdaptive.out
#SBATCH -e planesAdaptive.err

#SBATCH -n 16
#SBATCH -p development
#SBATCH -t 02:00:00

#SBATCH --mail-user=apetri@phys.columbia.edu
#SBATCH --mail-type=all

ibrun -n 16 -o 0 /opt/apps/intel14/mvapich2_2_0/python/2.7.6/lib/python2.7/site-packages/mpi4py/bin/python-mpi make_planes_adaptive.py
