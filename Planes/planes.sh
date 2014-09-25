#!/bin/bash

#SBATCH -A TG-AST140041

#SBATCH -J Planes
#SBATCH -o planes.out
#SBATCH -e planes.err

#SBATCH -n 16
#SBATCH -p development
#SBATCH -t 01:00:00

#SBATCH --mail-user=apetri@phys.columbia.edu
#SBATCH --mail-type=all

ibrun -n 16 -o 0 python-mpi make_planes.py
