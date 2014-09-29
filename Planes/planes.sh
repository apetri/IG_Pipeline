#!/bin/bash

#SBATCH -A TG-AST140041

#SBATCH -J Planes
#SBATCH -o planes.out
#SBATCH -e planes.err

#SBATCH -n 960
#SBATCH -p normal
#SBATCH -t 00:10:00

#SBATCH --mail-user=apetri@phys.columbia.edu
#SBATCH --mail-type=all

ibrun -n 960 -o 0 python-mpi make_planes.py
