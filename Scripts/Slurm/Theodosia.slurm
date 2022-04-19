#!/bin/bash
#SBATCH -J FPr1Ra9
#SBATCH --time=0:15:00
#SBATCH -n 8
#SBATCH -o FPr1Ra9-%j.out
#SBATCH --signal=10@60

mpirun ./BoussinesqFPlane3DQGModel
