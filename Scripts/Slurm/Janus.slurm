#!/bin/bash
#SBATCH -J NewMod
#SBATCH --time=0:15:00
#SBATCH -N 1
#SBATCH --tasks-per-node 12
#SBATCH -o slurm-%j.out
#SBATCH --qos=janus-debug
#SBATCH --signal=10@60

# Intel
module load intel
module load impi
# GCC
#module load gcc
#module load openmpi

module load petsc/3.6.0
module load hdf5/1.8.15
module load fftw/3.3.4
module load python/3.4.3
#module load python/2.7.10
module load numpy/1.9.2
module load scipy/0.16.0

mpirun ./BoussinesqFPlane3DQGModel
