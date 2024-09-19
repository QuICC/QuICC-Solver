#!/bin/bash -l
#SBATCH --job-name="DynSphere"
#SBATCH --partition=normal.4h
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_out.%j
#SBATCH --error=slurm_err.%j

module load stack
module load openmpi fftw hdf5 boost python openblas cmake

srun ./BoussinesqShellDynamoExplicitModel
