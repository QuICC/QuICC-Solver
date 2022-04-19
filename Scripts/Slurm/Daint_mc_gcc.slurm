#!/bin/bash -l
#SBATCH --job-name="DynSphere"
#SBATCH --partition=normal
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --hint=nomultithread
#SBATCH --output=slurm_out.%j
#SBATCH --error=slurm_err.%j

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load daint-mc
module switch PrgEnv-cray PrgEnv-gnu
module load cray-fftw
module load cray-hdf5-parallel
module load cray-python
module load cray-tpsl
module load CMake

srun ./BoussinesqSphereDynamoExplicitModel
