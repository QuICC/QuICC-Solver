#!/bin/bash -l
#SBATCH --job-name="QuICCJob"
##SBATCH --account="s872"
#SBATCH --output=slurm_out.%j
#SBATCH --error=slurm_err.%j
#SBATCH --mail-type=ALL
##SBATCH --mail-user=forgot@tosetemail.com
#SBATCH --time=00:20:00
#SBATCH --nodes=12
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

module load daint-gpu
module switch PrgEnv-cray PrgEnv-gnu
module load cray-fftw
module load cray-hdf5-parallel
module load cray-python
module load cray-tpsl
module load CMake

srun ./BoussinesqSphereDynamoExplicitModel
