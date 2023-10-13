#!/bin/bash -l
#SBATCH --job-name=JOBNAME   # Job name
#SBATCH --output=JOBNAME.o%j # Name of stdout output file
#SBATCH --error=JOBNAME.e%j  # Name of stderr error file
#SBATCH --partition=standard    # Partition (queue) name
#SBATCH --nodes=1               # Total number of nodes
#SBATCH --ntasks-per-node=32    # Number of mpi tasks per node
#SBATCH --cpus-per-task=1       # Number of cores (threads) per task
#SBATCH --distribution=block:cyclic # consecutive task on the same node, round robin between sockets
#SBATCH --mem=0                 # Allocate all the memory on the node
#SBATCH --time=0-01:00:00       # Run time (d-hh:mm:ss)
#SBATCH --mail-type=all         # Send email at begin and end of job
#SBATCH --account=ACCOUNT  # Project for billing
#SBATCH --mail-user=EMAIL

module purge
module load LUMI
module load PrgEnv-gnu
module load cray-python
module load cray-fftw
module load cray-hdf5-parallel
module load Boost
module load buildtools
module load craype-x86-milan

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun EXECUTABLE
