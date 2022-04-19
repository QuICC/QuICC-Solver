#!/bin/bash
#SBATCH -J ShellC1	   # job name
#SBATCH -o ShellC1.o%j     # output and error file name (%j expands to jobID)
#SBATCH -n 16            # total number of mpi tasks requested
#SBATCH -p development     	   # queue (partition) -- normal, development, etc.
#SBATCH -t 01:00:00        # run time (hh:mm:ss)
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

ibrun ./BoussinesqDynamoShellStdModel
