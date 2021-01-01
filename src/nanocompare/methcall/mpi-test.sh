#!/bin/bash
#SBATCH --job-name=MPIJob    	# job name (shows up in the queue)
#SBATCH -q batch
#SBATCH --cpus-per-task=4       # 2 Physical cores per task.
#SBATCH --ntasks=2          	# number of tasks (e.g. MPI)
#SBATCH --mem-per-cpu=512MB  	# memory/cpu in MB (half the actual required memory)
#SBATCH --time=00:01:00     	# Walltime (HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

srun pwd                        # Prints  working directory