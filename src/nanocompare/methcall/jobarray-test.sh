#!/bin/bash
#SBATCH --job-name=ArrayJob             # job name (shows up in the queue)
#SBATCH --partition=compute
#SBATCH --time=00:01:00                 # Walltime (HH:MM:SS)
#SBATCH --mem=512MB                     # Memory
#SBATCH --array=1-2                     # Array jobs
#SBATCH --time=00:01:00     	# Walltime (HH:MM:SS)
#SBATCH --output=log/run_%a-slurm_output.out
#SBATCH --error=log/run_%a-slurm_error.err

pwd

echo "This is result ${SLURM_ARRAY_TASK_ID}"

echo "Completed task ${SLURM_ARRAY_TASK_ID} / ${SLURM_ARRAY_TASK_COUNT} successfully"