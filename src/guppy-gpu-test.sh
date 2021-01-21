#!/bin/bash
#SBATCH --job-name=array.test
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=70g # memory pool for all cores
#SBATCH --time=00:50:00 # time
##SBATCH -o log/%x.%j.out # STDOUT
##SBATCH -e log/%x.%j.err # STDERR
#SBATCH -o log/%x.batch%a.%j.out # STDOUT
#SBATCH -e log/%x.batch%a.%j.err # STDERR
#SBATCH --array=297-300      # job array index

##SBATCH -p gpu
##SBATCH --gres=gpu:1
##SBATCH -q inference


echo SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}
echo SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}
echo SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}
echo SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX}
echo SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN}



#source /home/liuya/.bash_profile
#
#conda env list
#
#conda activate nanoai
#
#conda env list
#
#
#
#conda deactivate