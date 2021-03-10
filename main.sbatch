#!/bin/bash
#SBATCH --job-name=nanoc.pipe
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem-per-cpu=170G # memory pool for all cores
#SBATCH --time=06:00:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#cd $SLURM_SUBMIT_DIR
#date;hostname;pwd
#
#module load singularity

#curl -fsSL get.nextflow.io | bash

#rm -rf work

./nextflow run main.nf -resume

