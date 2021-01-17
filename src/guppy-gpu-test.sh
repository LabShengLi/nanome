#!/bin/bash
#SBATCH --job-name=guppysh
#SBATCH --partition=compute
##SBATCH -p gpu
##SBATCH --gres=gpu:1
##SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=70g # memory pool for all cores
#SBATCH --time=00:50:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


conda env remove --name luna16
conda env remove --name nmfw
conda env remove --name NoduleX
