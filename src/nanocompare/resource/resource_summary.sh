#!/bin/bash
#SBATCH --job-name=res.summary
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 150g # memory pool for all cores
#SBATCH --time=02:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

##SBATCH -p gpu
##SBATCH --gres=gpu:1             # number of gpus per node
##SBATCH -q inference

set -x

python resource_summary.py $@


