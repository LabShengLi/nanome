#!/bin/bash
#SBATCH --job-name=meth-common
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=300g # memory pool for all cores
#SBATCH -t 72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


set -x

pythonFile=eval_common.py

python ${pythonFile} $@
