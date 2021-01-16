#!/bin/bash
#SBATCH --job-name=rm
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 250g # memory pool for all cores
#SBATCH -t 23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.figures # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
set -x

pwd

echo ${abc}
echo ${abc1}
