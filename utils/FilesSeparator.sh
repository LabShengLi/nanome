#!/bin/bash
#SBATCH --job-name=seperate_files
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=80G # memory pool for all cores
#SBATCH --time=05:00:00 # time
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

FilesSeparator.py $@
