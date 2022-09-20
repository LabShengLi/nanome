#!/bin/bash
#SBATCH --job-name=comb.allchrs.NA12878
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem 50g # memory pool for all cores
#SBATCH -t 20:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


set -x

pythonFile=combine_output_na12878.py
python ${pythonFile} $@
