#!/bin/bash
#SBATCH --job-name=density.plot
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM:SS)
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

python density_plot.py