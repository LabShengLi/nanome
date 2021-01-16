#!/bin/bash
#SBATCH --job-name=install-methykit
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.figures # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

module load R/3.6.0

Rscript install-methykit.R

module unload R
