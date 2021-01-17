#!/bin/bash
#SBATCH --job-name=density.plot
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


# bash density_plot.sh -i results/K562_WGBS_Joined/

python density_plot.py $@