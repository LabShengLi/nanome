#!/bin/bash
#SBATCH --job-name=rm
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 150g # memory pool for all cores
#SBATCH -t 13:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

cd /projects/li-lab/ctimage/yang/results-ctimage-early-phase

#rm -rf 08-* 09-* 10-* 11-*
#rm -rf 07-2*
#rm -rf 07-15 07-16 07-18 07-19


rm -rf 07-22  09-22  10-28



