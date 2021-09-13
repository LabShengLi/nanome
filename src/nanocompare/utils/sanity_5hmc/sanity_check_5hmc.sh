#!/bin/bash
#SBATCH --job-name=santi.5hmc
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300G # memory pool for all cores
#SBATCH -t 72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

NanoCompareDir=/projects/li-lab/yang/workspace/nano-compare
pythonFile=${NanoCompareDir}/src/nanocompare/utils/sanity_5hmc/sanity_check_5hmc.py

PHTHONPATH=${NanoCompareDir}/src python ${pythonFile} $@
