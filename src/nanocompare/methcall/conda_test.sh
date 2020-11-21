#!/bin/bash
#SBATCH --job-name=basecall.npmeth
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 5-14:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


source /projects/liuya/workspace/tcgajax/nanocompare/methcall/conda_setup.sh

set -x

set +x; conda activate nanoai; set -x

read_fast5_basecaller.py

set +x; conda deactivate; set -x


