#!/bin/bash
#SBATCH --job-name=Tombo.extract
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 100g # memory pool for all cores
#SBATCH -t 02:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

prj_dir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${prj_dir}/src/nanocompare/methcall/Tombo_extract_per_read_stats.py

mkdir -p log

python ${pythonFile} $@