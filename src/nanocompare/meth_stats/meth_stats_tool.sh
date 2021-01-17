#!/bin/bash
#SBATCH --job-name=meth-tool
#SBATCH -q batch
##SBATCH -p gpu
##SBATCH --gres=gpu:1             # number of gpus per node
##SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 6 # number of cores
#SBATCH --mem 280 # memory pool for all cores
#SBATCH -t 16:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

prj_dir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${prj_dir}/src/nanocompare/meth_stats/meth_stats_tool.py

mkdir -p log

python ${pythonFile} $@
