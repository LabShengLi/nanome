#!/bin/bash
#SBATCH --job-name=meth-tool
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 30 # number of cores
#SBATCH --mem=150g # memory pool for all cores
#SBATCH --time=20:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

prj_dir=/projects/li-lab/yang/workspace/nano-compare
pythonFile=${prj_dir}/src/nanocompare/meth_stats/meth_stats_tool.py

mkdir -p log

python ${pythonFile} tombo-add-seq -i /fastscratch/liuya/nanocompare/K562-Runs/K562-Tombo-N50/K562-Tombo-N50-meth-call/K562.tombo.perReadsStats.combine.tsv --mpi --processors 30

#python ${pythonFile} $@
