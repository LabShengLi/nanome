#!/bin/bash
#SBATCH --job-name=benchdata
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300G # memory pool for all cores
#SBATCH -t 1-23:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

projectDir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${projectDir}/src/nanocompare/benchmarkingDataConstruct.py

mkdir -p log

#rm -rf /projects/li-lab/yang/results/2021-04-15/BenchData

python ${pythonFile} $@