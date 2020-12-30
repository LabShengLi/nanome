#!/bin/bash
#SBATCH --job-name=combine.tombo
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

cat ${methCallsDir}/*.bed > ${methCallsDir}/${analysisPrefix}.tombo.combine.bed