#!/bin/bash
#SBATCH --job-name=Univer.MethStat
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 100g # memory pool for all cores
#SBATCH -t 8:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#cd $PBS_O_WORKDIR # change the drectory for the one, from which the job was submitted

set -x

echo "DeepSignal_calls: $DeepSignal_calls"
echo "Tombo_calls: $Tombo_calls"
echo "Nanopolish_calls: $Nanopolish_calls"
echo "DeepMod_calls: $DeepMod_calls"
echo "bgTruth: $bgTruth"
echo "RunPrefix: $RunPrefix"
echo "parser: $parser"
echo "minCov: $minCov"
echo "tool: $tool"

mkdir -p $(dirname $RunPrefix)

baseDir=/projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats
pythonFn=Universal_meth_stats_evaluation.py

prjBaseDir=/projects/li-lab/yang/workspace/nano-compare/src

PYTHONPATH=${prjBaseDir} python ${baseDir}/${pythonFn} $DeepSignal_calls $Tombo_calls $Nanopolish_calls \
			$DeepMod_calls $bgTruth $RunPrefix $parser $minCov $tool


