#!/bin/bash
#SBATCH --job-name=Univer.MethStat
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem 280g # memory pool for all cores
#SBATCH --time=1-23:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x
echo "DeepSignal_calls: $DeepSignal_calls"
echo "Tombo_calls: $Tombo_calls"
echo "Nanopolish_calls: $Nanopolish_calls"
echo "DeepMod_calls: $DeepMod_calls"
echo "Megalodon_calls: ${Megalodon_calls}"
echo "bgTruth: $bgTruth"
echo "RunPrefix: $RunPrefix"
echo "parser: $parser"
echo "minCov: $minCov"
echo "dsname:" $dsname

prjBaseDir=/projects/li-lab/yang/workspace/nano-compare

pythonFn=${prjBaseDir}/src/nanocompare/meth_stats/Universal_meth_stats_evaluation.py

PYTHONPATH=${prjBaseDir}/src python ${pythonFn} --calls DeepSignal:$DeepSignal_calls \
	Tombo:$Tombo_calls Nanopolish:$Nanopolish_calls DeepMod:$DeepMod_calls \
	Megalodon:${Megalodon_calls} --bgtruth $parser:$bgTruth --runid MethPerf-$RunPrefix \
	--dsname $dsname --min-cov $minCov --report-joined
