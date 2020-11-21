#!/bin/bash
#SBATCH --job-name=meth-corr
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 100g # memory pool for all cores
#SBATCH -t 2-10:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

## Pass all following vars to the py script
set -x

echo "DeepSignal_calls: $DeepSignal_calls"
echo "Tombo_calls: $Tombo_calls"
echo "Nanopolish_calls: $Nanopolish_calls"
echo "DeepMod_calls: $DeepMod_calls"
echo "DeepMod_cluster_calls: $DeepMod_cluster_calls"
echo "bgTruth: $bgTruth"
echo "parser: $parser"
echo "RunPrefix: $RunPrefix"


# mkdir $(dirname $RunPrefix)
set -x

baseDir=/projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats
pythonFn=Methylation_correlation_plotting.py

python ${baseDir}/${pythonFn} $DeepSignal_calls $Tombo_calls $Nanopolish_calls \
			$DeepMod_calls $DeepMod_cluster_calls $bgTruth $parser $RunPrefix


