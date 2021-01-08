#!/bin/bash
#SBATCH --job-name=meth-corr
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 02:00:00 # time (D-HH:MM:SS)
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

prj_dir=/projects/li-lab/yang/workspace/nano-compare
baseDir=${prj_dir}/src/nanocompare/meth_stats
pythonFn=Methylation_correlation_plotting.py

python ${baseDir}/${pythonFn} --calls DeepSignal:$DeepSignal_calls \
	Tombo:$Tombo_calls Nanopolish:$Nanopolish_calls DeepMod:$DeepMod_calls \
	--bgtruth $parser:$bgTruth --runid $RunPrefix


