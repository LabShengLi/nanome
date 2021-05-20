#!/bin/bash
#SBATCH --job-name=meth-corr
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=190g
#SBATCH --time=72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

## Pass all following vars to py script
set -x

echo "DeepSignal_calls: $DeepSignal_calls"
echo "Tombo_calls: $Tombo_calls"
echo "Nanopolish_calls: $Nanopolish_calls"
echo "DeepMod_calls: $DeepMod_calls"
echo "DeepMod_cluster_calls: $DeepMod_cluster_calls"
echo "Megalodon_calls": ${Megalodon_calls}
echo "Guppy_calls": ${Guppy_calls}
echo "bgTruth: $bgTruth"
echo "parser: $parser"
echo "RunPrefix: $RunPrefix"
echo "otherOptions:" $otherOptions

set -x

pythonFn=nanocompare/site_level_eval.py

PYTHONPATH=. python ${pythonFn} --calls DeepSignal:${DeepSignal_calls} \
	Tombo:${Tombo_calls} Nanopolish:${Nanopolish_calls} DeepMod.Cluster:${DeepMod_calls} \
	Megalodon.ZW:${Megalodon_calls} Guppy:${Guppy_calls} --bgtruth ${parser}:${bgTruth} \
	--runid MethCorr-${RunPrefix} --dsname ${Dataset} \
	--toolcov-cutoff 3 --min-bgtruth-cov 5  \
	--baseFormat 1 $otherOptions #--enable-cache --using-cache #--enable-cache --using-cache # cutoff=3 is reasonable now
