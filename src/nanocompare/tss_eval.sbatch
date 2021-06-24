#!/bin/bash
#SBATCH --job-name=tss-eval
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=290g
#SBATCH --time=72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

##SBATCH --mem 290g # memory pool for all cores
##SBATCH -p gpu
##SBATCH --gres=gpu:1             # number of gpus per node
##SBATCH -q inference

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

echo ${NanoCompareDir}


pythonFn=nanocompare/tss_eval.py

PYTHONPATH=. python ${pythonFn} \
	--calls Nanopolish:${Nanopolish_calls} \
			Megalodon.ZW:${Megalodon_calls} \
			DeepSignal:${DeepSignal_calls} \
			Guppy:${Guppy_calls} \
			Tombo:${Tombo_calls} \
			DeepMod.Cluster:${DeepMod_calls} \
	--bgtruth ${parser}:${bgTruth} \
	--runid TSS-${RunPrefix} --dsname ${dsname} \
	--baseFormat 1 $otherOptions #--toolcov-cutoff 3 --min-bgtruth-cov 5  --enable-cache --using-cache #--enable-cache --using-cache # cutoff=3 is reasonable now
