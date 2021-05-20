#!/bin/bash
#SBATCH --job-name=meth.perf
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=150g
#SBATCH --time=72:00:00
#SBATCH -o log/%x.%j.out
#SBATCH -e log/%x.%j.err

set -x
echo "DeepSignal_calls: $DeepSignal_calls"
echo "Tombo_calls: $Tombo_calls"
echo "Nanopolish_calls: $Nanopolish_calls"
echo "DeepMod_calls: $DeepMod_calls"
echo "Megalodon_calls: ${Megalodon_calls}"
echo "Guppy_calls: ${Guppy_calls}"
echo "bgTruth: $bgTruth"
echo "RunPrefix: $RunPrefix"
echo "parser: $parser"
echo "minCov: $minCov"
echo "dsname:" $dsname
echo "otherOptions:" $otherOptions

echo ${NanoCompareDir}

pythonFn=nanocompare/read_level_eval.py

PYTHONPATH=. python ${pythonFn} --calls DeepSignal:${DeepSignal_calls} \
	Tombo:${Tombo_calls} Nanopolish:${Nanopolish_calls} DeepMod.C:${DeepMod_calls} \
	Megalodon.ZW:${Megalodon_calls} Guppy:${Guppy_calls} --bgtruth ${parser}:${bgTruth} --runid MethPerf-${RunPrefix} \
	--dsname ${dsname} --min-bgtruth-cov ${minCov} \
	--report-joined $otherOptions ###--enable-cache --using-cache --mpi # --test
