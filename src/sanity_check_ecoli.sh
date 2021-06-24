#!/bin/bash
#SBATCH --job-name=sanity.ecoli.metropaper
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=200g
#SBATCH --time=72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

callsDir=${1:-"/projects/li-lab/Nanopore_compare/suppdata/ecoli-sanity-check/ECOLI_METROPAPER-methylation-callings"}

DeepSignal_calls=$(find $callsDir -name "*.DeepSignal.combine.tsv.gz")
Tombo_calls=$(find $callsDir -name "*.Tombo.combine.tsv.gz")
Nanopolish_calls=$(find $callsDir -name "*.Nanopolish.combine.tsv.gz")
DeepMod_calls=$(find $callsDir -name "*.DeepModC.combine.*.gz")
Megalodon_calls=$(find $callsDir -name "*.Megalodon.combine.*.gz")
Guppy_calls=$(find $callsDir -name "*.guppy.fast5mod_site_level.combine.*.gz")

pythonFn=nanocompare/read_level_eval.py

python ${pythonFn} \
	--calls \
		Nanopolish:${Nanopolish_calls} \
		Megalodon:${Megalodon_calls} \
		DeepSignal:${DeepSignal_calls} \
		Guppy:${Guppy_calls} \
    	Tombo:${Tombo_calls} \
    	DeepMod.C:${DeepMod_calls} \
	--runid MethPerf-ECOLI_SANITY \
    --dsname ECOLI_METROPAPER \
    --chrSet NC_000913.3 \
    --analysis "ecoli ecoli_metropaper_sanity"
