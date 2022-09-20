#!/bin/bash
#SBATCH --job-name=sanity.ecoli.metropaper
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=300g
#SBATCH --time=72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#sbatch sanity_check_ecoli.sh

set -xe

callsDir=${1:-"/projects/li-lab/Nanopore_compare/nanome_paper_result/supp_data/ecoli-sanity-check/ECOLI_METROPAPER-methylation-callings"}

DeepSignal_calls=$(find $callsDir -maxdepth 1 -name "*.DeepSignal.combine.tsv.gz")
Tombo_calls=$(find $callsDir  -maxdepth 1 -name "*.Tombo.combine.tsv.gz")
Nanopolish_calls=$(find $callsDir  -maxdepth 1 -name "*.Nanopolish.combine.tsv.gz")
DeepMod_calls=$(find $callsDir  -maxdepth 1 -name "*.DeepModC.combine.*.gz")
Megalodon_calls=$(find $callsDir  -maxdepth 1 -name "*.Megalodon.combine.*.gz")
Guppy_calls=$(find $callsDir  -maxdepth 1 -name "*.guppy.fast5mod_site_level.combine.*.gz")
METEORE_calls=$(find $callsDir  -maxdepth 1 -name "Ecoli.METEORE.megalodon_deepsignal-optimized-model-perRead.combine.tsv.gz")

pythonFn=read_level_eval.py

${pythonFn} \
	--calls \
		Nanopolish:${Nanopolish_calls} \
		Megalodon:${Megalodon_calls} \
		DeepSignal:${DeepSignal_calls} \
		Guppy:${Guppy_calls} \
    	Tombo:${Tombo_calls} \
    	METEORE:${METEORE_calls} \
    	DeepMod:${DeepMod_calls} \
	--runid MethPerf-ECOLI_SANITY \
    --dsname ECOLI_Test \
    --chrSet NC_000913.3 \
    --analysis "ecoli ecoli_metropaper_sanity"
