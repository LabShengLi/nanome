#!/bin/bash
#SBATCH --job-name=sanity.ecoli.metropaper
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem-per-cpu=100g
#SBATCH --time=23:59:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

callsDir=${1:-/projects/li-lab/yang/workspace/nano-compare/outputs/ECOLI_METROPAPER-methylation-callings}

DeepSignal_calls=$(find $callsDir -name "*.DeepSignal.combine.tsv.gz")
Tombo_calls=$(find $callsDir -name "*.Tombo.combine.tsv.gz")
Nanopolish_calls=$(find $callsDir -name "*.Nanopolish.combine.tsv.gz")
DeepMod_calls=$(find $callsDir -name "*.DeepModC.combine.tsv.gz")
Megalodon_calls=$(find $callsDir -name "*.Megalodon.combine.tsv.gz")

#prjBaseDir=/projects/li-lab/yang/workspace/nano-compare
pythonFn=nanocompare/read_level_eval.py

python ${pythonFn} --calls DeepSignal:${DeepSignal_calls} \
    Tombo:${Tombo_calls} Nanopolish:${Nanopolish_calls} DeepMod.C:${DeepMod_calls} \
    Megalodon:${Megalodon_calls} --runid MethPerf-ECOLI_SANITY \
    --dsname ECOLI_METROPAPER --analysis "ecoli ecoli_metropaper_sanity"
