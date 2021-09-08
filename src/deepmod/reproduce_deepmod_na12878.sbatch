#!/bin/bash
#SBATCH --job-name=reproduce_na12878
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=90g
#SBATCH --time=72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x
set -u
set -e

chrName=${1:-"chr20"}
predThreshold=${2:-"0.5"}

DeepModNA19240Dir='/projects/li-lab/Nanopore_compare/suppdata/deepmod-albacore/na12878_pred/ecoli_org'
DeepModChrResults="${DeepModNA19240Dir}/cpredecoli_org_clusterCpG.${chrName}.C.bed"
bgTruth="/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz;/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz"

pythonFn=reproduce_deepmod_na12878.py
python ${NanoCompareDir}/src/deepmod/${pythonFn} \
    --chr ${chrName} \
    --deepmod-input ${DeepModChrResults} \
	--bgtruth ${bgTruth} \
	--pred-threshold ${predThreshold}
