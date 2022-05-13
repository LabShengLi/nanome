#!/bin/bash
#SBATCH --job-name=sanity_check_deepmod
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem=190g
#SBATCH --time=72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

## Pass all following vars to py script
set -x
set -u
set -e

echo ${NanoCompareDir}

DeepModNA19240Dir='/projects/li-lab/Nanopore_compare/suppdata/deepmod-albacore/na12878_pred/ecoli_org'

Dataset="NA12878"
chrName=${1:-"chr20"}
predThreshold=${2:-"0.5"}
otherOptions=${3:-"--enable-cache --using-cache --deepmod-paper-results"}

DeepModChrResults="${DeepModNA19240Dir}/cpredecoli_org_clusterCpG.${chrName}.C.bed"
callList="DeepMod.Cluster:${DeepModChrResults}"
bgTruth="/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz;/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz"
parser="encode"

pythonFn=sanity_check_deepmod.py

PYTHONPATH=${NanoCompareDir}/src python ${pythonFn} \
    --calls ${callList} \
	--chr ${chrName} \
	--bgtruth ${parser}:${bgTruth} \
	--runid SanityCheck-${Dataset} \
	--dsname ${Dataset} \
	--processors 10 \
	--pred-threshold ${predThreshold}  ${otherOptions}
