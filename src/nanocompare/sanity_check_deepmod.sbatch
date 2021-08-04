#!/bin/bash
#SBATCH --job-name=sanity_check
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

Dataset=${1:-"NA12878-CHR20"}
chrName=${2:-"chr20"}
callList=${3:-"/projects/li-lab/Nanopore_compare/data/NA12878/CHR20/NA12878CHR20.deepmod.C.combine.bed.gz"}
bgTruth=${4:-"/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz;/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz"}
parser=${5:-"encode"}
predThreshold=${6:-"0.5"}
otherOptions=${7:-"--enable-cache --using-cache"}

pythonFn=sanity_check_deepmod.py

PYTHONPATH=${NanoCompareDir}/src python ${pythonFn} --calls ${callList} \
	--chr ${chrName} --bgtruth ${parser}:${bgTruth} \
	--runid SanityCheck-${Dataset} --dsname ${Dataset} \
	--baseFormat 1 ${otherOptions} --processors 10 --pred-threshold ${predThreshold} #--toolcov-cutoff 3 --min-bgtruth-cov 5  --enable-cache --using-cache #--enable-cache --using-cache # cutoff=3 is reasonable now
