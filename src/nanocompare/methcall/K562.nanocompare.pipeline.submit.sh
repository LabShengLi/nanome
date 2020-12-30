#!/bin/bash
#SBATCH --job-name=nano.pipeline.submit
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 250g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

###################################################################################
### Input parameter settings for K562 dataset                                   ###
###################################################################################
# set -x

### Input parameters prepared for pipeline###
dsname=K562
targetNum=60
inputDataDir=/projects/li-lab/yang/workspace/nano-compare/data/raw-fast5/K562/K562-Nanopore_GT18-07372.fast5.tar

### which nanopore tools can be used, such as ToolList=(Tombo DeepSignal)
ToolList=(Tombo)

### Running config, if 'true' means running this step
### true or false
run_preprocessing=true
run_basecall=true
run_resquiggling=true
run_methcall=true
run_combine=true

### Reference file path configuration, used by each base or meth calling
correctedGroup="RawGenomeCorrected_000"
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
chromSizesFile="/projects/li-lab/yang/workspace/nano-compare/data/genome-annotation/hg38.chrom.sizes"

###################################################################################
###################################################################################
###################################################################################


###################################################################################
### Preserve followings to run pipeline                                         ###
###################################################################################
for Tool in ${ToolList[@]}; do
	### Building output folder configuration
	analysisPrefix=${dsname}-${Tool}-N${targetNum}
	untaredInputDir=/fastscratch/liuya/nanocompare/${analysisPrefix}/${analysisPrefix}-untar
	septInputDir=/fastscratch/liuya/nanocompare/${analysisPrefix}/${analysisPrefix}-sept
	basecallOutputDir=/fastscratch/liuya/nanocompare/${analysisPrefix}/${analysisPrefix}-basecall
	methCallsDir=/fastscratch/liuya/nanocompare/${analysisPrefix}/${analysisPrefix}-meth-call

	### Start script for submiting jobs
	echo "Start pipeline submit for analysisPrefix=${analysisPrefix}"
	pipelinefn=/projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.pipeline.submit.sh
	source ${pipelinefn}

done
###################################################################################
###################################################################################
###################################################################################