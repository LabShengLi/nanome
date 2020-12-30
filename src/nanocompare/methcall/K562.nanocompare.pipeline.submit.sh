#!/bin/bash

# set -x
###################################################################################
### Input parameter settings for each dataset                                   ###
###################################################################################
### Input parameters prepared for pipeline###
dsname=K562
ToolList=(Tombo DeepSignal)
targetNum=1
inputDataDir=/projects/li-lab/yang/results/2020-12-28/K562-Nanopore_GT18-07372.fast5.tar

### Running config, if TRUE means running this step
run_preprocessing=true
run_basecall=true
run_methcall=true
run_combine=true

### Reference file configuration
correctedGroup="RawGenomeCorrected_000"
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
chromSizesFile="/projects/li-lab/yang/results/2020-12-28/k562-tombo/hg38.chrom.sizes"

###################################################################################
###################################################################################
###################################################################################


for Tool in ${ToolList[@]}; do
    echo ${Tool}
	### Building output folder configuration
	analysisPrefix=${dsname}-${Tool}-N${targetNum}
	untaredInputDir=/fastscratch/liuya/nanocompare/${analysisPrefix}/${dsname}-untar
	septInputDir=/fastscratch/liuya/nanocompare/${analysisPrefix}/${dsname}-sept
	basecallOutputDir=/fastscratch/liuya/nanocompare/${analysisPrefix}/${dsname}-basecalled
	methCallsDir=/fastscratch/liuya/nanocompare/${analysisPrefix}/${analysisPrefix}-meth-call
	echo analysisPrefix=${analysisPrefix}

	###################################################################################
	### Preserve followings to run pipeline                                         ###
	###################################################################################
	pipelinefn=/projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.pipeline.submit.sh

	source ${pipelinefn}
	###################################################################################
	###################################################################################
	###################################################################################
done

