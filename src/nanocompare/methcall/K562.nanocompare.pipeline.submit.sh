#!/bin/bash
#SBATCH --job-name=nanoc.pipeline.submit
#SBATCH --partition=compute
#SBATCH --mem=50g # memory pool for all cores
#SBATCH --time=03:00:00 # time (DD-HH:MM:SS)
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

###################################################################################
### Settings for K562 dataset experimentation                                   ###
### This is the only file we need to modify for different data                  ###
### Nanopore tools tested for DeepSignal, Nanopolish, DeepMod and Tombo         ###
### Nano-compare project by Yang Liu                                            ###
###################################################################################
# set -x

### Input parameters prepared for pipeline###
dsname=K562
targetNum=1
inputDataDir=/projects/li-lab/yang/workspace/nano-compare/data/raw-fast5/K562/K562-Nanopore_GT18-07372.fast5.tar

### Running configurations
### which nanopore tools can be used, such as ToolList=(Tombo DeepSignal)

ToolList=(Tombo)

#ToolList=(DeepSignal Tombo DeepMod Nanopolish)


### Which step is going to run, true or false, if 'true' means running this step

run_preprocessing=false
run_basecall=true
run_resquiggling=true
run_methcall=true
run_combine=true
run_clean=false

# true if inputDataDir is a folder contains *.tar or *.tar.gz
multipleInputs=false

# which kind of intermediate file we want to clean
clean_preprocessing=true
clean_basecall=false

# The output base dir
outbasedir=/fastscratch/liuya/nanocompare/${dsname}-Runs
mkdir -p ${outbasedir}

### Reference file path configuration, used by each base or meth calling

correctedGroup="RawGenomeCorrected_000"
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
chromSizesFile="/projects/li-lab/yang/workspace/nano-compare/data/genome-annotation/hg38.chrom.sizes"

deepsignalModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/model.CpG.R9.4_1D.human_hx1.bn17.sn360/bn_17.sn_360.epoch_7.ckpt"
deepModModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/rnn_conmodC_P100wd21_f7ne1u0_4/mod_train_conmodC_P100wd21_f3ne1u0"
clusterDeepModModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/na12878_cluster_train_mod-keep_prob0.7-nb25-chr1/Cg.cov5.nb25"

isGPU="no"

# Number of processes for basecall, alignment, and methlation nanopore tool
processors=64

###################################################################################
###################################################################################
###################################################################################

###################################################################################
### Preserve followings to run pipeline                                         ###
###################################################################################

source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.pipeline.submit.sh

###################################################################################
###################################################################################
###################################################################################