#!/bin/bash
#SBATCH --job-name=nanoc.pipeline.submit.APL
#SBATCH --partition=compute
#SBATCH --mem=50g # memory pool for all cores
#SBATCH --time=03:00:00 # time (DD-HH:MM:SS)
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

###################################################################################
### Settings for APL dataset experimentation                                   ###
### This is the only file we need to modify for different data                  ###
### Nanopore tools tested for DeepSignal, Nanopolish, DeepMod and Tombo         ###
### Nano-compare project by Yang Liu                                            ###
###################################################################################
# set -x

### Input dataset parameters prepared for pipeline###
# dsname    -   data set name
# targetNum -   number of tasks splited to run
dsname=APL
targetNum=50

### Nanopore raw signal fast5 files, K562 is from: /projects/li-lab/AML-Nanopore/20180517_180508-18-li-001-GXB01186-001/APL-1750_GT18-06409.fast5.tar
# inputDataDir  -   input of Nanopore reads file/files, can be tar file or a directory contains files
inputDataDir=/fastscratch/liuya/nanocompare/Nanopore-reads/APL/APL-1750_GT18-06409.fast5.tar

### Running configurations
### which nanopore tools can be used, such as ToolList=(Tombo DeepSignal)
# ToolList  -   a list of Nanopore tools prepared to run
#ToolList=(DeepMod)

ToolList=(DeepSignal Tombo DeepMod Nanopolish)


### Which step is going to run, true or false, if 'true' means running this step

basecall_name=Albacore

run_preprocessing=true
run_basecall=true
run_resquiggling=true
run_methcall=true
run_combine=true
run_clean=false

### true if inputDataDir is a folder contains *.tar or *.tar.gz
multipleInputs=false

### which kind of intermediate file we want to backup or clean, these options are used in final stage of combine step
tar_basecall=false
tar_methcall=false
clean_basecall=true
clean_preprocessing=false

### The output base dir
outbasedir=/fastscratch/liuya/nanocompare/${dsname}-Runs
mkdir -p ${outbasedir}

### Reference file path configuration, used by each base or meth calling
correctedGroup="RawGenomeCorrected_000"
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
chromSizesFile="/projects/li-lab/yang/workspace/nano-compare/data/genome-annotation/hg38.chrom.sizes"
deepsignalModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/model.CpG.R9.4_1D.human_hx1.bn17.sn360/bn_17.sn_360.epoch_7.ckpt"
deepModModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/rnn_conmodC_P100wd21_f7ne1u0_4/mod_train_conmodC_P100wd21_f3ne1u0"
clusterDeepModModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/na12878_cluster_train_mod-keep_prob0.7-nb25-chr1/Cg.cov5.nb25"

### Number of processes for basecall, alignment, and methlation nanopore tool
processors=16
#processors=64

isGPU="no"
###################################################################################
###################################################################################
###################################################################################


###################################################################################
### Preserve followings to run Base Modified Prediction pipeline                ###
###################################################################################
# Please put this file at nano-compare/src dir, or it need modify following paths
# change working path to script path
cd "$(dirname "$0")"/nanocompare/methcall

source nanocompare.pipeline.submit.sh

###################################################################################
###################################################################################
###################################################################################