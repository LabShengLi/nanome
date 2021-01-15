#!/bin/bash
#SBATCH --job-name=nanoc.pipeline.submit.APL
##SBATCH --partition=compute
#SBATCH -p gpu
#SBATCH -q inference
#SBATCH --gres=gpu:1
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

##Run on sumner: sbatch --partition=compute -q batch --gres=   APL.nanocompare.pipeline.submit.sh

# NanoCompareDir should be env var set to project dir
source ${NanoCompareDir}/src/utils.common.sh

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
#ToolList=(Guppy)

ToolList=(DeepSignal Tombo DeepMod Nanopolish)


### Which step is going to run, true or false, if 'true' means running this step

basecall_name=Guppy

run_preprocessing=false
run_basecall=false
run_resquiggling=true
run_methcall=false
run_combine=false
run_clean=false

### true if inputDataDir is a folder contains *.tar or *.tar.gz
multipleInputs=false

### which kind of intermediate file we want to backup or clean, these options are used in final stage of combine step
tar_basecall=false
tar_methcall=false
clean_basecall=false
clean_preprocessing=false

### The output base dir
outbasedir=/fastscratch/liuya/nanocompare/${dsname}-Runs
mkdir -p ${outbasedir}

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
cd ${NanoCompareDir}/src/nanocompare/methcall

source nanocompare.pipeline.submit.sh

###################################################################################
###################################################################################
###################################################################################