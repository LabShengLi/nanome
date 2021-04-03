#!/bin/bash
#SBATCH --job-name=nanoc.pipeline.submit.NA19240
#SBATCH --partition=compute
##SBATCH -p gpu
##SBATCH -q inference
##SBATCH --gres=gpu:1
#SBATCH --mem=50g # memory pool for all cores
#SBATCH --time=03:00:00 # time (DD-HH:MM:SS)
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err

###################################################################################
### Settings for NA19240 dataset experimentation                                ###
### This is the only file we need to modify for different data                  ###
### Nanopore tools tested for DeepSignal, Nanopolish, DeepMod and Tombo         ###
### Nano-compare project by Yang Liu                                            ###
###################################################################################

## Sumner run: sbatch --partition=compute -q batch --gres= NA19240.nanocompare.pipeline.submit.sh

### Input dataset parameters prepared for pipeline###
# dsname    -   data set name
# targetNum -   number of tasks splited to run
dsname=NA19240
targetNum=300 # due to too large number of fast5 files

# inputDataDir  -   input of Nanopore reads file/files, can be tar file or a directory contains files
inputDataDir=/fastscratch/liuya/nanocompare/input/NA19240

### Running configurations
### which nanopore tools can be used, such as ToolList=(Tombo DeepSignal)
# ToolList  -   a list of Nanopore tools prepared to run
ToolList=(DeepMod)
#ToolList=(DeepSignal Tombo)


### Which step is going to run, true or false, if 'true' means running this step
basecall_name=Albacore
#basecall_name=Guppy

### true if inputDataDir is a folder contains *.tar or *.tar.gz
multipleInputs=true

### The output base dir
outbasedir=/fastscratch/liuya/nanocompare/${dsname}-Runs
mkdir -p ${outbasedir}

### Number of processes for basecall, alignment, and methlation nanopore tool
processors=16
#processors=64

isGPU="no"
isGPU="yes"
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