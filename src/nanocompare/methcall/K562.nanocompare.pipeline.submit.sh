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
targetNum=50
inputDataDir=/projects/li-lab/yang/workspace/nano-compare/data/raw-fast5/K562/K562-Nanopore_GT18-07372.fast5.tar

### Running configurations
### which nanopore tools can be used, such as ToolList=(Tombo DeepSignal)
ToolList=(DeepSignal)
#ToolList=(Tombo)

### Which step is going to run, true or false, if 'true' means running this step
#run_preprocessing=false
#run_basecall=false
#run_resquiggling=false
#run_methcall=false
#run_combine=true

run_preprocessing=true
run_basecall=true
run_resquiggling=true
run_methcall=true
run_combine=true

### Reference file path configuration, used by each base or meth calling
correctedGroup="RawGenomeCorrected_000"
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
chromSizesFile="/projects/li-lab/yang/workspace/nano-compare/data/genome-annotation/hg38.chrom.sizes"

deepsignalModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/model.CpG.R9.4_1D.human_hx1.bn17.sn360/bn_17.sn_360.epoch_7.ckpt"
isGPU="no"

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