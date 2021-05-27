#!/bin/bash
#SBATCH --job-name=nano.nxfl
#SBATCH -p compute
#SBATCH -q long
##SBATCH -p gpu
##SBATCH --gres=gpu:1
##SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=70G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

date;hostname;pwd

# Get nextflow and install it
# curl -s https://get.nextflow.io | bash

#workingDir=/fastscratch/liuya/nanocompare/work
#rm -rf outputs
#rm -rf $workingDir

###############################################################
###############################################################
###############################################################
## Option 1:
# Running nanome pipeline on test data (from NA19240) in Human Genome
#rm -rf ${workingDir}-test4
nextflow run main.nf -profile sumner --input inputs/test.demo.filelist.txt \
    -work-dir ${workingDir}-test6
#    -config conf/test.config


#./nextflow run main.nf -profile winter \
#    -work-dir ${workingDir}-test3 \
#    -config conf/NA12878.config

###############################################################
###############################################################
###############################################################
## Option 2:
## Running nanome pipeline on E. coli data
#./nextflow run main.nf -profile winter \
#    -work-dir ${workingDir}-ecoli \
#    -config conf/ecoli_demo.config -resume


###############################################################
###############################################################
###############################################################
## Option 3:
## Run on benchmarking data
#rm -rf ${workingDir}-benchmark3
#./nextflow run main.nf -profile winter  \
#	-with-report -with-timeline -with-trace -with-dag \
#	-work-dir ${workingDir}-benchmark3 \
#	-resume -config conf/benchmarking.config

###############################################################
###############################################################
###############################################################
## Option 4:
## Run on APL
#./nextflow run main.nf -profile winter --input test/APL.filelist.txt -resume --runMethcall false

## Run only basecall on APL
#./nextflow run main.nf -profile winter --input test/APL.filelist.txt -resume --runMethcall false

###############################################################
###############################################################
###############################################################
#./nextflow run main.nf -profile winter,test -with-singularity nanome.sif


# module load singularity
# builder="singularity run http://s3-far.jax.org/builder/builder"
# $builder nanocompare.def nanocompare.sif

