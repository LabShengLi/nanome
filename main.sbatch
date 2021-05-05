#!/bin/bash
#SBATCH --job-name=nano.nxfl
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=70G # memory pool for all cores
#SBATCH --time=08:00:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


set -x

#cd $SLURM_SUBMIT_DIR
#module load singularity

# Get nextflow from internet
#curl -fsSL get.nextflow.io | bash
date;hostname;pwd


workingDir=/fastscratch/liuya/nanocompare/work
#rm -rf outputs
#rm -rf $workingDir

###############################################################
###############################################################
###############################################################
## Option 1:
# Running nanome pipeline on test data (from NA19240) in Human Genome
#rm -rf ${workingDir}-test
./nextflow run main.nf -profile winter \
    -work-dir ${workingDir}-test \
    -config conf/test.config


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
#rm -rf ${workingDir}-benchmark
#./nextflow run main.nf -profile winter  \
#	-with-report -with-timeline -with-trace -with-dag \
#	-work-dir ${workingDir}-benchmark \
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





























