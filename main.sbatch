#!/bin/bash
#SBATCH --job-name=nanoc.pipe
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem-per-cpu=170G # memory pool for all cores
#SBATCH --time=00:50:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#cd $SLURM_SUBMIT_DIR
#date;hostname;pwd
#
#module load singularity

#curl -fsSL get.nextflow.io | bash
set -x

indata=/fastscratch/liuya/nanocompare/input/HL60-Nanopore_GT18-07373.fast5.tar
ntask=50
workdir=/fastscratch/liuya/nanocompare/work-HL60

#rm -rf ${workdir}

#./nextflow run main.nf -profile winter

#/nextflow run main.nf -resume -with-report -with-timeline -with-trace
#
#./nextflow run main.nf -work-dir ${workdir} \
#	--indata ${indata} --ntask ${ntask} \
#	-cache -resume -with-report -with-timeline \
#	-with-trace