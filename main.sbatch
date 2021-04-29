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

date;hostname;pwd

#module load singularity

# Get nextflow from internet
#curl -fsSL get.nextflow.io | bash


# Test default test datasets on winter HPC
#./nextflow run main.nf -profile winter,test --ntask 1

#./nextflow run main.nf -profile winter,test -resume -with-dag -bg

#./nextflow run main.nf -profile winter,test -with-report -with-timeline -with-trace  -with-dag

#./nextflow run main.nf -profile winter,test -with-singularity nanome.sif


# Run on benchmarking data
#rm -rf /fastscratch/liuya/nanocompare/work
#./nextflow run main.nf -profile winter,benchmarking  \
#	-with-report -with-timeline -with-trace -with-dag \
#	-work-dir /fastscratch/liuya/nanocompare/work2
#


#./nextflow run main.nf -profile winter \
#  --input_fast5_tar http://s3.amazonaws.com/nanopore-human-wgs/rel3-fast5-chr20.part05.tar \
#  --online true \
#  -work-dir /fastscratch/liuya/nanocompare/work \
#  -resume


rm -rf /fastscratch/liuya/nanocompare/work

./nextflow run main.nf -profile winter \
    -work-dir /fastscratch/liuya/nanocompare/work \
    -config conf/Ecoli_METROPAPER.config -resume

#./nextflow run main.nf -profile winter \
#    -work-dir /fastscratch/liuya/nanocompare/work \
#    --input /fastscratch/c-panz/ecoli_meteore.tar.gz \
#    --dsname ECOLI_METROPAPER \
#    --referenceGenome reference_genome/ecoli/Ecoli_k12_mg1655.fasta \
#    --chromSizesFile reference_genome/ecoli/Ecoli_k12_mg1655.fasta.genome.sizes

#./nextflow run main.nf -profile winter \
#    -work-dir /fastscratch/liuya/nanocompare/work -resume \
#    --input https://github.com/liuyangzzu/nano-compare/raw/6684b71587c19d191ccf87832f25536be10e372f/test_data/demo.fast5.reads.tar.gz \
#    --dsname NA12878 --online true




#./nextflow run main.nf -profile winter,benchmarking

# Test with output report timelines
#/nextflow run main.nf -profile winter -with-report -with-timeline -with-trace

# Running resource test
#./nextflow main.nf run --eval false -profile winter \
#   -with-report -with-timeline -with-trace -work-dir /fastscratch/liuya/nanocompare/work

# Test HL60 as input, $ntask subjobs on winter, output to $workdir
#indata=/fastscratch/liuya/nanocompare/input/HL60-Nanopore_GT18-07373.fast5.tar
#ntask=50
#workdir=/fastscratch/liuya/nanocompare/work-HL60

#./nextflow run main.nf -work-dir ${workdir} \
#	--indata ${indata} --ntask ${ntask} \
#	-cache -resume -with-report -with-timeline \
#	-with-trace

# module load singularity
# builder="singularity run http://s3-far.jax.org/builder/builder"
# $builder nanocompare.def nanocompare.sif





























