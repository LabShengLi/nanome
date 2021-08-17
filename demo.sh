#!/bin/bash
#SBATCH --job-name=nanome.demo
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=20G # memory pool for all cores
#SBATCH --time=06:00:00 # time
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

outDir=${1:-/fastscratch/li-lab/nanome}
mkdir -p $outDir

date; hostname; pwd

# Get nextflow and install it
curl -s https://get.nextflow.io | bash

# Pull nanome singularity
module load singularity

export SINGULARITY_CACHEDIR="/fastscratch/yang/singularity-cache"
mkdir -p  $SINGULARITY_CACHEDIR

sif_dir=/fastscratch/yang/sif
mkdir -p $sif_dir
cd $sif_dir
singularity pull docker://quay.io/liuyangzzu/nanome:v1.4
nanome_singularity="${sif_dir}/nanome_v1.4.sif"

cd -
./nextflow run main.nf \
    -profile winter2 -resume \
    -with-report -with-timeline -with-trace -with-dag \
    -work-dir ${outDir}/work \
    -with-singularity ${nanome_singularity} \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt \
    --outputDir ${outDir}/outputs
