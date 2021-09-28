#!/bin/bash
#SBATCH --job-name=nanome.ecoli.hpc
#SBATCH -p gpu
#SBATCH --gres=gpu:v100:1
#SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=20G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

date; hostname; pwd

# Base directory of running and output for nanome
baseDir=${1:-/fastscratch/li-lab/nanome}

workDir=${baseDir}/work-ecoli
outputsDir=${baseDir}/outputs-ecoli

########################################
########################################
# Get nextflow and install it
if [ ! -f "nextflow" ]; then
    curl -s https://get.nextflow.io | bash
fi

########################################
########################################
# Clean old results
rm -rf ${workDir} ${outputsDir}

########################################
########################################
# Running pipeline for E. coli data
module load singularity
set -x
./nextflow run main.nf -resume\
    -profile singularity,hpc \
    -work-dir ${workDir} \
    --outputDir ${outputsDir} \
    -config conf/jax_hpc.config,conf/ecoli_demo.config \
    --singularity_cache "${baseDir}/singularity-cache" \
    --cleanCache false

# Report
tree ${workDir} > ${baseDir}/work_ecoli_filetree.txt
tree ${outputsDir} > ${baseDir}/outputs_ecoli_filetree.txt

echo "### nanome pipeline for ecoli data on HPC DONE"
