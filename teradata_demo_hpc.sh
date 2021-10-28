#!/bin/bash -e
#SBATCH --job-name=nanome.teradata.demo
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=10G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR
#SBATCH --mail-user=yang.liu@jax.org
#SBATCH --mail-type=END
date; hostname; pwd
# run on winter: chr22
#             sbatch teradata_demo_hpc.sh

# Base directory of running and output for nanome
baseDir=${1:-/fastscratch/li-lab/nanome}
workDir=${baseDir}/work-na12878
outputsDir=${baseDir}/outputs-na12878

########################################
# Clean old results
rm -rf ${workDir} ${outputsDir} &&\
    mkdir -p ${workDir} ${outputsDir}

########################################
########################################
# Running pipeline for tera data, more options:  -with-report -with-timeline -with-trace -with-dag
module load singularity
echo "### Start test on teradata"
set -ex
nextflow run main.nf -resume\
        -profile singularity,hpc \
        -config conf/executors/jaxhpc_input.config,conf/executors/na12878_hpc.config\
        -work-dir ${workDir} \
        --outdir ${outputsDir} \
        --dsname NA12878_CHR22 \
        --input 'https://raw.githubusercontent.com/TheJacksonLaboratory/nanome/master/inputs/na12878_chr22.filelist.txt'

echo "### Run pass on teradata for NANOME"
