#!/bin/bash
#SBATCH --job-name=nanome.ecoli
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=20G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

date; hostname; pwd

baseDir=${1:-/fastscratch/li-lab/nanome}

workDir=${baseDir}/work-ecoli
outputsDir=${baseDir}/outputs-ecoli

########################################
########################################
# Ensure directories
mkdir -p ${baseDir}; chmod ugo+w ${baseDir}
export SINGULARITY_CACHEDIR="${baseDir}/singularity-cache"
mkdir -p  $SINGULARITY_CACHEDIR; chmod ugo+w $SINGULARITY_CACHEDIR


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
mkdir -p ${workDir}; chmod ugo+w ${workDir}
mkdir -p ${outputsDir}; chmod ugo+w ${outputsDir}


########################################
########################################
# Running pipeline for E. coli data
module load singularity

set -x
./nextflow run main.nf \
    -profile singularity,hpc \
    -work-dir ${workDir} \
    --outputDir ${outputsDir} \
    -config conf/jax_hpc.config,conf/ecoli_demo.config \
    --cleanCache false \
    --singularity_cache_dir '/fastscratch/li-lab/nanome/singularity-cache'

# Report
tree ${workDir} > ${baseDir}/work_ecoli.tree.txt
tree ${outputsDir} > ${baseDir}/outputs_ecoli.tree.txt

echo "### nanome pipeline for ecoli data DONE"
