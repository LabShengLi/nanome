#!/bin/bash
#SBATCH --job-name=nanome.benchmark.hpc
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=20G # memory pool for all cores
#SBATCH --time=06:00:00 # time
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

date; hostname; pwd

baseDir=${1:-/fastscratch/li-lab/nanome}

sif_dir="${baseDir}/sif"
nanome_singularity="${sif_dir}/nanome_v1.4.sif"
workDir=${baseDir}/work-benchmark
outputsDir=${baseDir}/outputs-benchmark


########################################
########################################
# Ensure directories
mkdir -p ${baseDir}; chmod ugo+w ${baseDir}
export SINGULARITY_CACHEDIR="${baseDir}/singularity-cache"
mkdir -p  $SINGULARITY_CACHEDIR; chmod ugo+w $SINGULARITY_CACHEDIR
mkdir -p $sif_dir; chmod ugo+w $sif_dir


########################################
########################################
# Get nextflow and install it
if [ ! -f "nextflow" ]; then
    curl -s https://get.nextflow.io | bash
fi


########################################
########################################
# Get nextflow and install it
if [ ! -f "nextflow" ]; then
    curl -s https://get.nextflow.io | bash
fi

########################################
########################################
# Pull nanome singularity
module load singularity
if [ ! -f ${nanome_singularity} ]; then
    singularity pull ${nanome_singularity} docker://quay.io/liuyangzzu/nanome:v1.4
fi


########################################
########################################
# Clean old results
#rm -rf ${workDir} ${outputsDir}
mkdir -p ${workDir}; chmod ugo+w ${workDir}
mkdir -p ${outputsDir}; chmod ugo+w ${outputsDir}


########################################
########################################
# Running pipeline for benchmark data
set -x
./nextflow run main.nf \
    -profile winter -resume\
    -with-report -with-timeline -with-trace -with-dag \
    -work-dir ${workDir} \
    --outputDir ${outputsDir} \
    -config conf/benchmarking.config \
    --processors 8
