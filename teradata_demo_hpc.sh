#!/bin/bash
#SBATCH --job-name=nanome.teradata_hpc
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 12 # number of cores
#SBATCH --mem=10G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR
#SBATCH --mail-user=yang.liu@jax.org
#SBATCH --mail-type=END

# sbatch teradata_demo_hpc.sh chr22
set -e
date; hostname; pwd
chrName=${1:-"chr22"}
baseDir=${2:-"/fastscratch/$USER/nanome"}

pipelineDir=${baseDir}/na12878_${chrName}_test
rm -rf $pipelineDir
mkdir -p $pipelineDir
cd $pipelineDir

########################################
########################################
# Running pipeline for tera data, more options:
module load singularity
echo "### Start test on teradata ${chrName}"

set -x
nextflow run ${NANOME_DIR}\
        -resume -with-report -with-timeline -with-trace -with-dag\
        -profile singularity,hpc \
        -config ${NANOME_DIR}/conf/executors/jaxhpc_input.config,${NANOME_DIR}/conf/executors/na12878_hpc.config\
        --dsname NA12878_${chrName^^} \
        --input "https://raw.githubusercontent.com/TheJacksonLaboratory/nanome/master/inputs/na12878_${chrName}.filelist.txt"
echo "### Run pass on teradata ${chrName} for NANOME"
