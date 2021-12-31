#!/bin/bash
#SBATCH --job-name=nanome.ecoli_demo_hpc
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=12G # memory pool for all cores
#SBATCH --time=02:00:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR
#SBATCH --mail-user=yang.liu@jax.org
#SBATCH --mail-type=END
set -e
date; hostname; pwd

# Base directory of running and output for nanome
baseDir=${1:-/fastscratch/$USER/nanome}
pipelineName=${2:-"ecoli_demo"}

rm -rf ${baseDir}/${pipelineName}
mkdir -p ${baseDir}/${pipelineName}

########################################
########################################
# Running pipeline for E. coli data
module load singularity
set -x

cd ${baseDir}/${pipelineName}

nextflow run ${NANOME_DIR}\
    -resume -with-report -with-timeline -with-trace -with-dag\
    -profile singularity,hpc\
    --dsname EcoliDemo\
    -config ${NANOME_DIR}/conf/executors/jaxhpc_input.config,${NANOME_DIR}/conf/examples/ecoli_demo.config\
    --runDeepMod --runTombo --runMETEORE\
    --outputIntermediate --outputRaw\
    --outputGenomeBrowser --outputBam --outputONTCoverage\
    --deduplicate --sort\
    --processors 8

# Report
tree work >  ${pipelineName}_work_filetree.txt
tree results >  ${pipelineName}_results_filetree.txt

# save all session records
tar -czf ${pipelineName}.tar.gz  \
    ${pipelineName}_results_filetree.txt ${pipelineName}_work_filetree.txt \
    work/*/*/.command.* work/*/*/*.run.log \
    *trace/  .nextflow.log

find  work  -name '*.Report.run.log' -exec tail {} \;

echo "### NANOME pipeline for ecoli_demo on HPC DONE"
