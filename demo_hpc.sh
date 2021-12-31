#!/bin/bash
#SBATCH --job-name=nanome.human_demo_hpc
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=10G # memory pool for all cores
#SBATCH --time=01:30:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR
#SBATCH --mail-user=yang.liu@jax.org
#SBATCH --mail-type=END
set -e
date; hostname; pwd

# Base directory of running and output for nanome
baseDir=${1:-/fastscratch/$USER/nanome}
pipelineName=${2:-"human_demo"}

rm -rf ${baseDir}/${pipelineName}
mkdir -p ${baseDir}/${pipelineName}

########################################
########################################
# Running pipeline for demo human data
# More options: -with-report -with-timeline -with-trace -resume -with-dag nanome_dag.png
#
module load singularity
set -x

cd ${baseDir}/${pipelineName}
nextflow run ${NANOME_DIR}\
    -resume -with-report -with-timeline -with-trace -with-dag\
    -profile singularity,hpc\
    -config ${NANOME_DIR}/conf/executors/jaxhpc_input.config\
    --dsname TestData\
    --input https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/demo1_fast5_reads.tar.gz\
	--runTombo --runMETEORE --runDeepMod --useDeepModCluster\
	--outputIntermediate --outputRaw\
	--outputGenomeBrowser --outputBam --outputONTCoverage\
	--deduplicate --sort \
	--processors 8

# Report
tree work > ${pipelineName}_work_filetree.txt
tree results >  ${pipelineName}_results_filetree.txt

# save all session records
tar -czf ${pipelineName}.tar.gz  \
    ${pipelineName}_results_filetree.txt ${pipelineName}_work_filetree.txt \
    work/*/*/.command.* work/*/*/*.run.log \
    *trace/  .nextflow.log

find  work  -name '*.Report.run.log' -exec tail {} \;

echo "### NANOME pipeline for human_demo data on HPC winter DONE"
