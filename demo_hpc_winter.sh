#!/bin/bash -e
#SBATCH --job-name=nanome.demo.hpc
#SBATCH -p gpu
#SBATCH --gres=gpu:v100:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=20G # memory pool for all cores
#SBATCH --time=02:00:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR

date; hostname; pwd

# Base directory of running and output for nanome
baseDir=${1:-/fastscratch/li-lab/nanome}

workDir=${baseDir}/work
outputsDir=${baseDir}/outputs

########################################
# Clean old results
rm -rf ${workDir} ${outputsDir}

########################################
########################################
# Running pipeline for demo human data
# More options: -with-report -with-timeline -with-trace -resume -with-dag nanome_dag.png
# https://raw.githubusercontent.com/TheJacksonLaboratory/nanome/master/inputs/test.demo.filelist.txt
module load singularity
set -x
nextflow run main.nf\
    -profile singularity,hpc\
    -config conf/jax_hpc.config\
    -work-dir ${workDir}\
    --outputDir ${outputsDir}\
    --dsname TestData\
    --input https://raw.githubusercontent.com/TheJacksonLaboratory/nanome/master/inputs/test.demo.filelist.txt\
    --cleanCache false

# Report
tree ${workDir} > ${baseDir}/work_demo_filetree.txt
tree ${outputsDir} > ${baseDir}/outputs_demo_filetree.txt

echo "### nanome pipeline for demo data on HPC winter DONE"
