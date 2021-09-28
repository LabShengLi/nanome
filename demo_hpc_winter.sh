#!/bin/bash
#SBATCH --job-name=nanome.demo.hpc
#SBATCH -p gpu
#SBATCH --gres=gpu:v100:1
#SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=20G # memory pool for all cores
#SBATCH --time=1-00:00:00 # time
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

date; hostname; pwd

# Base directory of running and output for nanome
baseDir=${1:-/fastscratch/li-lab/nanome}

workDir=${baseDir}/work
outputsDir=${baseDir}/outputs


########################################
########################################
# Get nextflow and install it
if [ ! -f "nextflow" ]; then
    curl -s https://get.nextflow.io | bash
fi


########################################
# Clean old results
rm -rf ${workDir} ${outputsDir}


########################################
########################################
# Running pipeline for demo human data
# More options: -with-report -with-timeline -with-trace -with-dag -resume
# https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt
module load singularity
set -x
./nextflow run main.nf -resume \
    -with-dag nanome_dag.png \
    -profile singularity,hpc \
    -config conf/jax_hpc.config \
    -work-dir ${workDir} \
    --outputDir ${outputsDir} \
    --singularity_cache "${baseDir}/singularity-cache" \
    --dsname TestData \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt \
    --cleanCache false

# Report
tree ${workDir} > ${baseDir}/work_demo_filetree.txt
tree ${outputsDir} > ${baseDir}/outputs_demo_filetree.txt

echo "### nanome pipeline for demo data on HPC winter DONE"
