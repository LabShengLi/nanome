#!/bin/bash -e
#SBATCH --job-name=nanome.demo.hpc
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 3 # number of cores
#SBATCH --mem=10G # memory pool for all cores
#SBATCH --time=01:30:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR
#SBATCH --mail-user=yang.liu@jax.org
#SBATCH --mail-type=END

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
#
module load singularity
set -ex
nextflow run main.nf\
    -profile singularity,hpc\
    -config conf/executors/jaxhpc_input.config\
    -work-dir ${workDir}\
    --outdir ${outputsDir}\
    --dsname TestData\
    --input https://raw.githubusercontent.com/TheJacksonLaboratory/nanome/master/inputs/test.demo.filelist.txt\
	--runTombo --runMETEORE --runDeepMod --useDeepModCluster\
	--outputIntermediate --outputRaw\
	--outputGenomeBrowser --outputBam --outputONTCoverage\
	--deduplicate --sort \
	--processors 8

# Report
tree ${workDir} > ${baseDir}/work_demo_filetree.txt
tree ${outputsDir} > ${baseDir}/outputs_demo_filetree.txt

echo "### nanome pipeline for demo data on HPC winter DONE"
