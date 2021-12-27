#!/bin/bash
#SBATCH --job-name=nanome.google_demo
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=6G # memory pool for all cores
#SBATCH --time=3-00:00:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR
#SBATCH --mail-user=yang.liu@jax.org
#SBATCH --mail-type=END
set -e
baseDir=${1:-/fastscratch/$USER/nanome}
gcpProjectName=${2:-"jax-nanopore-01"}
WORK_DIR_BUCKET=${3:-"gs://jax-nanopore-01-project-data/NANOME-TestData-work"}
OUTPUT_DIR_BUCKET=${4:-"gs://jax-nanopore-01-export-bucket/NANOME-TestData-ouputs"}

pipelineName="gcp_nanome_demo"

rm -rf $baseDir/$pipelineName
mkdir -p $baseDir/$pipelineName
cd $baseDir/$pipelineName

set -ex
date;hostname;pwd

###########################################
###########################################
###########################################
### Run Test pipeline on google cloud
set +x
source $(conda info --base)/etc/profile.d/conda.sh
conda activate py39
gsutil -m rm -rf ${WORK_DIR_BUCKET}  ${OUTPUT_DIR_BUCKET} >/dev/null 2>&1 || true
set -x

## Run test demo on google cloud
echo "### NANOME pipeline for demo data on google START"
nextflow run ${NANOME_DIR}\
    -profile docker,google \
	-w ${WORK_DIR_BUCKET} \
	--outdir ${OUTPUT_DIR_BUCKET} \
	--googleProjectName ${gcpProjectName}\
	--dsname TestData \
	--input https://storage.googleapis.com/jax-nanopore-01-project-data/nanome-input/demo1_fast5_reads.tar.gz

echo "### NANOME pipeline for demo data on google DONE"
