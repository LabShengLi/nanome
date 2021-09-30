#!/bin/bash
#SBATCH --job-name=nanome.google.demo
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=30G # memory pool for all cores
#SBATCH --time=72:00:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR

set -x

date;hostname;pwd

###########################################
###########################################
###########################################
### Run Test pipeline on google cloud
## working and outputs dir
WORK_DIR_BUCKET=${1:-"gs://jax-nanopore-01-project-data/NANOME-TestData-work"}
OUTPUT_DIR_BUCKET=${2:-"gs://jax-nanopore-01-export-bucket/NANOME-TestData-ouputs"}

gsutil rm -rf ${WORK_DIR_BUCKET}  ${OUTPUT_DIR_BUCKET}

## Run test demo on google cloud
nextflow run main.nf -resume\
    -profile docker,google \
	-w ${WORK_DIR_BUCKET} \
	--outputDir ${OUTPUT_DIR_BUCKET} \
	--dsname TestData \
	--input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt

echo "### nanome pipeline for demo data on google DONE"
