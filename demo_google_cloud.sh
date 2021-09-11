#!/bin/bash
#SBATCH --job-name=nanome.google
#SBATCH -p compute
#SBATCH -q long
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=30G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

date;hostname;pwd

###########################################
###########################################
###########################################
### Run Test pipeline on google cloud
## working and outputs dir
WORK_DIR_BUCKET=gs://jax-nanopore-01-project-data/TestData-work
OUTPUT_DIR_BUCKET=gs://jax-nanopore-01-export-bucket/TestData-ouputs

gsutil rm -rf ${WORK_DIR_BUCKET}
gsutil rm -rf ${OUTPUT_DIR_BUCKET}

## Run test demo on google cloud
nextflow run main.nf -resume\
    -profile docker,google \
	-w ${WORK_DIR_BUCKET} \
	--outputDir ${OUTPUT_DIR_BUCKET} \
	--dsname TestData \
	--input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt
