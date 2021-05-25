#!/bin/bash
#SBATCH --job-name=nxf.google
#SBATCH -p compute
#SBATCH -q long
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=50G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

date;hostname;pwd

# Get nextflow and install
# curl -s https://get.nextflow.io | bash

########################################################################
########################################################################
########################################################################
### Run NA12878 pipeline on google cloud
## working and outputs dir
WORK_DIR_BUCKET=gs://jax-nanopore-01-project-data/NA12878-work
OUTPUT_DIR_BUCKET=gs://jax-nanopore-01-export-bucket/NA12878-ouputs

## Run NA12878 datasets
./nextflow run main.nf -profile gls -config conf/NA12878.config -resume \
	-w ${WORK_DIR_BUCKET} \
	--outputDir ${OUTPUT_DIR_BUCKET}

exit 0

########################################################################
########################################################################
########################################################################
### Run Test pipeline on google cloud
## working and outputs dir
WORK_DIR_BUCKET=gs://jax-nanopore-01-project-data/TestData-work
OUTPUT_DIR_BUCKET=gs://jax-nanopore-01-export-bucket/TestData-ouputs

## Run test demo
./nextflow run main.nf -profile gls -resume \
	-w ${WORK_DIR_BUCKET} \
	--input inputs/test.demo.filelist.txt \
	--outputDir ${OUTPUT_DIR_BUCKET} \
	--midDiskSize "50 GB" --highDiskSize "50 GB"

exit 0
