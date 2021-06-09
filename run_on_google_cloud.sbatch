#!/bin/bash
#SBATCH --job-name=nxf.google
#SBATCH -p compute
#SBATCH -q long
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=20G # memory pool for all cores
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
### Run Test pipeline on google cloud
## working and outputs dir
WORK_DIR_BUCKET=gs://jax-nanopore-01-project-data/TestData-work
OUTPUT_DIR_BUCKET=gs://jax-nanopore-01-export-bucket/TestData-ouputs

## Run test demo
nextflow run main.nf -profile gls -resume \
	-w ${WORK_DIR_BUCKET} \
	--input inputs/test.demo.filelist.txt \
	--outputDir ${OUTPUT_DIR_BUCKET} \
	--midDiskSize "150 GB" --highDiskSize "180 GB" \
	--runCombine false \
	-with-report -with-timeline -with-trace -with-dag

exit 0


chrName=${1:-"CHR1_3"}
inputFileListName=${2:-"inputs/na12878.google_cloud.chr1-3.filelist.txt"}

########################################################################
########################################################################
########################################################################
### Run NA12878 chr1-3 pipeline on google cloud
## working and outputs dir
WORK_DIR_BUCKET=gs://jax-nanopore-01-project-data/NA12878-${chrName}-work
OUTPUT_DIR_BUCKET=gs://jax-nanopore-01-export-bucket/NA12878-${chrName}-ouputs
resumeOptions="-resume"

## Run NA12878 datasets, only DeepMod, save time
./nextflow run main.nf -profile gls -config conf/NA12878.config ${resumeOptions} \
	-w ${WORK_DIR_BUCKET} \
	--outputDir ${OUTPUT_DIR_BUCKET} \
	--input ${inputFileListName} \
	-with-report -with-timeline -with-trace -with-dag \
	--runGuppy false --runResquiggle false --runMegalodon false --runNanopolish false
exit 0



########################################################################
########################################################################
########################################################################
### Run NA12878 pipeline on google cloud
## working and outputs dir
WORK_DIR_BUCKET=gs://jax-nanopore-01-project-data/NA12878CHR20-work
OUTPUT_DIR_BUCKET=gs://jax-nanopore-01-export-bucket/NA12878CHR20-ouputs
resumeOptions="-resume"

## Run NA12878 datasets
./nextflow run main.nf -profile gls -config conf/NA12878.config ${resumeOptions} \
	-w ${WORK_DIR_BUCKET} \
	--outputDir ${OUTPUT_DIR_BUCKET} \
	-with-report -with-timeline -with-trace -with-dag

exit 0

