#!/bin/bash
#SBATCH --job-name=nanome.google_tera
#SBATCH -p compute
#SBATCH -q long
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=16G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR
#SBATCH --mail-user=yang.liu@jax.org
#SBATCH --mail-type=END

set -e
date;hostname;pwd
chrName=${1:-"chr22"}
baseDir=${2:-"/fastscratch/$USER/nanome"}

pipelineDir=${baseDir}/na12878_${chrName}_gcptest
rm -rf $pipelineDir
mkdir -p $pipelineDir
cd $pipelineDir

###########################################
###########################################
###########################################
### Run pipeline on google cloud
## working and outputs dir
WORK_DIR_BUCKET=${1:-"gs://jax-nanopore-01-project-data/NANOME_tera-work"}
OUTPUT_DIR_BUCKET=${2:-"gs://jax-nanopore-01-export-bucket/NANOME_tera_ouputs"}

gsutil -m rm -rf ${WORK_DIR_BUCKET}  ${OUTPUT_DIR_BUCKET} >/dev/null 2>&1 || true

###########################################
###########################################
###########################################
### Run pipeline on google for NA12878
set -x
echo "### nanome pipeline for NA12878 some chr and part file on google START"
nextflow run ${NANOME_DIR}\
    -profile docker,google -resume -with-report -with-timeline -with-trace -with-dag\
	-w ${WORK_DIR_BUCKET} \
	--outdir ${OUTPUT_DIR_BUCKET}\
	--dsname NA12878_${chrName^^}\
	--input ${NANOME_DIR}/inputs/na12878_${chrName}_gs.filelist.txt\
	--cleanAnalyses true\
	--tomboResquiggleOptions '--signal-length-range 0 500000  --sequence-length-range 0 50000'\
	--midDiskSize "850.GB" --highDiskSize "1024.GB"\
	--machineType n1-highmem-16

echo "### nanome pipeline for NA12878 some chr and part file on google DONE"
