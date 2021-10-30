#!/bin/bash -e
#SBATCH --job-name=nanome.google.tera
#SBATCH -p compute
#SBATCH -q long
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=6G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR
#SBATCH --mail-user=yang.liu@jax.org
#SBATCH --mail-type=END

set -ex
date;hostname;pwd

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
echo "### nanome pipeline for NA12878 some chr and part file on google START"
nextflow run main.nf\
    -profile docker,google -resume\
    -config conf/executors/gcp_input.config\
	-w ${WORK_DIR_BUCKET} \
	--outdir ${OUTPUT_DIR_BUCKET} \
	--dsname NA12878_CHR22_P3 \
	--input 'gs://jax-nanopore-01-project-data/na12878/s3.amazonaws.com/nanopore-human-wgs/rel3-fast5-chr22.part03.tar'\
	--cleanAnalyses true\
	--tomboResquiggleOptions '--signal-length-range 0 500000  --sequence-length-range 0 50000'\
	--reduceProcTimes  0.5

echo "### nanome pipeline for NA12878 some chr and part file on google DONE"
