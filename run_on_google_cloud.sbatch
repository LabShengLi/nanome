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


## Run NA12878 datasets
./nextflow run main.nf -profile gls -config conf/NA12878.config -resume \
	-w gs://jax-nanopore-01-project-data/nanome-work-NA12878 \
	--outputDir gs://jax-nanopore-01-project-data/nanome-outputs-NA12878

exit 0

## Run test demo
./nextflow run main.nf -profile gls -resume \
	-w gs://jax-nanopore-01-project-data/nanome-work-test \
	--outputDir gs://jax-nanopore-01-project-data/nanome-outputs

exit 0
