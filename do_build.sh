#!/bin/bash
#SBATCH --job-name=build.docker
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

## Build Docker for pipeline on google cloud

set -x

gcloud builds submit --tag us.gcr.io/jax-nanopore-01/nanome:v1.4 --timeout=4000s
