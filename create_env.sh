#!/bin/bash
#SBATCH --job-name=create.env
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem-per-cpu=170G # memory pool for all cores
#SBATCH --time=00:50:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#conda create -n nanocompare python=3.6.8

#conda activate nanocompare

#conda info --envs

# Get nextflow from internet
#curl -fsSL get.nextflow.io | bash

conda install -c bioconda nanopolish=0.13.2

conda install -c bioconda deepmod

conda install -c bioconda ont-tombo

conda install -c bioconda megalodon

pip install deepsignal

conda install -c bioconda samtools

conda install -c bioconda minimap2

pip install ont-pyguppy-client-lib==4.2.2

conda install -c bioconda pybedtools

conda env export > environment.yml

docker build -t nanome .

# How to run script/command using docker images and mapping folders
docker run -v `pwd`:/usr/src/nanocompare -w /usr/src/nanocompare -t nanome conda run -n nanocompare python /usr/src/nanocompare/src/plot_figure.py


# Singularity build in JAX HPC
module load singularity
builder="singularity run http://s3-far.jax.org/builder/builder"
$builder nanome.def nanome.sif

singularity shell nanocompare.sif

singularity exec nanocompare.sif /opt/conda/envs/nanocompare/bin/nanopolish --version
