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

conda install -c bioconda nanopolish=0.13.2

conda install -c bioconda deepmod

conda install -c bioconda ont-tombo

conda install -c bioconda megalodon

pip install deepsignal

conda install -c bioconda samtools

conda install -c bioconda minimap2


conda env export > nanocompare.yml

