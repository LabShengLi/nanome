#!/bin/bash
#SBATCH --job-name=build.singularity
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem-per-cpu=130g
#SBATCH --time=03:59:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

##SBATCH --mem=300g # memory pool for all cores
##SBATCH -p gpu
##SBATCH --gres=gpu:1             # number of gpus per node
##SBATCH -q inference

module load singularity
builder="singularity run http://s3-far.jax.org/builder/builder"

rm -f nanome.sif
$builder nanome.def nanome.sif
