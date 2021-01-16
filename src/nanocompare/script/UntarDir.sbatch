#!/bin/bash
#SBATCH --job-name=UntarDir
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 2-14:00 # time (D-HH:MM)
#SBATCH -o %x.%j.figures # STDOUT
#SBATCH -e %x.%j.err # STDERR

#####################################
#
# Untar fast5.tar from dir
#
# Args:
#   $1  -   input dir
#   $2  -   untared dir#
#
#####################################

set -u

inputDataDir=$1
untarDir=$2

tarTag="-xf"

mkdir -p ${untarDir}

filelist=$(find ${inputDataDir} -name "*.fast5.tar")

for fast5tar in $filelist; do
    tar ${tarTag} $fast5tar -C ${untarDir}
done

echo "### untar done. ###"
