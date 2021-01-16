#!/bin/bash
#SBATCH --job-name=Sept
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 2-14:00 # time (D-HH:MM)
#SBATCH -o %x.%j.figures # STDOUT
#SBATCH -e %x.%j.err # STDERR

#####################################
#
# Seperate fast5 files into # targetNum of folders (0-N-1). Note: job2 waiting for job1 finished and returned ok.
#
# Args:
#   $1  -   untared fast5 dir
#   $2  -   seperated dir
#   $3  -   # of seperation folders (0 to N-1)
#
#####################################

set -u

untarDir=$1
septDir=$2
targetNum=$3

FileSeptScript=/projects/liuya/workspace/tcgajax/nanocompare/methcall/FilesSeparatorNew.py

# Seperate fast5 files into $targetNum
mkdir -p ${septDir}

python ${FileSeptScript} ${untarDir} $targetNum ${septDir}

echo "### seperate done. ###"
