#!/bin/bash
#SBATCH --job-name=prep.fast5
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

# How may seperate file folders used to seperate all fast5 files to # of groups
# 1-N: will create 0 -- (N-1) folder names
#targetNum=50

# Input
#inputDataDir=/projects/li-lab/yang/results/2020-12-28/K562-Nanopore_GT18-07372.fast5.tar

#dsname=K562

# Output
#untaredInputDir=/fastscratch/liuya/nanocompare/${dsname}_untar
#septInputDir=/fastscratch/liuya/nanocompare/${dsname}_sept
scriptSeptFile=/projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/FilesSeparatorNew.py

tar -xf ${inputDataDir} -C ${untaredInputDir}

# Seperate fast5 files into $targetNum
python ${scriptSeptFile} ${untaredInputDir} $targetNum ${septInputDir}


