#!/bin/bash
#SBATCH --job-name=apl.prep
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 5-14:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

# How may seperate file folders used to seperate all fast5 files to # of groups
targetNum=50

# Input
inputDataDir=/projects/li-lab/AML-Nanopore/20180517_180508-18-li-001-GXB01186-001/APL-1750_GT18-06409.fast5.tar


dsname=APL
#inputDataFile=testData_resquiggled1test.tar.gz

# Output
untaredInputDir=/fastscratch/liuya/nanocompare/${dsname}_untar

septInputDir=/fastscratch/liuya/nanocompare/${dsname}_sept


scriptSeptFile=/projects/liuya/workspace/tcgajax/nanocompare/methcall/FilesSeparatorNew.py


mkdir -p ${untaredInputDir}
mkdir -p ${septInputDir}


tar -xf ${inputDataDir} -C ${untaredInputDir}


#
#filelist=$(find ${inputDataDir} -name "*.fast5.tar")
#
#for fast5tar in $filelist; do
#    tar -xf $fast5tar -C ${untaredInputDir}
#done
#exit 0


# Seperate fast5 files into $targetNum
python ${scriptSeptFile} ${untaredInputDir} $targetNum ${septInputDir}

