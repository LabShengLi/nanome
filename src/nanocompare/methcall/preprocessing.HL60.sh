#!/bin/bash
#SBATCH --job-name=hl60.untar
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
inputDataDir=/projects/li-lab/NanoporeData/Leukemia_ONT/20180612_180601-18-li-004-GXB01102-002/HL60-Nanopore_GT18-07373.fast5.tar
#inputDataFile=testData_resquiggled1test.tar.gz

# Output
untaredInputDir=/fastscratch/liuya/nanocompare/HL60_untar

septInputDir=/fastscratch/liuya/nanocompare/HL60_sept

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
#
#exit 0

#tar xf ${inputDataDir}/${inputDataFile} -C ${untarInputDir}
#
## tar xzf ${inputDataDir} -C ${preprocessedInputDir}


# Seperate fast5 files into $targetNum
python ${scriptSeptFile} ${untaredInputDir} $targetNum ${septInputDir}


