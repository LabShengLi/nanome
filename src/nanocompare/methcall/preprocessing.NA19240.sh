#!/bin/bash
#SBATCH --job-name=na19240
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 5-14:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

# How may seperate file folders used to seperate all fast5 files to # of groups

# Input
inputDataDir=/fastscratch/liuya/nanocompare/NA19240
#inputDataFile=testData_resquiggled1test.tar.gz

# Output
untarDir=/fastscratch/liuya/nanocompare/NA19240_untar2
septDir=/fastscratch/liuya/nanocompare/NA19240_sept2
basecalledDir=/fastscratch/liuya/nanocompare/NA19240_basecalled_guppy2

targetNum=5

dsname=NA19240

#/projects/liuya/workspace/tcgajax/nanocompare/script/PreprocessingTarDir.sh ${inputDataDir} ${untarDir} ${septDir} ${targetNum} ${dsname}


#/projects/liuya/workspace/tcgajax/nanocompare/script/Basecall.Albacore.Array.Submit.sh ${septDir} ${targetNum} ${basecalledDir}

/projects/liuya/workspace/tcgajax/nanocompare/script/Basecall.Guppy.Array.Submit.sh ${septDir} ${targetNum} ${basecalledDir}
