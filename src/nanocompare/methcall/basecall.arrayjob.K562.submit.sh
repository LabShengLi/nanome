#!/bin/bash

set -x

#Specify input params:
dsname=K562
inputDataDir=/projects/li-lab/yang/results/2020-12-28/K562-Nanopore_GT18-07372.fast5.tar

untaredInputDir=/fastscratch/liuya/nanocompare/${dsname}_untar-1
septInputDir=/fastscratch/liuya/nanocompare/${dsname}_sept-1
basecallOutputDir=/fastscratch/liuya/nanocompare/${dsname}_basecalled-1

targetNum=50

################################################################################
# Step 1: Preprocessing fast5 files, seperate into 0-(N-1) n folders
################################################################################
prep_ret=$(sbatch --job-name=prep.fast5.${dsname} --export=targetNum=${targetNum},inputDataDir=${inputDataDir},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir}, /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/preprocessing.fast5.sh)
prep_taskid=$(echo ${prep_ret} |grep -Eo '[0-9]+$')


################################################################################
# Step 2: Basecalling with Albacore
################################################################################
rm -rf ${basecallOutputDir}
mkdir -p ${basecallOutputDir}

sbatch --job-name=albacore.${dsname} --array=1-${targetNum} --dependency=afterok:${prep_taskid} --export=septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/basecall.sh

echo "Submitted all jobs finished."
