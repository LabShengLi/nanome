#!/bin/bash

dsname=APL
septInputDir=/fastscratch/liuya/nanocompare/${dsname}_sept
basecallOutputDir=/fastscratch/liuya/nanocompare/${dsname}_basecalled

targetNum=50

mkdir -p ${basecallOutputDir}

sbatch --job-name=bascal.${dsname} --array=1-${targetNum} --export=septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir} basecall.sh

