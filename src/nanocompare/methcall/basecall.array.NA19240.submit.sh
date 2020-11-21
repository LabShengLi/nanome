#!/bin/bash

septInputDir=/fastscratch/liuya/nanocompare/NA19240_sept
basecallOutputDir=/fastscratch/liuya/nanocompare/NA19240_basecalled

targetNum=100

mkdir -p ${basecallOutputDir}

sbatch --array=1-${targetNum} --export=septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir} basecall.sh