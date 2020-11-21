#!/bin/bash

septInputDir=/fastscratch/liuya/nanocompare/HL60_sept
basecallOutputDir=/fastscratch/liuya/nanocompare/HL60_basecalled

targetNum=50

mkdir -p ${basecallOutputDir}

sbatch --array=1-${targetNum} --export=septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir} basecall.sh