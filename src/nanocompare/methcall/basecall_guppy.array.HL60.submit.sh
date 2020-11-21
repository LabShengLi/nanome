#!/bin/bash

septInputDir=/fastscratch/liuya/nanocompare/HL60_sept
basecallOutputDir=/fastscratch/liuya/nanocompare/HL60_basecalled_guppy

targetNum=50

mkdir -p ${basecallOutputDir}

sbatch --array=1-${targetNum} --export=septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir} basecall_guppy.sh