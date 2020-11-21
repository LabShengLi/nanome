#!/bin/bash
#SBATCH --job-name=BaseCall.Albacore
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 2-20:00 # time (D-HH:MM)
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR
##SBATCH --array=1-11
bash

set -u

processors=8
numk=$((SLURM_ARRAY_TASK_ID-1))

septDir=$1
basecallDir=$2

numKSeptDir=${septDir}/${numk}
numKBasecallDir=${basecallDir}/${numk}

echo numKSeptDir=${numKSeptDir}
echo numKBasecallDir=${numKBasecallDir}

mkdir -p ${numKBasecallDir}

#source /projects/liuya/workspace/tcgajax/nanocompare/methcall/conda_setup.sh

conda activate nanoai

## Run Basecalling with Albacore:
read_fast5_basecaller.py -o fastq,fast5 -t $processors -s ${numKBasecallDir} -i ${numKSeptDir} -c r94_450bps_linear.cfg -n 20000000000

conda deactivate

echo "###   Basecalling DONE    ###"
