#!/bin/bash
#SBATCH --job-name=basecall
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 250g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

set +x
source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh

set -x

processors=8

numk=$((SLURM_ARRAY_TASK_ID-1))

echo JOBINDEX=${numk}

numKSeptInputDir=${septInputDir}/${numk}
numKBasecallOutputDir=${basecallOutputDir}/${numk}

echo NumKSeptInputDir=${numKSeptInputDir}
echo NumKBasecallOutputDir=${numKBasecallOutputDir}

mkdir -p ${numKBasecallOutputDir}

set +x
conda activate nanoai
set -x

## Run Basecalling with Albacore:
read_fast5_basecaller.py -o fastq,fast5 -t $processors -s ${numKBasecallOutputDir} -i ${numKSeptInputDir} -c r94_450bps_linear.cfg -n 20000000000

set +x
conda deactivate
set -x

echo "###   Basecalling DONE    ###"