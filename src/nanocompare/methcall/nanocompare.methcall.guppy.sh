#!/bin/bash
#SBATCH --job-name=guppy.methcall
##SBATCH --partition=compute
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem-per-cpu=100G # memory pool for all cores
#SBATCH --time=06:00:00 # time
#SBATCH -o log/%x.%j.figures # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# Guppy methylation call workflow
# Need to populate the parameters into this script
################################################################################

set -e
set +x
source ../../utils.common.sh

set -x

# Only need raw fast5 files as input
job_index=$((SLURM_ARRAY_TASK_ID))
jobkSeptInputDir=${septInputDir}/${job_index}

#jobkBasecallOutputDir=${basecallOutputDir}/${job_index}
jobkMethCallDir=${methCallsDir}/${job_index}

rm -rf ${jobkMethCallDir}
mkdir -p ${jobkMethCallDir}

## Modify directory for processed files after basecalling:
#processedFast5DIR=${jobkBasecallOutputDir}/workspace

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "targetNum: ${targetNum}"
echo "analysisPrefix: ${analysisPrefix}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "jobkSeptInputDir: ${jobkSeptInputDir}"
echo "jobkMethCallDir: ${jobkMethCallDir}"
echo "methCallsDir: ${methCallsDir}"
echo "refGenome: ${refGenome}"
echo "isGPU: ${isGPU}"
echo "processors: ${processors}"
echo "##################"
set +u

set +x
conda activate nanoai
set -x


time ${GuppyDir}/bin/guppy_basecaller --input_path ${jobkSeptInputDir} \
    --save_path ${jobkMethCallDir} --config dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg \
    --gpu_runners_per_device ${processors} --num_callers 3 --fast5_out --verbose_logs --device auto

echo "###   Guppy methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   Guppy Meth-call task is DONE    ###"
