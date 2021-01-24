#!/bin/bash
#SBATCH --job-name=deepsignal.methcall
##SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem=250g # memory pool for all cores
#SBATCH --time=06:00:00
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-50%20
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH --mem=120g

################################################################################
# DeepMod methylation call workflow
# Need to populate the parameters into this script
################################################################################
#cd "$(dirname "$0")"

set -e
set +x
source ../../utils.common.sh

set -x

job_index=$((SLURM_ARRAY_TASK_ID))

## Modify directory for raw batch files
processedFast5DIR=${septInputDir}/${job_index}

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "targetNum: ${targetNum}"
echo "analysisPrefix: ${analysisPrefix}"
echo "septInputDir: ${septInputDir}"
echo "processedFast5DIR: ${processedFast5DIR}"
echo "methCallsDir: ${methCallsDir}"
echo "refGenome: ${refGenome}"
echo "isGPU: ${isGPU}"
echo "processors: ${processors}"
echo "##################"
set +u


set +x
conda activate nanoai
set -x



set +x
conda deactivate
set -x

echo "###   Megalodon Meth-call task is DONE    ###"
