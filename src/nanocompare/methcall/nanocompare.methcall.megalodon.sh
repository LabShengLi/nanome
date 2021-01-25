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
# Megalodon methylation call workflow
# Need to populate the parameters into this script
################################################################################

set -e
set +x
source ../../utils.common.sh

set -x

SAMTOOLS_PATH=/home/liuya/anaconda3/envs/nanoai/bin/samtools
GUPPY_MOD_CONFIG="res_dna_r941_min_modbases_5mC_v001.cfg"
GUPPY_TIMEOUT=240
READS_PER_GUPPY_BATCH=100

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

mkdir -p ${methCallsDir}/${job_index}
rm -rf ${methCallsDir}/${job_index}/*
megalodon \
    ${processedFast5DIR} \
    --output-directory ${methCallsDir}/${job_index} \
    --overwrite \
    --outputs basecalls mod_basecalls mappings \
    per_read_mods mods mod_mappings \
    per_read_refs \
    --guppy-server-path ${GuppyDir}/bin/guppy_basecall_server \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    --guppy-params "--num_callers 5 --ipc_threads 80" \
    --reads-per-guppy-batch ${READS_PER_GUPPY_BATCH} \
    --guppy-timeout ${GUPPY_TIMEOUT} \
    --samtools-executable ${SAMTOOLS_PATH} \
    --sort-mappings \
    --mappings-format bam \
    --reference ${refGenome} \
    --mod-motif m CG 0 \
    --mod-output-formats bedmethyl wiggle \
    --write-mods-text \
    --write-mod-log-probs \
    --devices 0 \
    --processes ${processors}

set +x
conda deactivate
set -x

echo "###   Megalodon Meth-call task is DONE    ###"
