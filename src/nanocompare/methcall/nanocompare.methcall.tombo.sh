#!/bin/bash
#SBATCH --job-name=tombo.methcall
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem=150g
#SBATCH --time=1-16:00:00
#SBATCH -o log/%x.%j.out
#SBATCH -e log/%x.%j.err
##SBATCH --array=1-11

##SBATCH -p gpu
##SBATCH --gres=gpu:1
##SBATCH -q inference
################################################################################
# Tombo methylation call workflow
# Need to populate the parameters into this script
################################################################################

set +e
set +u
set +x
source ../../utils.common.sh
set -x

job_index=$((SLURM_ARRAY_TASK_ID))

## Original basecalled results dir, we do not want resquiggling modify it
#jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

## Resquiggle results, firstly copy from basecall, then do resquiggling
#jobkResquiggleOutputDir=${resquiggleDir}/${job_index}
processedFast5DIR=${resquiggleDir}/${job_index}/workspace

echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "targetNum: ${targetNum}"
echo "analysisPrefix: ${analysisPrefix}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "resquiggleDir: ${resquiggleDir}"
echo "processedFast5DIR: ${processedFast5DIR}"
echo "methCallsDir: ${methCallsDir}"
echo "refGenome: ${refGenome}"
echo "run_resquiggling: ${run_resquiggling}"
echo "correctedGroup: ${correctedGroup}"
echo "chromSizesFile: ${chromSizesFile}"
echo "processors: ${processors}"
echo "##################"

set +x
conda activate nanoai
set -x

rm -rf ${methCallsDir}/${analysisPrefix}.batch_${job_index}*

## Call methylation from resquiggled fast5 files:
# TODO: what is difference of 5mC and CpG, Ref: https://nanoporetech.github.io/tombo/modified_base_detection.html#
time tombo detect_modifications alternative_model --fast5-basedirs ${processedFast5DIR} \
		--dna --standard-log-likelihood-ratio --statistics-file-basename \
		${methCallsDir}/${analysisPrefix}.batch_${job_index} \
		--per-read-statistics-basename ${methCallsDir}/${analysisPrefix}.batch_${job_index} \
		--alternate-bases CpG --processes ${processors} --corrected-group ${correctedGroup} #--quiet

echo "###   Tombo methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   Tombo Meth-call task is DONE    ###"
