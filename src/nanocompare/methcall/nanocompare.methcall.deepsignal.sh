#!/bin/bash
#SBATCH --job-name=deepsignal.methcall
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=150g # memory pool for all cores
#SBATCH --time=1-16:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# DeepSignal methylation call workflow
# Need to populate the parameters into this script
# Step 1: Tombo resquiglling; Step 2: Deepsignal call_mods
################################################################################

set -e
set +x
source ../../utils.common.sh

set -x

job_index=$((SLURM_ARRAY_TASK_ID))

jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

## Modify directory for processed files after basecalling:
#processedFast5DIR=${jobkBasecallOutputDir}/workspace/pass/0
## Resquiggle results, firstly copy from basecall, then do resquiggling
jobkResquiggleOutputDir=${resquiggleDir}/${job_index}
processedFast5DIR=${jobkResquiggleOutputDir}/workspace

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "targetNum: ${targetNum}"
echo "analysisPrefix: ${analysisPrefix}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "jobkBasecallOutputDir: ${jobkBasecallOutputDir}"
echo "processedFast5DIR: ${processedFast5DIR}"
echo "methCallsDir: ${methCallsDir}"
echo "refGenome: ${refGenome}"
echo "run_resquiggling: ${run_resquiggling}"
echo "correctedGroup: ${correctedGroup}"
echo "chromSizesFile: ${chromSizesFile}"
echo "deepsignalModel: ${deepsignalModel}"
echo "isGPU: ${isGPU}"
echo "processors: ${processors}"
echo "##################"
set +u

set +x
conda activate nanoai
set -x

## Call methylation from processed fast5 files:
#date
time deepsignal call_mods --input_path ${processedFast5DIR} --model_path ${deepsignalModel} \
		--result_file ${methCallsDir}/${analysisPrefix}.batch_${job_index}.CpG.deepsignal.call_mods.tsv \
		--reference_path ${refGenome} --corrected_group ${correctedGroup} --nproc ${processors} --is_gpu ${isGPU}

echo "###   DeepSignal methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   DeepSignal Meth-call task is DONE    ###"
