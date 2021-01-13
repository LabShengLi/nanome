#!/bin/bash
#SBATCH --job-name=deepsignal.methcall
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=250g # memory pool for all cores
#SBATCH --time=15:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# DeepSignal methylation call workflow
# Need to populate the parameters into this script
# Step 1: Tombo resquiglling; Step 2: Deepsignal call_mods
################################################################################
#cd "$(dirname "$0")"

set -e
set +x
#source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh
#source /home/liuya/.bash_profile
source ../../utils.common.sh

set -x

#processors=8

job_index=$((SLURM_ARRAY_TASK_ID))

jobkBasecallOutputDir=${basecallOutputDir}/${job_index}
echo jobkBasecallOutputDir=${jobkBasecallOutputDir}

## Modify directory for processed files after basecalling:
processedFast5DIR=${jobkBasecallOutputDir}/workspace/pass/0

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

if [ "${run_resquiggling}" = true ] ; then
	## Re-squiggle the data:
	time tombo resquiggle $processedFast5DIR $refGenome --processes $processors --corrected-group $correctedGroup \
			--basecall-group Basecall_1D_000 --overwrite #--quiet

	echo "###   Re-squiggling DONE"
fi

## Call methylation from processed fast5 files:
#date
time deepsignal call_mods --input_path $processedFast5DIR --model_path $deepsignalModel \
		--result_file $methCallsDir/$analysisPrefix.batch_${job_index}.CpG.deepsignal.call_mods.tsv \
		--reference_path $refGenome --corrected_group $correctedGroup --nproc $processors --is_gpu $isGPU

echo "###   DeepSignal methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   DeepSignal Meth-call task is DONE    ###"
