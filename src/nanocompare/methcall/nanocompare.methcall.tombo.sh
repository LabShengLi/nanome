#!/bin/bash
#SBATCH --job-name=tombo.methcall
##SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH --mem=150g
#SBATCH --time=06:00:00
#SBATCH -o log/%x.%j.out
#SBATCH -e log/%x.%j.err
##SBATCH --array=1-11

################################################################################
# Tombo methylation call workflow
# Need to populate the parameters into this script
# Step 1: Tombo resquiglling; Step 2: Tombo detect_modifications; Step 3: Extract Tombo per_read_stats
################################################################################
#cd "$(dirname "$0")"

set +e
set +u
set +x
source ../../utils.common.sh
set -x

job_index=$((SLURM_ARRAY_TASK_ID))

## Modify directory for processed files after basecalling:
resquiggleDir=${outbasedir}/${dsname}-N${targetNum}-resquiggle

mkdir -p ${resquiggleDir}

## Original basecalled results dir, we do not want resquiggling modify it
jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

## Resquiggle results, firstly copy from basecall, then do resquiggling
jobkResquiggleOutputDir=${resquiggleDir}/${job_index}
processedFast5DIR=${jobkResquiggleOutputDir}/workspace

echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "targetNum: ${targetNum}"
echo "analysisPrefix: ${analysisPrefix}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "jobkBasecallOutputDir: ${jobkBasecallOutputDir}"
echo "jobkBasecallOutputDir: ${jobkResquiggleOutputDir}"
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

if [ "${run_resquiggling}" = true ] ; then
	## Re-squiggle the data:
	## TODO: Only can be done once currently? Why --overwrite is not allow multiple runs? Ref https://github.com/nanoporetech/tombo/issues/5

	# See ref of tombo resquiggle by tombo resquiggle -h

	# Firstly clean it, then copy basecall results into resquiggle dir
	rm -rf ${jobkResquiggleOutputDir}
	mkdir -p ${jobkResquiggleOutputDir}

	cp -rf ${jobkBasecallOutputDir}/* ${jobkResquiggleOutputDir}/

	time tombo resquiggle --processes ${processors} --corrected-group ${correctedGroup} \
		--basecall-group Basecall_1D_000 --overwrite ${processedFast5DIR} ${refGenome} #--quiet
	echo "###   Re-squiggling DONE"
fi

## Call methylation from resquiggled fast5 files:
# TODO: what is difference of 5mC and CpG, Ref: https://nanoporetech.github.io/tombo/modified_base_detection.html#
time tombo detect_modifications alternative_model --fast5-basedirs ${processedFast5DIR} \
		--dna --standard-log-likelihood-ratio --statistics-file-basename ${methCallsDir}/${analysisPrefix}.batch_${job_index} \
		--per-read-statistics-basename ${methCallsDir}/${analysisPrefix}.batch_${job_index} \
		--alternate-bases CpG --processes ${processors} --corrected-group ${correctedGroup} #--quiet

echo "###   Tombo methylation calling DONE"


set +x
conda deactivate
set -x

echo "###   Tombo Meth-call task is DONE    ###"
