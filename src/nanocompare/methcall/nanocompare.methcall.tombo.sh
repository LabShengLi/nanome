#!/bin/bash
#SBATCH --job-name=tombo.methcall
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=250g # memory pool for all cores
#SBATCH --time=2-19:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# Tombo methylation call workflow
# Need to populate the parameters into this script
# Step 1: Tombo resquiglling; Step 2: Tombo detect_modifications; Step 3: Extract Tombo per_read_stats
################################################################################
set -e
set +x

#source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh
source /home/liuya/.bash_profile
set -x

#processors=8

job_index=$((SLURM_ARRAY_TASK_ID-1))
jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

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
echo "processors: ${processors}"
echo "##################"
set +u

set +x
conda activate nanoai
set -x

if [ "${run_resquiggling}" = true ] ; then
	## Re-squiggle the data:
	## TODO: Only can be done once currently? Why --overwrite is not allow multiple runs? Ref https://github.com/nanoporetech/tombo/issues/5

	# See ref of tombo resquiggle by tombo resquiggle -h
	#date
	time tombo resquiggle $processedFast5DIR ${refGenome} --processes ${processors} \
			--corrected-group ${correctedGroup} --basecall-group Basecall_1D_000 --overwrite #--quiet
	#date
	echo "###   Re-squiggling DONE"
fi

## Call methylation from processed fast5 files:


# TODO: what is difference of 5mC and CpG, Ref: https://nanoporetech.github.io/tombo/modified_base_detection.html#
date
time tombo detect_modifications alternative_model --fast5-basedirs $processedFast5DIR \
		--statistics-file-basename $methCallsDir/$analysisPrefix.batch_${job_index} \
		--per-read-statistics-basename $methCallsDir/$analysisPrefix.batch_${job_index} \
		--alternate-bases CpG --processes $processors --corrected-group $correctedGroup \
		--multiprocess-region-size 10000 #--quiet
date

echo "###   Tombo methylation calling DONE"


set +x
conda deactivate
set -x

echo "###   Tombo Meth-call task is DONE    ###"
