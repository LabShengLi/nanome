#!/bin/bash
################################################################################
# Pipeline for each nanopore methylation tool running on a dataset
# Need to populate the parameters into this script
################################################################################
# This script will be sourced by nanocompare.pipeline.submit.sh, so working dir is same with it

#set -x
#set -e
#set -u

################################################################################
# Step 4.1: Methylation call
################################################################################
if [ "${run_methcall}" = true ] ; then
	echo "Step4.1: methylation calling"

#	rm -rf ${methCallsDir}
	mkdir -p ${methCallsDir}/log

	# If previous step need to depend on
	depend_param=""
	if [ "$run_basecall" = true ] ; then
		depend_param="afterok${base_taskids}"
	fi

	exp_param="dsname=${dsname},Tool=${Tool},targetNum=${targetNum},outbasedir=${outbasedir},analysisPrefix=${analysisPrefix},septInputDir=${septInputDir},basecall_name=${basecall_name},basecallOutputDir=${basecallOutputDir},methCallsDir=${methCallsDir},correctedGroup=${correctedGroup},refGenome=${refGenome},chromSizesFile=${chromSizesFile},run_resquiggling=${run_resquiggling},isGPU=${isGPU},deepsignalModel=${deepsignalModel},deepModModel=${deepModModel},processors=${processors},resquiggleDir=${resquiggleDir}"

	out_param="${methCallsDir}/log/%x.batch%a.%j.out"
	err_param="${methCallsDir}/log/%x.batch%a.%j.err"

	if [ "${Tool}" = "Tombo" ] ; then
		# Tombo methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=tmb.mcal.${analysisPrefix} --ntasks=${processors} --output=${out_param} --error=${err_param} --array=1-${targetNum} --export=ALL,${exp_param} --dependency=${depend_param} nanocompare.methcall.tombo.sh)
	elif [ "${Tool}" = "DeepSignal" ] ; then
		# DeepSignal methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=dpsig.mcal.${analysisPrefix} --ntasks=${processors} --output=${out_param} --error=${err_param} --array=1-${targetNum} --export=ALL,${exp_param} --dependency=${depend_param} nanocompare.methcall.deepsignal.sh)
	elif [ "${Tool}" = "DeepMod" ] ; then # TODO: we use only 2 cores for speed up
		# DeepMod methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=dpmod.mcal.${analysisPrefix} --ntasks=6 --output=${out_param} --error=${err_param} --array=1-${targetNum}%50 --export=ALL,${exp_param} --dependency=${depend_param} nanocompare.methcall.deepmod.sh)

#		meth_arrayjob_ret=$(sbatch --job-name=dpmod.mcal.${analysisPrefix} --ntasks=6 --output=${out_param} --error=${err_param} --array=147,148,1,166,267,294,297,298 --export=ALL,${exp_param} --dependency=${depend_param} nanocompare.methcall.deepmod.sh)

	elif [ "${Tool}" = "Nanopolish" ] ; then
		# Nanopolish methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=napol.mcal.${analysisPrefix} --ntasks=${processors} --output=${out_param} --error=${err_param} --array=1-${targetNum} --export=ALL,${exp_param} --dependency=${depend_param} nanocompare.methcall.nanopolish.sh)
	elif [ "${Tool}" = "Guppy" ] ; then
		# Must run on GPU now
		meth_arrayjob_ret=$(sbatch --job-name=guppy.mcal.${analysisPrefix} --ntasks=${processors} --output=${out_param} --error=${err_param} --array=1-${targetNum} --export=ALL,${exp_param} --dependency=${depend_param} nanocompare.methcall.guppy.sh)
	elif [ "${Tool}" = "Megalodon" ] ; then
		echo "will add Megalodon"
	else
		echo "Tool=${Tool} is not support now"
		exit -1
	fi

	# Get the job id as :123_1:123_2, etc.
#	meth_taskids=$(get_arrayjob_ids "${meth_arrayjob_ret}" "${targetNum}")

	meth_arrayjob_id=${meth_arrayjob_ret##* }
	meth_taskids="afterok:${meth_arrayjob_id}_[1-${targetNum}]"

	echo ${meth_arrayjob_ret}
	echo "Submitted ${analysisPrefix} methylation calling array-job finished."
fi
################################################################################
################################################################################
################################################################################


################################################################################
# Step 4.2: Combining results together after all batches methylation call finished
################################################################################
if [ "${run_combine}" = true ] ; then
	echo "Step4.2: combing results for ${analysisPrefix}"

	# If previous step need to depend on
	depend_param=""
	if [ "$run_methcall" = true ] ; then
		depend_param="--dependency=${meth_taskids}"
	fi
	combine_ret=$(sbatch --job-name=cmbi.${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err ${depend_param} --export=ALL,dsname=${dsname},Tool=${Tool},targetNum=${targetNum},analysisPrefix=${analysisPrefix},outbasedir=${outbasedir},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},methCallsDir=${methCallsDir},clusterDeepModModel=${clusterDeepModModel},refGenome=${refGenome},chromSizesFile=${chromSizesFile},clean_preprocessing=${clean_preprocessing},clean_basecall=${clean_basecall},tar_basecall=${tar_basecall},tar_methcall=${tar_methcall},run_clean=${run_clean} nanocompare.combine.results.sh)
	combine_taskid=$(echo ${combine_ret} |grep -Eo '[0-9]+$')

	combine_taskids=${combine_taskids}:${combine_taskid}
	echo ${combine_ret}

	echo "Submitted combine and clean results task for ${analysisPrefix}."
fi

