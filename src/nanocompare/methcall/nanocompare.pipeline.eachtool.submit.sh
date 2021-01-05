#!/bin/bash
################################################################################
# Pipeline for each nanopore methylation tool running on a dataset
# Need to populate the parameters into this script
################################################################################

# set -x

source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/utils.common.sh

################################################################################
# Step 2: Basecalling with Albacore
################################################################################
if [ "${run_basecall}" = true ] ; then
	echo "Step2: basecalling"

	# If previous step need to depend on
	depend_param=""
	if [ "${run_preprocessing}" = true ] ; then
		depend_param="--dependency=afterok:${prep_taskid}"
	fi

	rm -rf ${basecallOutputDir}
	mkdir -p ${basecallOutputDir}
	mkdir -p ${basecallOutputDir}/log

	basecall_task_ret=$(sbatch --job-name=albacore.${analysisPrefix} --ntasks=${processors} --array=1-${targetNum} --output=${basecallOutputDir}/log/%x.batch%a.%j.out --error=${basecallOutputDir}/log/%x.batch%a.%j.err ${depend_param} --export=dsname=${dsname},Tool=${Tool},targetNum=${targetNum},analysisPrefix=${analysisPrefix},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},processors=${processors} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.basecall.albacore.sh)

	# Get the job id as :123_1:123_2, etc.
	base_taskids=$(get_arrayjob_ids "${basecall_task_ret}" "${targetNum}")

	# set -x
	echo ${basecall_task_ret}
	echo "### Submitted all basecalling for ${analysisPrefix} by Albacore array-job finished."
fi
################################################################################
################################################################################
################################################################################


################################################################################
# Step 3: Methylation call
################################################################################
if [ "$run_methcall" = true ] ; then
	echo "Step3: methylation calling"

	rm -rf ${methCallsDir}
	mkdir -p ${methCallsDir}
	mkdir -p ${methCallsDir}/log

	# If previous step need to depend on
	depend_param=""
	if [ "$run_basecall" = true ] ; then
		depend_param="afterok${base_taskids}"
	fi

	exp_param="dsname=${dsname},Tool=${Tool},targetNum=${targetNum},analysisPrefix=${analysisPrefix},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},methCallsDir=${methCallsDir},correctedGroup=${correctedGroup},refGenome=${refGenome},chromSizesFile=${chromSizesFile},run_resquiggling=${run_resquiggling},isGPU=${isGPU},deepsignalModel=${deepsignalModel},deepModModel=${deepModModel},processors=${processors}"

	out_param="${methCallsDir}/log/%x.batch%a.%j.out"
	err_param="${methCallsDir}/log/%x.batch%a.%j.err"

	if [ "${Tool}" = "Tombo" ] ; then
		# Tombo methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=tmb.mcal.${analysisPrefix} --ntasks=${processors} --output=${out_param} --error=${err_param} --array=1-${targetNum} --export=${exp_param} --dependency=${depend_param} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.methcall.tombo.sh)
	fi

	if [ "${Tool}" = "DeepSignal" ] ; then
		# DeepSignal methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=dpsig.mcal.${analysisPrefix} --ntasks=${processors} --output=${out_param} --error=${err_param} --array=1-${targetNum} --export=${exp_param} --dependency=${depend_param} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.methcall.deepsignal.sh)
	fi

	if [ "${Tool}" = "DeepMod" ] ; then
		# DeepSignal methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=dpmod.mcal.${analysisPrefix} --ntasks=${processors} --output=${out_param} --error=${err_param} --array=1-${targetNum} --export=${exp_param} --dependency=${depend_param} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.methcall.deepmod.sh)
	fi

	if [ "${Tool}" = "Nanopolish" ] ; then
		# Nanopolish methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=napol.mcal.${analysisPrefix} --ntasks=${processors} --output=${out_param} --error=${err_param} --array=1-${targetNum} --export=${exp_param} --dependency=${depend_param} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.methcall.nanopolish.sh)
	fi

	# Get the job id as :123_1:123_2, etc.
	meth_taskids=$(get_arrayjob_ids "${meth_arrayjob_ret}" "${targetNum}")

	echo ${meth_arrayjob_ret}
	echo "Submitted ${analysisPrefix} methylation calling array-job finished."
fi
################################################################################
################################################################################
################################################################################


################################################################################
# Step 4: Combining results together
################################################################################
if [ "$run_combine" = true ] ; then
	echo "Step4: combing results for ${analysisPrefix}"

	# If previous step need to depend on
	depend_param=""
	if [ "$run_methcall" = true ] ; then
		depend_param="--dependency=afterok${meth_taskids}"
	fi
	combine_ret=$(sbatch --job-name=cmbi.${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err ${depend_param} --export=dsname=${dsname},Tool=${Tool},targetNum=${targetNum},analysisPrefix=${analysisPrefix},outbasedir=${outbasedir},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},methCallsDir=${methCallsDir},clusterDeepModModel=${clusterDeepModModel},refGenome=${refGenome},chromSizesFile=${chromSizesFile},clean_preprocessing=${clean_preprocessing},clean_basecall=${clean_basecall},tar_basecall=${tar_basecall},tar_methcall=${tar_methcall},run_clean=${run_clean} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.combine.results.sh)
	combine_taskid=$(echo ${combine_ret} |grep -Eo '[0-9]+$')

	combine_taskids=${combine_taskids}:${combine_taskid}
	echo ${combine_ret}

	echo "Submitted combine and clean results task for ${analysisPrefix}."
fi

