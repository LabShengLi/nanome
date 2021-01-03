#!/bin/bash
################################################################################
# Pipeline for each nanopore methylation tool running on a dataset
# Need to populate the parameters into this script
################################################################################

# set -x

source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/utils.common.sh

################################################################################
# Step 1: Pre-processing (Untar, seperate files)
################################################################################
if [ "$run_preprocessing" = true ] ; then
	echo "Step1: pre-processing"

	# Remove dataset running dir if pre-processing
	rm -rf ${outbasedir}/${analysisPrefix}
	mkdir -p ${outbasedir}/${analysisPrefix}

	mkdir -p ${untaredInputDir}
	mkdir -p ${septInputDir}
	mkdir -p ${septInputDir}/log

	prep_ret=$(sbatch --job-name=prep.fast5.${analysisPrefix} --output=${septInputDir}/log/%x.%j.out --error=${septInputDir}/log/%x.%j.err --export=dsname=${dsname},Tool=${Tool},targetNum=${targetNum},inputDataDir=${inputDataDir},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},analysisPrefix=${analysisPrefix},multipleInputs=${multipleInputs} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.preprocessing.fast5.sh)
	prep_taskid=$(echo ${prep_ret} |grep -Eo '[0-9]+$')
	echo ${prep_ret}
	echo "### Submitted preprocessing task for ${analysisPrefix}."
fi

################################################################################
################################################################################
################################################################################


################################################################################
# Step 2: Basecalling with Albacore
################################################################################
if [ "$run_basecall" = true ] ; then
	echo "Step2: basecalling"

	# If previous step need to depend on
	depend_param=""
	if [ "$run_preprocessing" = true ] ; then
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

	# If previous step need to depend on
	depend_param=""
	if [ "$run_basecall" = true ] ; then
		depend_param="afterok${base_taskids}"
	fi

	exp_param="dsname=${dsname},Tool=${Tool},targetNum=${targetNum},analysisPrefix=${analysisPrefix},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},methCallsDir=${methCallsDir},correctedGroup=${correctedGroup},refGenome=${refGenome},chromSizesFile=${chromSizesFile},run_resquiggling=${run_resquiggling},isGPU=${isGPU},deepsignalModel=${deepsignalModel},deepModModel=${deepModModel},processors=${processors}"

	out_param="${methCallsDir}/log/%x.batch%a.%j.out"
	err_param="${methCallsDir}/log/%x.batch%a.%j.err"

	rm -rf ${methCallsDir}
	mkdir -p ${methCallsDir}
	mkdir -p ${methCallsDir}/log

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
	echo "Submitted all ${analysisPrefix} methylation calling array-job finished."
fi
################################################################################
################################################################################
################################################################################


################################################################################
# Step 4: Combining results together
################################################################################
if [ "$run_combine" = true ] ; then
	echo "Step4: combing results"

	# If previous step need to depend on
	depend_param=""
	if [ "$run_methcall" = true ] ; then
		depend_param="--dependency=afterok${meth_taskids}"
	fi
	combine_ret=$(sbatch --job-name=cmbi.${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err ${depend_param} --export=analysisPrefix=${analysisPrefix},methCallsDir=${methCallsDir},Tool=${Tool},dsname=${dsname},clusterDeepModModel=${clusterDeepModModel},refGenome=${refGenome} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.combine.results.sh)
	combine_taskid=$(echo ${combine_ret} |grep -Eo '[0-9]+$')
	echo ${combine_ret}

	echo "Submitted combine results task for ${analysisPrefix}."
fi

################################################################################
# Step 5: Clean intermediate dirs
################################################################################
if [ "${run_clean}" = true ] ; then
	echo "Step5: clean dirs"

	# If previous step need to depend on
	depend_param=""
	if [ "${run_combine}" = true ] ; then
		depend_param="afterok:${combine_taskid}"
	fi
	sbatch --job-name=clen.${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err --dependency=${depend_param} --export=dsname=${dsname},Tool=${Tool},targetNum=${targetNum},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},analysisPrefix=${analysisPrefix},basecallOutputDir=${basecallOutputDir},clean_preprocessing=${clean_preprocessing},clean_basecall=${clean_basecall} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.clean.intermediate.sh
	echo "Submitted clean dirs task for ${analysisPrefix}."
fi
