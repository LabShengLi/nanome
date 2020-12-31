#!/bin/bash
################################################################################
# Pipeline for each nanopore methylation tool running on a dataset
# Need to populate the parameters into this script
################################################################################

# set -x

source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/common-utils.sh

################################################################################
# Step 1: Pre-processing (Untar, seperate files)
################################################################################
if [ "$run_preprocessing" = true ] ; then
	echo Step1: pre-processing

	# Remove basedir if pre-processing
	rm -rf /fastscratch/liuya/nanocompare/${analysisPrefix}
	mkdir -p /fastscratch/liuya/nanocompare/${analysisPrefix}

	rm -rf ${untaredInputDir}
	rm -rf ${septInputDir}

	mkdir -p ${untaredInputDir}
	mkdir -p ${septInputDir}
	mkdir -p ${septInputDir}/log

	prep_ret=$(sbatch --job-name=prep.fast5.${analysisPrefix} --output=${septInputDir}/log/%x.%j.out --error=${septInputDir}/log/%x.%j.err --export=targetNum=${targetNum},inputDataDir=${inputDataDir},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir}, /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/preprocessing.fast5.sh)
	prep_taskid=$(echo ${prep_ret} |grep -Eo '[0-9]+$')
	echo ${prep_ret}
fi

################################################################################
################################################################################
################################################################################


################################################################################
# Step 2: Basecalling with Albacore
################################################################################
if [ "$run_basecall" = true ] ; then
	echo Step2: basecalling

	depend_param=""
	if [ "$run_preprocessing" = true ] ; then
		depend_param="--dependency=afterok:${prep_taskid}"
	fi

	rm -rf ${basecallOutputDir}
	mkdir -p ${basecallOutputDir}
	mkdir -p ${basecallOutputDir}/log

	basecall_task_ret=$(sbatch --job-name=albacore.${analysisPrefix} --array=1-${targetNum} --output=${basecallOutputDir}/log/%x.%j.out --error=${basecallOutputDir}/log/%x.%j.err ${depend_param} --export=septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/basecall.sh)

	base_taskids=$(get_arrayjob_ids "${basecall_task_ret}" "${targetNum}")

	# set -x
	echo ${basecall_task_ret}
	echo "Submitted all basecalling by Albacore array-job finished."
fi
################################################################################
################################################################################
################################################################################


################################################################################
# Step 3: Methylation call
################################################################################
if [ "$run_methcall" = true ] ; then
	echo Step3: methylation calling

	depend_param=""
	if [ "$run_basecall" = true ] ; then
		depend_param="--dependency=afterok${base_taskids}"
	fi

	rm -rf ${methCallsDir}
	mkdir -p ${methCallsDir}
	mkdir -p ${methCallsDir}/log

	if [ "${Tool}" = "Tombo" ] ; then
		# Tombo methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=tombo-methcall-${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err --array=1-${targetNum} --export=septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},dataname=${dsname},methCallsDir=${methCallsDir},analysisPrefix=${analysisPrefix},correctedGroup=${correctedGroup},refGenome=${refGenome},chromSizesFile=${chromSizesFile},run_resquiggling=${run_resquiggling} ${depend_param} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/methcall.tombo.sh)
	fi

	if [ "${Tool}" = "DeepSignal" ] ; then
		# DeepSignal methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=deepsignal-methcall-${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err --array=1-${targetNum} --export=septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},dataname=${dsname},methCallsDir=${methCallsDir},analysisPrefix=${analysisPrefix},correctedGroup=${correctedGroup},refGenome=${refGenome},chromSizesFile=${chromSizesFile},run_resquiggling=${run_resquiggling},deepsignalModel=${deepsignalModel},isGPU=${isGPU} ${depend_param} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/methcall.deepsignal.sh)
	fi

	if [ "${Tool}" = "DeepMod" ] ; then
		# DeepSignal methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=deepmod-methcall-${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err --array=1-${targetNum} --export=septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},dataname=${dsname},methCallsDir=${methCallsDir},analysisPrefix=${analysisPrefix},refGenome=${refGenome},deepModModel=${deepModModel},isGPU=${isGPU} ${depend_param} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/methcall.deepmod.sh)
	fi

	if [ "${Tool}" = "Nanopolish" ] ; then
		# Nanopolish methylation call pipeline
		meth_arrayjob_ret=$(sbatch --job-name=nanopolish-methcall-${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err --array=1-${targetNum} --export=basecallOutputDir=${basecallOutputDir},dataname=${dsname},methCallsDir=${methCallsDir},analysisPrefix=${analysisPrefix},refGenome=${refGenome},isGPU=${isGPU},run_resquiggling=${run_resquiggling} ${depend_param} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/methcall.nanopolish.sh)
	fi

	meth_taskids=$(get_arrayjob_ids "${meth_arrayjob_ret}" "${targetNum}")

	echo ${meth_arrayjob_ret}
	echo "Submitted all ${Tool} methylation calling array-job finished."
fi
################################################################################
################################################################################
################################################################################


################################################################################
# Step 4: Combining results together
################################################################################
if [ "$run_combine" = true ] ; then
	echo Step4: combing results
	depend_param=""
	if [ "$run_methcall" = true ] ; then
		depend_param="--dependency=afterok${meth_taskids}"
	fi
	sbatch --job-name=combine.${Tool}.${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err ${depend_param} --export=analysisPrefix=${analysisPrefix},methCallsDir=${methCallsDir},Tool=${Tool},dsname=${dsname},clusterDeepModModel=${clusterDeepModModel},refGenome=${refGenome} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/combine.results.sh

	echo "Submitted combine results task for ${Tool}."
fi