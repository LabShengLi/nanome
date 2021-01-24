#!/bin/bash

###################################################################################
### Run pipeline for multiple nanopore tool in ToolList                         ###
###################################################################################

# Working dir now must be ${NanoCmpDir}/src/nanocompare/methcall
set +x
#source ../../utils.common.sh
set -x

## Currently we use command arg to do four kind of actions
cmd=${1:-noaction} # Preprocess, Basecall, Resquiggle, Methcall

run_preprocessing=false
run_basecall=false
run_resquiggling=false
run_methcall=false
run_combine=false
run_clean=false

if [ "Preprocess" = "${cmd}" ] ; then ## bash APL.nanocompare.pipeline.submit.sh Preprocess
	run_preprocessing=true
elif [ "Basecall" = "${cmd}" ] ; then
	run_basecall=true
elif [ "Resquiggle" = "${cmd}" ] ; then ## bash APL.nanocompare.pipeline.submit.sh Resquiggle
	run_resquiggling=true
elif [ "Methcall" = "${cmd}" ] ; then
	run_methcall=true
#	run_combine=true
elif [ "Combine" = "${cmd}" ] ; then
	run_combine=true
else
	echo "No action needed."
	exit -1
fi

mkdir -p ${outbasedir}

# Do preprocessing only once for all tools
untaredInputDir=${outbasedir}/${dsname}-N${targetNum}-untar
septInputDir=${outbasedir}/${dsname}-N${targetNum}-sept

basecallOutputDir=${outbasedir}/${dsname}-N${targetNum}-basecall

if [ "${basecall_name}" = "Albacore" ] ; then
	basecallOutputDir=${outbasedir}/${dsname}-N${targetNum}-basecall-Albacore
fi

## Modify directory for processed files after basecalling:
resquiggleDir=${outbasedir}/${dsname}-N${targetNum}-resquiggle

################################################################################
# Step 1: Pre-processing (Untar, seperate files)
################################################################################
if [ "${run_preprocessing}" = true ] ; then
	echo "Step1: pre-processing"

	rm -rf ${septInputDir}/log/*
	mkdir -p ${septInputDir}/log

	prep_ret=$(sbatch --job-name=prep.fast5.${dsname}.N${targetNum} --output=${septInputDir}/log/%x.%j.out --error=${septInputDir}/log/%x.%j.err --export=ALL,dsname=${dsname},Tool=${Tool},targetNum=${targetNum},inputDataDir=${inputDataDir},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},multipleInputs=${multipleInputs} nanocompare.preprocessing.sh)
	prep_taskid=$(echo ${prep_ret} |grep -Eo '[0-9]+$')
	echo ${prep_ret}
	echo "### Submitted preprocessing task for ${dsname}.N${targetNum}."
fi

################################################################################
# Step 2: Basecalling
################################################################################
if [ "${run_basecall}" = true ] ; then
	echo "Step2: basecalling"
	# If previous step need to depend on
	depend_param=""
	if [ "${run_preprocessing}" = true ] ; then
		depend_param="--dependency=afterok:${prep_taskid}"
	fi

	mkdir -p ${basecallOutputDir}/log

	if [ "${basecall_name}" = "Albacore" ] ; then
		cpuCFG="--partition=compute --qos=batch --ntasks=8 --time=1-23:00:00 --array=1-${targetNum}%50 --mem=200g"
	elif [ "${basecall_name}" = "Guppy" ] ; then
		cpuCFG="--gres=gpu:1 --partition=gpu --qos=inference --ntasks=4 --time=06:00:00 --array=1-${targetNum}%20 --mem=150g"
	fi

	basecall_task_ret=$(sbatch --job-name=bascal.${basecall_name}.${dsname}.N${targetNum} --output=${basecallOutputDir}/log/%x.batch%a.%j.out --error=${basecallOutputDir}/log/%x.batch%a.%j.err ${cpuCFG} ${depend_param} --export=ALL,dsname=${dsname},Tool=${Tool},targetNum=${targetNum},analysisPrefix=${analysisPrefix},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},processors=${processors},basecall_name=${basecall_name} nanocompare.basecall.sh)

	basecall_arrayjob_id=${basecall_task_ret##* }
	base_taskids="afterok:${basecall_arrayjob_id}_[1-${targetNum}]"

	echo ${basecall_task_ret}
	echo "### Submitted all basecalling for by basecall_name=${basecall_name} array-job finished."
fi
################################################################################
################################################################################
################################################################################

################################################################################
# Step 3: Resquiggling
################################################################################
if [ "${run_resquiggling}" = true ] ; then
	echo "Step3: resquiggling"

	rm -rf ${resquiggleDir}
	mkdir -p ${resquiggleDir}/log

	# If previous step need to depend on
	depend_param=""
	if [ "${run_basecall}" = true ] ; then
		depend_param="--dependency=afterok:${base_taskids}"
	fi

	resquiggling_task_ret=$(sbatch --job-name=rsquigl.${dsname}.N${targetNum} --ntasks=${processors} --array=1-${targetNum} --output=${resquiggleDir}/log/%x.batch%a.%j.out --error=${resquiggleDir}/log/%x.batch%a.%j.err ${depend_param} --export=ALL,dsname=${dsname},targetNum=${targetNum},basecallOutputDir=${basecallOutputDir},processors=${processors},basecall_name=${basecall_name},resquiggleDir=${resquiggleDir} nanocompare.resquiggle.sh)

	# Get the job id as :123_1:123_2, etc.
#	resquiggling_taskids=$(get_arrayjob_ids "${resquiggling_task_ret}" "${targetNum}")
	resquiggling_arrayjob_id=${resquiggling_task_ret##* }
	resquiggling_taskids="afterok:${resquiggling_arrayjob_id}_[1-${targetNum}]"

	echo ${resquiggling_task_ret}
	echo "### Submitted all resquiggling for ${dsname}.N${targetNum} array-job finished."

fi

################################################################################
################################################################################
################################################################################

################################################################################
# Step 4: For each Tool, do methylation call
################################################################################
if [ "${run_methcall}" = true ] || [ "${run_combine}" = true ] ; then
	combine_taskids=""
	for Tool in ${ToolList[@]}; do
		### Building output folder configuration, such as
#		tree /fastscratch/liuya/nanocompare/HL60-Runs/ -L 1
#		/fastscratch/liuya/nanocompare/HL60-Runs/
#		├── HL60-DeepMod-N50
#		├── HL60-DeepSignal-N50
#		├── HL60-N50-basecall
#		├── HL60-N50-resquiggle
#		├── HL60-N50-sept
#		├── HL60-N50-untar
#		├── HL60-Nanopolish-N50
#		└── HL60-Tombo-N50
		analysisPrefix=${dsname}-${Tool}-N${targetNum}
		methCallsDir=${outbasedir}/${analysisPrefix}/${analysisPrefix}-methcall

		if [ "${basecall_name}" = "Albacore" ] ; then
			methCallsDir=${outbasedir}/${analysisPrefix}/${analysisPrefix}-methcall-Albacore
		fi

		### Start script for submiting jobs
		echo "########################################################################"
		echo "Start pipeline submit for tool with analysisPrefix=${analysisPrefix}"
		echo "Basedir=${outbasedir}/${analysisPrefix}"
		source nanocompare.pipeline.eachtool.submit.sh
		echo ""
	done
fi
###################################################################################
###################################################################################
###################################################################################

#set -x

################################################################################
# Step 5: Clean after each tool finished
################################################################################

# After all tools combine task finished, we can clean preprocessing files safely
if [ "${run_clean}" = true ] ; then
	echo "Step5: clean pre-processing results for ${dsname}-N${targetNum}"

	# If previous step need to be depended on
	depend_param=""
	if [ -n "${combine_taskids}" ] ; then
		depend_param="afterok${combine_taskids}"
	fi

	mkdir -p ${outbasedir}/log

	# Use "" wrap space string like '"DeepSignal Tombo"', then can be passed to export vars
	ToolListStr='"'${ToolList[@]}'"'

	sbatch --job-name=clen.${dsname}.N${targetNum} --output=${outbasedir}/log/%x.%j.out --error=${outbasedir}/log/%x.%j.err --dependency=${depend_param} --export=ALL,dsname=${dsname},ToolList="${ToolListStr}",targetNum=${targetNum},analysisPrefix=${analysisPrefix},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},methCallsDir=${methCallsDir},outbasedir=${outbasedir},clean_preprocessing=${clean_preprocessing},clean_basecall=${clean_basecall},clean_methcall=${clean_methcall},tar_basecall=${tar_basecall},tar_methcall=${tar_methcall},run_clean=${run_clean} nanocompare.clean.sh
	echo "Submitted clean preprocessing dirs task for ${dsname}-N${targetNum}."

fi

echo "### After finished pipeline for all tools, use following to get results of each Nanopore tool"

echo "find ${outbasedir} -type f -name "*.combine.tsv" -exec ls -lh {} \;"
#find ${outbasedir} -name "*.combine.tsv" -exec ls -lh {} \;
