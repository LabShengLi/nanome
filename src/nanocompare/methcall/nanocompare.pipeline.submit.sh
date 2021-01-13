#!/bin/bash

###################################################################################
### Run pipeline for multiple nanopore tool in ToolList                         ###
###################################################################################

# Working dir now must be ${NanoCmpDir}/src/nanocompare/methcall
set +x
source ../../utils.common.sh
set -x

set -e
set -x

mkdir -p ${outbasedir}

# Do preprocessing only once for all tools
untaredInputDir=${outbasedir}/${dsname}-N${targetNum}-untar
septInputDir=${outbasedir}/${dsname}-N${targetNum}-sept

################################################################################
# Step 1: Pre-processing (Untar, seperate files)
################################################################################
if [ "${run_preprocessing}" = true ] ; then
	echo "Step1: pre-processing"

	# Remove previous pre-processing dataset dir if we redo it
	rm -rf ${untaredInputDir} ${septInputDir}
	mkdir -p ${untaredInputDir}
	mkdir -p ${septInputDir}
	mkdir -p ${septInputDir}/log

	prep_ret=$(sbatch --job-name=prep.fast5.${dsname}.N${targetNum} --output=${septInputDir}/log/%x.%j.out --error=${septInputDir}/log/%x.%j.err --export=ALL,dsname=${dsname},Tool=${Tool},targetNum=${targetNum},inputDataDir=${inputDataDir},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},multipleInputs=${multipleInputs} nanocompare.preprocessing.sh)
	prep_taskid=$(echo ${prep_ret} |grep -Eo '[0-9]+$')
	echo ${prep_ret}
	echo "### Submitted preprocessing task for ${analysisPrefix}."
fi



################################################################################
# Step 2: Basecalling
################################################################################
if [ "${run_basecall}" = true ] ; then
	echo "Step2: basecalling"
	basecallOutputDir=${outbasedir}/${dsname}-N${targetNum}-basecall

	# If previous step need to depend on
	depend_param=""
	if [ "${run_preprocessing}" = true ] ; then
		depend_param="--dependency=afterok:${prep_taskid}"
	fi

	rm -rf ${basecallOutputDir}
	mkdir -p ${basecallOutputDir}
	mkdir -p ${basecallOutputDir}/log

	if [ "${basecall_name}" = "Albacore" ] ; then
		# Albacore basecall
		basecall_task_ret=$(sbatch --job-name=bascal.${basecall_name} --ntasks=${processors} --array=1-${targetNum} --output=${basecallOutputDir}/log/%x.batch%a.%j.out --error=${basecallOutputDir}/log/%x.batch%a.%j.err ${depend_param} --export=ALL,dsname=${dsname},Tool=${Tool},targetNum=${targetNum},analysisPrefix=${analysisPrefix},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},processors=${processors},basecall_name=${basecall_name} nanocompare.basecall.sh)

	elif [ "${basecall_name}" = "Guppy" ] ; then
		# Use GPU of Winter server to do Guppy basecall
		basecall_task_ret=$(sbatch --job-name=bascal.${basecall_name}.${analysisPrefix} --ntasks=${processors} --array=1-${targetNum} --output=${basecallOutputDir}/log/%x.batch%a.%j.out --error=${basecallOutputDir}/log/%x.batch%a.%j.err ${depend_param} --export=ALL,dsname=${dsname},Tool=${Tool},targetNum=${targetNum},analysisPrefix=${analysisPrefix},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},processors=${processors},basecall_name=${basecall_name} nanocompare.basecall.sh)
	fi

	# Get the job id as :123_1:123_2, etc.
	base_taskids=$(get_arrayjob_ids "${basecall_task_ret}" "${targetNum}")

	# set -x
	echo ${basecall_task_ret}
	echo "### Submitted all basecalling for ${analysisPrefix} by basecall_name=${basecall_name} array-job finished."
fi
################################################################################
################################################################################
################################################################################


################################################################################
# Step 3: For each Tool, resquiggling and methylation call
################################################################################
combine_taskids=""
for Tool in ${ToolList[@]}; do
	### Building output folder configuration, such as
	#	/fastscratch/liuya/nanocompare/K562-Runs/
	#       K562-N50-untar
	#       K562-N50-sept
	#	├── K562-Tombo-N50
	#	    ├── K562-Tombo-N50-basecall
	#	    ├── K562-Tombo-N50-methcall

	analysisPrefix=${dsname}-${Tool}-N${targetNum}

	methCallsDir=${outbasedir}/${analysisPrefix}/${analysisPrefix}-methcall
	### Start script for submiting jobs
	echo "########################################################################"
	echo "Start pipeline submit for tool with analysisPrefix=${analysisPrefix}"
	echo "Basedir=${outbasedir}/${analysisPrefix}"
	source nanocompare.pipeline.eachtool.submit.sh
	echo ""
done
###################################################################################
###################################################################################
###################################################################################

set -x

################################################################################
# Step 3: Clean after each tool finished
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
