#!/bin/bash

###################################################################################
### Run pipeline for multiple nanopore tool in ToolList                         ###
###################################################################################
# Change from base src/ dir to methcall dir
#cd "$(dirname "$0")"/nanocompare/methcall

set +x
source utils.common.sh
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

	# Remove pre-processing dataset dir if needed
	rm -rf ${untaredInputDir} ${septInputDir}
	mkdir -p ${untaredInputDir}
	mkdir -p ${septInputDir}
	mkdir -p ${septInputDir}/log

	prep_ret=$(sbatch --job-name=prep.fast5.${dsname}.N${targetNum} --output=${septInputDir}/log/%x.%j.out --error=${septInputDir}/log/%x.%j.err --export=dsname=${dsname},Tool=${Tool},targetNum=${targetNum},inputDataDir=${inputDataDir},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},multipleInputs=${multipleInputs} nanocompare.preprocessing.fast5.sh)
	prep_taskid=$(echo ${prep_ret} |grep -Eo '[0-9]+$')
	echo ${prep_ret}
	echo "### Submitted preprocessing task for ${analysisPrefix}."
fi

################################################################################
# Step 2: For each Tool
################################################################################
combine_taskids=""
for Tool in ${ToolList[@]}; do
	### Building output folder configuration, such as
	#	/fastscratch/liuya/nanocompare/K562-Runs/
	#       K562-N50-untar
	#       K562-N50-sept
	#	├── K562-Tombo-N50
	#	    ├── K562-Tombo-N50-base-call
	#	    ├── K562-Tombo-N50-meth-call

	analysisPrefix=${dsname}-${Tool}-N${targetNum}
	basecallOutputDir=${outbasedir}/${analysisPrefix}/${analysisPrefix}-basecall
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

	# If previous step need to depend on
	depend_param=""
	if [ -n "${combine_taskids}" ] ; then
		depend_param="afterok${combine_taskids}"
	fi

	mkdir -p ${outbasedir}/log

	# Use "" wrap space string like "DeepSignal Tombo", then can be passed to export vars
	ToolListStr='"'${ToolList[@]}'"'

	sbatch --job-name=clen.${dsname}.N${targetNum} --output=${outbasedir}/log/%x.%j.out --error=${outbasedir}/log/%x.%j.err --dependency=${depend_param} --export=dsname=${dsname},ToolList="${ToolListStr}",targetNum=${targetNum},analysisPrefix=${analysisPrefix},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},methCallsDir=${methCallsDir},outbasedir=${outbasedir},clean_preprocessing=${clean_preprocessing},clean_basecall=${clean_basecall},clean_methcall=${clean_methcall},tar_basecall=${tar_basecall},tar_methcall=${tar_methcall},run_clean=${run_clean} nanocompare.clean.intermediate.sh
#	sbatch --job-name=clen.${dsname}.N${targetNum} --output=${outbasedir}/log/%x.%j.out --error=${outbasedir}/log/%x.%j.err --dependency=${depend_param} --export=ALL nanocompare.clean.preprocessing.sh
	echo "Submitted clean preprocessing dirs task for ${dsname}-N${targetNum}."

fi

echo "### After finished pipeline for all tools, use following to get results of each Nanopore tool"

echo "find ${outbasedir} -type f -name "*.combine.tsv" -exec ls -lh {} \;"
#find ${outbasedir} -name "*.combine.tsv" -exec ls -lh {} \;
