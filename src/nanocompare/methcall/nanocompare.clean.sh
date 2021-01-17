#!/bin/bash
#SBATCH --job-name=clean
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=100g # memory pool for all cores
#SBATCH --time=01:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

################################################################################
# Clean nanopore pipeline pre-processing file dirs
################################################################################
set -e
set -x

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "ToolList: ${ToolList}"
echo "outbasedir: ${outbasedir}"
echo "untaredInputDir: ${untaredInputDir}"
echo "septInputDir: ${septInputDir}"
echo "septInputDir: ${basecallOutputDir}"
echo "septInputDir: ${methCallsDir}"
echo "clean_preprocessing: ${clean_preprocessing}"
echo "tar_basecall: ${tar_basecall}"
echo "tar_methcall: ${tar_methcall}"
echo "clean_basecall: ${clean_basecall}"
echo "clean_methcall: ${clean_methcall}"
echo "##################"

set +u
set +e

ToolList="${ToolList%\"}"
ToolList="${ToolList#\"}"

ToolList=( ${ToolList} )
echo "${ToolList[@]}"

for Tool in ${ToolList[@]}; do
	echo "Tool=${Tool}"
	analysisPrefix=${dsname}-${Tool}-N${targetNum}
	if [ "${tar_basecall}" = true ] ; then
		tar -czf ${outbasedir}/${analysisPrefix}/${analysisPrefix}-basecall.tar.gz ${outbasedir}/${analysisPrefix}/${analysisPrefix}-basecall
		echo "### tar basecall dir OK"
	fi

	if [ "${tar_methcall}" = true ] ; then
		tar -czf ${outbasedir}/${analysisPrefix}/${analysisPrefix}-meth-call.tar.gz ${outbasedir}/${analysisPrefix}/${analysisPrefix}-methcall
		echo "### tar methcall dir OK"
	fi

	if [ "${clean_basecall}" = true ] ; then
		rm -rf ${outbasedir}/${analysisPrefix}/${analysisPrefix}-basecall
		echo "### clean basecalling dir OK"
	fi

	if [ "${clean_methcall}" = true ] ; then
		rm -rf ${outbasedir}/${analysisPrefix}/${analysisPrefix}-methcall
		echo "### clean methcalling dir OK"
	fi

done


if [ "${clean_preprocessing}" = true ] ; then
	rm -rf ${untaredInputDir} ${septInputDir}
	echo "### clean preprocessing dir OK"
fi