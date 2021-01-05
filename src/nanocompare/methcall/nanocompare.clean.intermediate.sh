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
# Clean nanopore pipeline intermediate file dirs
################################################################################
set -e
set -x

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "analysisPrefix: ${analysisPrefix}"
echo "outbasedir: ${outbasedir}"
echo "untaredInputDir: ${untaredInputDir}"
echo "septInputDir: ${septInputDir}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "tar_basecall: ${tar_basecall}"
echo "tar_methcall: ${tar_methcall}"
echo "clean_preprocessing: ${clean_preprocessing}"
echo "clean_basecall: ${clean_basecall}"
echo "##################"
set +u

if [ "${tar_basecall}" = true ] ; then
	tar -czf ${outbasedir}/${analysisPrefix}/${analysisPrefix}-basecall.tar.gz ${outbasedir}/${analysisPrefix}/${analysisPrefix}
	echo "### tar basecall dir OK"
fi

if [ "${tar_methcall}" = true ] ; then
	tar -czf ${outbasedir}/${analysisPrefix}/${analysisPrefix}-meth-call.tar.gz ${outbasedir}/${analysisPrefix}/${analysisPrefix}-meth-call
	echo "### tar methcall dir OK"
fi

if [ "${clean_basecall}" = true ] ; then
	rm -rf ${basecallOutputDir}
	echo "### clean basedcalling dir OK"
fi