#!/bin/bash
#SBATCH --job-name=clean
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 100g # memory pool for all cores
#SBATCH -t 1-23:00:00 # time (D-HH:MM:SS)
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
echo "untaredInputDir: ${untaredInputDir}"
echo "septInputDir: ${septInputDir}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "clean_preprocessing: ${clean_preprocessing}"
echo "clean_basecall: ${clean_basecall}"
echo "##################"
set +u

if [ "${clean_preprocessing}" = true ] ; then
	rm -rf ${untaredInputDir} ${septInputDir}
	echo "### clean preprocessing dir OK"
fi

if [ "${clean_basecall}" = true ] ; then
	rm -rf ${basecallOutputDir}
	echo "### clean basedcalling dir OK"
fi