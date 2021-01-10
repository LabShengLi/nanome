#!/bin/bash
#SBATCH --job-name=prep.fast5
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=200g # memory pool for all cores
#SBATCH --time=1-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

################################################################################
# Pre-processing workflow
# Need to populate the parameters into this script
################################################################################
#cd "$(dirname "$0")"

set -x

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "targetNum: ${targetNum}"
echo "inputDataDir: ${inputDataDir}"
echo "untaredInputDir: ${untaredInputDir}"
echo "septInputDir: ${septInputDir}"
echo "multipleInputs: ${multipleInputs}"
echo "##################"
set +u

# How may seperate file folders used to seperate all fast5 files to # of groups
# ${targetNum}=N, 1-N: will create 0 -- (N-1) folder names

if [ "${multipleInputs}" = true ] ; then
	echo "Multiple fast5.tar need to be untared"

#	filelist=$(find ${inputDataDir} -name "*.fast5.tar")

	filelist=$(find ${inputDataDir} -type f \( -iname \*.fast5.tar -o -iname \*.fast5.tar.gz -o -iname \*.fast5 \))

	echo "Input file list:${filelist}"

	for fast5tarfn in ${filelist}; do
		if [ "${fast5tarfn##*.}" = "tar" ]; then
			tar -xf ${fast5tarfn} -C ${untaredInputDir}
        elif [ "${fast5tarfn##*.}" = "gz" ]; then
			tar -xzf ${fast5tarfn} -C ${untaredInputDir}
        elif [ "${fast5tarfn##*.}" = "fast5" ]; then
			cp ${fast5tarfn} ${untaredInputDir}
        fi
	done
else
	if [ "${inputDataDir##*.}" = "tar" ]; then
		tar -xf ${inputDataDir} -C ${untaredInputDir}
	elif [ "${inputDataDir##*.}" = "gz" ]; then
		tar -xzf ${inputDataDir} -C ${untaredInputDir}
	fi
fi

echo "### Untar input done. ###"

# Seperate fast5 files into N=$targetNum, 1-N folders
time python FilesSeparatorNew.py ${untaredInputDir} ${targetNum} ${septInputDir}

echo "### Seperation fast5 files done. ###"


