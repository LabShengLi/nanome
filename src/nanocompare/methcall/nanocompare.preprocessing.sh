#!/bin/bash
#SBATCH --job-name=prep.fast5
#SBATCH --partition=compute
##SBATCH -p gpu
##SBATCH -q inference
##SBATCH --gres=gpu:1
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=250g # memory pool for all cores
##SBATCH --time=06:00:00 # time (D-HH:MM)
#SBATCH --time=1-06:00:00 # time (D-HH:MM) # for NA19240 large input
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
echo "targetNum: ${targetNum}"
echo "inputDataDir: ${inputDataDir}"
echo "untaredInputDir: ${untaredInputDir}"
echo "septInputDir: ${septInputDir}"
echo "multipleInputs: ${multipleInputs}"
echo "##################"
set +u


# Remove previous pre-processing dataset dir if we redo it
rm -rf ${untaredInputDir}/*

mkdir -p ${untaredInputDir}

rm -rf "${septInputDir}/!(log)"

# How may seperate file folders used to seperate all fast5 files to # of groups
# ${targetNum}=N, 1-N: will create 0 -- (N-1) folder names

if [ "${multipleInputs}" = true ] ; then
	echo "Multiple fast5.tar need to be untared"

	filelist=$(find ${inputDataDir} -type f \( -iname \*.fast5.tar -o -iname \*.fast5.tar.gz -o -iname \*.fast5 \))

	for fast5tarfn in ${filelist}; do
		echo "fn=${fast5tarfn}"
		if [ "${fast5tarfn##*.}" = "tar" ]; then
			tar -xf ${fast5tarfn} -C ${untaredInputDir} & # bg running
        elif [ "${fast5tarfn##*.}" = "gz" ]; then
			tar -xzf ${fast5tarfn} -C ${untaredInputDir} & # bg running
        elif [ "${fast5tarfn##*.}" = "fast5" ]; then
			cp ${fast5tarfn} ${untaredInputDir} & # bg running
        fi
	done
	wait # wait until all background task finished
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


