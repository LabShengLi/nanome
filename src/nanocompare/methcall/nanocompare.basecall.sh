#!/bin/bash
#SBATCH --job-name=basecall
##SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH -p gpu
##SBATCH --gres=gpu:1
##SBATCH -q inference
##SBATCH --mem=100g
#SBATCH --time=01:00:00
#SBATCH -o log/%x.%j.out
#SBATCH -e log/%x.%j.err
##SBATCH --array=1-11

################################################################################
# Basecalling workflow
# Need to populate the parameters into this script
################################################################################
# Working dir now must be at file location now

set +x
source ../../utils.common.sh
set -x


job_index=$((SLURM_ARRAY_TASK_ID))
jobkSeptInputDir=${septInputDir}/M${job_index}
jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "targetNum: ${targetNum}"
echo "basecall_name: ${basecall_name}"
echo "septInputDir: ${septInputDir}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "jobkSeptInputDir: ${jobkSeptInputDir}"
echo "jobkBasecallOutputDir: ${jobkBasecallOutputDir}"
echo "processors: ${processors}"
echo "##################"
set +u

rm -rf ${jobkBasecallOutputDir}/*

mkdir -p ${jobkBasecallOutputDir}

set +x
conda activate nanoai
set -x

if [ "${basecall_name}" = "Albacore" ] ; then
	## Run Basecalling with Albacore:
	time read_fast5_basecaller.py -o fastq,fast5 -t ${processors} \
			-s ${jobkBasecallOutputDir} -i ${jobkSeptInputDir} -c r94_450bps_linear.cfg \
			-n 20000000000

elif [ "${basecall_name}" = "Guppy" ] ; then
	## Run Basecalling with Guppy:
	time ${GuppyDir}/bin/guppy_basecaller --input_path ${jobkSeptInputDir} \
	    --save_path ${jobkBasecallOutputDir} --config dna_r9.4.1_450bps_hac.cfg \
	    --gpu_runners_per_device ${processors} --num_callers 3 --fast5_out --verbose_logs --device auto
fi

set +x
conda deactivate
set -x

echo "###   Basecalling DONE    ###"