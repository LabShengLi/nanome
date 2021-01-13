#!/bin/bash
#SBATCH --job-name=basecall
##SBATCH --partition=compute
##SBATCH -N 1 # number of nodes
##SBATCH -n 16 # number of cores
#SBATCH -p gpu
#SBATCH --gres=gpu:3             # number of gpus per node
#SBATCH -q inference
#SBATCH --mem=125g # memory pool for all cores
#SBATCH --time=06:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
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
jobkSeptInputDir=${septInputDir}/${job_index}
jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "targetNum: ${targetNum}"
echo "analysisPrefix: ${analysisPrefix}"
echo "septInputDir: ${septInputDir}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "basecall_name: ${basecall_name}"
echo "jobkSeptInputDir: ${jobkSeptInputDir}"
echo "jobkBasecallOutputDir: ${jobkBasecallOutputDir}"
echo "processors: ${processors}"
echo "##################"
set +u

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
#	echo "Not implemented yet"
fi

set +x
conda deactivate
set -x

echo "###   Basecalling DONE    ###"