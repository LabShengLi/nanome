#!/bin/bash
#SBATCH --job-name=basecall
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=250g # memory pool for all cores
#SBATCH --time=2-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# Basecalling workflow
# Need to populate the parameters into this script
################################################################################
cd "$(dirname "$0")"

set +x
#source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh
#source /home/liuya/.bash_profile
#export PATH=/cm/shared/apps/slurm/18.08.8/bin:${PATH}
source utils.common.sh
set -x

#processors=8

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
echo "jobkSeptInputDir: ${jobkSeptInputDir}"
echo "jobkBasecallOutputDir: ${jobkBasecallOutputDir}"
echo "processors: ${processors}"
echo "##################"
set +u

mkdir -p ${jobkBasecallOutputDir}

set +x
conda activate nanoai
set -x

## Run Basecalling with Albacore:
time read_fast5_basecaller.py -o fastq,fast5 -t ${processors} \
		-s ${jobkBasecallOutputDir} -i ${jobkSeptInputDir} -c r94_450bps_linear.cfg \
		-n 20000000000

set +x
conda deactivate
set -x

echo "###   Basecalling DONE    ###"