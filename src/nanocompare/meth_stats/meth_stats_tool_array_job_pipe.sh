#!/bin/bash

# Usage:
# bash meth_stats_tool_array_job_pipe.sh tombo-add-seq 350 <input-file> <outdir>
# bash meth_stats_tool_array_job_pipe.sh deepmod-add-seq 100 <input-file> <outdir>
# $1 - command: tombo-add-seq, deepmod-add-seq
# $2 - num-tasks
# $3 - input file
# $4 - output dir

set -x

export PATH=/cm/shared/apps/slurm/18.08.8/bin:${PATH}
sbatch --version

# Command name
cmd=${1:-tombo-add-seq}

# Number of tasks
num_jobs=${2:-350}

# Input file of original Tombo or DeepMod
inputfn=${3:-/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/K562/K562.tombo_perReadsStats.bed}

output_dir=${4:-/projects/li-lab/yang/results/$(date +%F)}

dsname=$(basename ${inputfn})
dsname=$(echo "${dsname%%.*}")

mkdir -p ${output_dir}
mkdir -p ${output_dir}/log

echo cmd=${cmd}, N=${num_jobs}, output_dir=${output_dir}

#################################################
### Task 1: calculate tombo sequence using job array
#################################################
job_array_ret=$(sbatch --job-name=${dsname}.${cmd} --output=${output_dir}/log/%x.batch%a.%j.out --error=${output_dir}/log/%x.batch%a.%j.err --array=1-${num_jobs} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/meth_stats_tool_array_job.sh ${cmd} ${inputfn} ${output_dir})
job_array_id=${job_array_ret##* }

echo submit ${cmd} array job ok, jobid=${job_array_id}, num of jobs=$num_jobs

ARRAY_JOB_LIST=

for ((i=1;i<=$num_jobs;i++));
do
   ARRAY_JOB_LIST=${ARRAY_JOB_LIST}:${job_array_id}_$i
done

#echo ${ARRAY_JOB_LIST}

#################################################
### Task 2: combine tombo results together
#################################################
sbatch --job-name combine.${cmd} --output=${output_dir}/log/%x.%j.out --error=${output_dir}/log/%x.%j.err --dependency=afterok${ARRAY_JOB_LIST} /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/meth_stats_tool_array_job_combine.sh ${cmd} ${inputfn} ${output_dir} $num_jobs

echo submit combine results job ok