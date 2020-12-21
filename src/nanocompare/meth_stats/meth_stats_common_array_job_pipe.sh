#!/bin/bash

# Usage:
# bash meth_stats_common_array_job_pipe.sh tombo-add-seq 350
# bash meth_stats_common_array_job_pipe.sh deepmod-add-seq 100
# $1 - command
# $2 - num-tasks

#set -x

# Command name
cmd=${1:-tombo-add-seq}

# Number of tasks
N=${2:-350}

# Input file of original Tombo or DeepMod
inputfn=${3:-/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/K562/K562.tombo_perReadsStats.bed}

output_dir=/projects/li-lab/yang/results/$(date +%F)

echo cmd=${cmd}, N=$N, output_dir=${output_dir}

#################################################
### Task 1: calculate tombo sequence using job array
#################################################
job_array_ret=$(sbatch --job-name ${cmd}.arrjob --array=1-$N meth_stats_common_array_job.sh ${cmd} ${inputfn})
job_array_id=${job_array_ret##* }

echo submit ${cmd} array job ok, jobid=${job_array_id}, num of jobs=$N

ARRAY_JOB_LIST=

for ((i=1;i<=$N;i++));
do
   ARRAY_JOB_LIST=${ARRAY_JOB_LIST}:${job_array_id}_$i
done

#echo ${ARRAY_JOB_LIST}

#################################################
### Task 2: combine tombo results together
#################################################
sbatch --job-name combine.${cmd} --dependency=afterok${ARRAY_JOB_LIST} meth_stats_common_array_job_combine.sh ${cmd} ${output_dir}

echo submit combine results job ok