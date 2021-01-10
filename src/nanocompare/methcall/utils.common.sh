#!/bin/bash

# Specify your home bash, this required by conda setup
HOME_BASH_FN=/home/liuya/.bash_profile

source ${HOME_BASH_FN}
export PATH=/cm/shared/apps/slurm/18.08.8/bin:${PATH}

get_arrayjob_ids(){
	## Generate dependency pattern format from sbatch job submission return string
	## Input is <ret-str> <num>
	## Ouput is like: ":12345_1:12345_2:12345_3:12345_4:12345_5:12345_6:12345_7:12345_8:12345_9:12345_10"
	## Example:
	## ret='Submitted batch job 5812809'
	## num=10
	## task_ids=$(get_arrayjob_ids "${ret}" "${num}")
	## Later can be used by "--dependency=afterok${task_ids}"

	usage="get_arrayjob_ids '<ret-str>' '<num>'"
	local ret=${1:?"undefined '<ret-str>': $usage"}
	local num=${2:?"undefined '<num>': $usage"}

	local arrayjob_id=${ret##* }
	local taskids=
	for ((i=1;i<=${num};i++));
	do
		taskids=${taskids}:${arrayjob_id}_$i
	done
	echo ${taskids}
	return 0
}
