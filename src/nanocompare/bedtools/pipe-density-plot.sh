#!/bin/bash

set -x

# Submit meth-corr tasks for each row of config file
meth_corr_task_ret=$(bash /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/Methylation_correlation_plotting_submit.sh)

echo Submit the meth-corr tasks.

# For each task
while IFS= read -r sepline; do
	read -r line1
	read -r line2

    RunPrefix=${line1[0]##*=}
	meth_corr_taskid=$(echo ${line2} |grep -Eo '[0-9]+$')
	echo RunPrefix=$RunPrefix, taskid=${meth_corr_taskid}, and submit the closest analysis task.

	# Submit bedtool closest analysis for each meth-corr task
	bedtool_taskid=$(sbatch --dependency=afterok:${meth_corr_taskid} closest_bedtools.sh ${RunPrefix})

#    read -r sepline
done <<< "$meth_corr_task_ret"





