#!/bin/bash

#set -x

meth_corr_task_ret=$(bash /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/Methylation_correlation_plotting_submit.sh)

while IFS= read -r line1; do
	read -r line2

    RunPrefix=${line1[0]##*=}
	echo RunPrefix=$RunPrefix

	meth_corr_taskid=$(echo ${line2} |grep -Eo '[0-9]+$')
	echo taskid=${meth_corr_taskid}

	bedtool_taskid=$(sbatch --dependency=afterok:${meth_corr_taskid} closest_bedtools.sh ${RunPrefix})

    read -r sepline
done <<< "$meth_corr_task_ret"




