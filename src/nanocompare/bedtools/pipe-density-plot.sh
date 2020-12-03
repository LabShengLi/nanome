#!/bin/bash

set -x

plot_taskid=$(bash /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/Methylation_correlation_plotting_submit.sh)
plot_taskid=$(echo ${plot_taskid} |grep -Eo '[0-9]+$')


bedtool_taskid=$(sbatch --dependency=afterok:${plot_taskid} closest_bedtools.sh)
bedtool_taskid=$(echo ${bedtool_taskid} |grep -Eo '[0-9]+$')

sbatch --dependency=afterok:${bedtool_taskid} density_plot.sh
