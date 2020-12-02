#!/bin/bash

set -x

bedtool_taskid=$(sbatch closest_bedtools.sh)

bedtool_taskid=$(echo ${bedtool_taskid} |grep -Eo '[0-9]+$')

sbatch --dependency=afterok:${bedtool_taskid} density_plot.sh
