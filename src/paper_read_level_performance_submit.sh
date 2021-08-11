#!/bin/bash

#set -x

#bash paper_read_level_performance_submit.sh "HL60 K562 APL NA19240 NA12878"

#Running datasets "HL60 K562 APL NA19240 NA12878"
datasetList=${1:-"HL60 K562"}

if [[ "$datasetList" == *"HL60"* ]]; then
    echo "Running HL60"
    sbatch --job-name=meth.perf.HL60_seven \
        paper_read_level_performance.sbatch HL60 RRBS_2Reps seven --mpi
fi

if [[ "$datasetList" == *"K562"* ]]; then
    echo "Running K562"
    sbatch --job-name=meth.perf.K562_seven \
        paper_read_level_performance.sbatch K562 WGBS_2Reps seven --mpi
fi

if [[ "$datasetList" == *"APL"* ]]; then
    echo "Running APL"
    sbatch --job-name=meth.perf.APL_seven \
        paper_read_level_performance.sbatch APL WGBS seven
fi

if [[ "$datasetList" == *"NA19240"* ]]; then
    echo "Running NA19240"
    sbatch --job-name=meth.perf.NA19240_seven \
        paper_read_level_performance.sbatch NA19240 RRBS_2Reps seven
fi

if [[ "$datasetList" == *"NA12878"* ]]; then
    echo "Running NA12878"
    sbatch --job-name=meth.perf.NA12878_seven \
        paper_read_level_performance.sbatch NA12878 WGBS_2Reps seven
fi
