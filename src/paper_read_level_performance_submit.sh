#!/bin/bash

#set -x

#Running datasets "HL60 K562 APL NA19240 NA12878"
datasetList=${1:-"HL60 K562"}

if [[ "$datasetList" == *"NA12878"* ]]; then
  echo "Running NA12878"
  ### NA12878
  sbatch --job-name=meth.perf.NA12878 \
    paper_read_level_performance.sbatch NA12878 WGBS_2Reps NA12878 "  "
  sbatch --job-name=meth.perf.NA12878_seven \
    paper_read_level_performance.sbatch NA12878 WGBS_2Reps_seven NA12878-1 "  "
fi

if [[ "$datasetList" == *"HL60"* ]]; then
  echo "Running HL60"
  ### METEORE include
  sbatch --job-name=meth.perf.HL60_METEORE \
    paper_read_level_performance.sbatch HL60 RRBS_2Reps_METEORE METEORE
fi

if [[ "$datasetList" == *"K562"* ]]; then
  echo "Running K562"
  ### METEORE include
  sbatch --job-name=meth.perf.K562_METEORE \
    paper_read_level_performance.sbatch K562 WGBS_2Reps_METEORE METEORE
fi

if [[ "$datasetList" == *"APL"* ]]; then
  echo "Running APL"
  ### METEORE include
  sbatch --job-name=meth.perf.APL_METEORE \
    paper_read_level_performance.sbatch APL WGBS_METEORE METEORE "  "
fi

if [[ "$datasetList" == *"NA19240"* ]]; then
  echo "Running NA19240"
  ### METEORE include
  sbatch --job-name=meth.perf.NA19240 \
    paper_read_level_performance.sbatch NA19240 RRBS_2Reps six "   "
fi
