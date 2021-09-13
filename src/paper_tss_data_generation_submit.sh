#!/bin/bash

# bash paper_tss_data_generation_submit.sh NA12878
# bash paper_tss_data_generation_submit.sh "HL60 K562 APL NA19240 NA12878"
# Running datasets "HL60 K562 APL NA19240 NA12878"
datasetList=${1:-"HL60 K562"}

if [[ "$datasetList" == *"HL60"* ]]; then
  echo "Running HL60"
  sbatch --job-name=tss.HL60 \
    paper_tss_data_generation.sbatch HL60
fi

if [[ "$datasetList" == *"K562"* ]]; then
  echo "Running K562"
  sbatch --job-name=tss.K562 \
    paper_tss_data_generation.sbatch K562
fi

if [[ "$datasetList" == *"APL"* ]]; then
  echo "Running APL"
  sbatch --job-name=tss.APL \
    paper_tss_data_generation.sbatch APL
fi

if [[ "$datasetList" == *"NA19240"* ]]; then
  echo "Running NA19240"
  sbatch --job-name=tss.NA19240 \
    paper_tss_data_generation.sbatch NA19240
fi

if [[ "$datasetList" == *"NA12878"* ]]; then
  echo "Running NA12878"
  sbatch --job-name=tss.NA12878 \
    paper_tss_data_generation.sbatch NA12878
fi
