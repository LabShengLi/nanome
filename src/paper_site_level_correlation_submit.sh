#!/bin/bash

# bash paper_site_level_correlation_submit.sh "HL60 K562 APL NA19240 NA12878"
# Running datasets "HL60 K562 APL NA19240 NA12878"

datasetList=${1:-"HL60 K562"}
bedDir=${2:-"/projects/li-lab/yang/results/2021-07-07"} # folder contained MethPerf concordant and discordent bed files

if [[ "$datasetList" == *"HL60"* ]]; then
    echo "Running HL60"
    sbatch --job-name=meth.corr.HL60_seven \
        paper_site_level_correlation.sbatch HL60 RRBS_2Reps seven ${bedDir}
fi

if [[ "$datasetList" == *"K562"* ]]; then
    echo "Running K562"
    sbatch --job-name=meth.corr.K562_seven \
        paper_site_level_correlation.sbatch K562 WGBS_2Reps seven ${bedDir}
fi

if [[ "$datasetList" == *"APL"* ]]; then
    echo "Running APL"
    sbatch --job-name=meth.corr.APL_seven \
        paper_site_level_correlation.sbatch APL WGBS seven ${bedDir}
fi

if [[ "$datasetList" == *"NA19240"* ]]; then
    echo "Running NA19240"
    sbatch --job-name=meth.corr.NA19240_seven \
        paper_site_level_correlation.sbatch NA19240 RRBS_2Reps seven ${bedDir}
fi

if [[ "$datasetList" == *"NA12878"* ]]; then
    echo "Running NA12878"
    sbatch --job-name=meth.corr.NA12878_seven \
        paper_site_level_correlation.sbatch NA12878 WGBS_2Reps seven ${bedDir}
fi
