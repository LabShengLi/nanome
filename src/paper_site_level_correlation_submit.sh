#!/bin/bash

# bash paper_site_level_correlation_submit.sh
# bash paper_site_level_correlation_submit.sh "HL60 K562"
set -x

#Running datasets "HL60 K562 APL NA19240 NA12878"
datasetList=${1:-"HL60 K562"}

if [[ "$datasetList" == *"HL60"* ]]; then
    echo "Running HL60"
    ## include METEORE
    sbatch --job-name=meth.corr.HL60_METEORE \
        paper_site_level_correlation.sbatch HL60 RRBS_2Reps_METEORE METEORE
fi

if [[ "$datasetList" == *"K562"* ]]; then
    echo "Running K562"
    ### METEORE include
    sbatch --job-name=meth.corr.K562_METEORE \
        paper_site_level_correlation.sbatch K562 WGBS_2Reps_METEORE METEORE
fi

if [[ "$datasetList" == *"APL"* ]]; then
    echo "Running APL"
    ### METEORE include
    sbatch --job-name=meth.corr.APL_METEORE \
        paper_site_level_correlation.sbatch APL WGBS_METEORE METEORE
fi

if [[ "$datasetList" == *"NA19240"* ]]; then
    echo "Running NA19240"
    ### METEORE include
    sbatch --job-name=meth.corr.NA19240_METEORE \
        paper_site_level_correlation.sbatch NA19240 RRBS_2Reps_METEORE METEORE
fi

if [[ "$datasetList" == *"NA12878"* ]]; then
    echo "Running NA12878"
    ## run NA12878
    #  sbatch --job-name=meth.corr.NA12878 \
    #    paper_site_level_correlation.sbatch NA12878 WGBS_2Reps NA12878 /projects/li-lab/yang/results/2021-07-05
    sbatch --job-name=meth.corr.NA12878_seven \
        paper_site_level_correlation.sbatch NA12878 WGBS_2Reps_seven NA12878-1 /projects/li-lab/yang/results/2021-07-07
fi

exit 0

## six tool
#sbatch  --job-name=meth.corr.HL60 \
#	paper_site_level_correlation.sbatch HL60 RRBS_2Reps six
#
#sbatch  --job-name=meth.corr.K562 \
#	paper_site_level_correlation.sbatch K562 WGBS_2Reps six
#
#sbatch  --job-name=meth.corr.APL \
#	paper_site_level_correlation.sbatch APL WGBS six
#
#sbatch  --job-name=meth.corr.NA19240 \
#	paper_site_level_correlation.sbatch NA19240 RRBS_2Reps six

#exit 0
