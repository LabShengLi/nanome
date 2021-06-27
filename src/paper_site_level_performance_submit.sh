#!/bin/bash

set -x

## six tool
sbatch  --job-name=meth.corr.HL60 \
	paper_site_level_correlation.sbatch HL60 RRBS_2Reps

sbatch  --job-name=meth.corr.K562 \
	paper_site_level_correlation.sbatch K562 WGBS_2Reps

sbatch  --job-name=meth.corr.APL \
	paper_site_level_correlation.sbatch APL WGBS

sbatch  --job-name=meth.corr.NA19240 \
	paper_site_level_correlation.sbatch NA19240 RRBS_2Reps

#exit 0

## include METEORE
sbatch  --job-name=meth.corr.HL60_METEORE \
	paper_site_level_correlation.sbatch HL60 RRBS_2Reps_METEORE METEORE

sbatch  --job-name=meth.corr.K562_METEORE \
	paper_site_level_correlation.sbatch K562 WGBS_2Reps_METEORE METEORE

sbatch  --job-name=meth.corr.APL_METEORE \
	paper_site_level_correlation.sbatch APL WGBS_METEORE METEORE

sbatch  --job-name=meth.corr.NA19240_METEORE \
	paper_site_level_correlation.sbatch NA19240 RRBS_2Reps_METEORE METEORE

#exit 0

