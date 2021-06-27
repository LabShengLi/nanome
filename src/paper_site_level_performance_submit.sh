#!/bin/bash

set -x

sbatch  --job-name=meth.corr.paper.HL60 \
	paper_site_level_correlation.sbatch HL60 RRBS_2Reps

sbatch  --job-name=meth.corr.paper.K562 \
	paper_site_level_correlation.sbatch K562 WGBS_2Reps

sbatch  --job-name=meth.corr.paper.APL \
	paper_site_level_correlation.sbatch APL WGBS

sbatch  --job-name=meth.corr.paper.NA19240 \
	paper_site_level_correlation.sbatch NA19240 RRBS_2Reps
