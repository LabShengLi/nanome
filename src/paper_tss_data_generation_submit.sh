#!/bin/bash

set -x

## six tool
sbatch  --job-name=tss.paper.HL60 \
	paper_tss_data_generation.sbatch HL60

exit 0

sbatch  --job-name=tss.paper.K562 \
	paper_tss_data_generation.sbatch K562

sbatch  --job-name=tss.paper.APL \
	paper_tss_data_generation.sbatch APL

sbatch  --job-name=tss.paper.NA19240 \
	paper_tss_data_generation.sbatch NA19240
