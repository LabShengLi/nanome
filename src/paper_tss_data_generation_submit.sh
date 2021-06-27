#!/bin/bash

set -x

sbatch  --job-name=tss.HL60 \
	paper_tss_data_generation.sbatch HL60

sbatch  --job-name=tss.K562 \
	paper_tss_data_generation.sbatch K562

sbatch  --job-name=tss.APL \
	paper_tss_data_generation.sbatch APL

sbatch  --job-name=tss.NA19240 \
	paper_tss_data_generation.sbatch NA19240
