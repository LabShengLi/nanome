#!/bin/bash

set -x

sbatch  --job-name=meth.perf.paper.HL60 \
	paper_read_level_performance.sbatch HL60 RRBS_2Reps

sbatch  --job-name=meth.perf.paper.K562 \
	paper_read_level_performance.sbatch K562 WGBS_2Reps

sbatch  --job-name=meth.perf.paper.APL \
	paper_read_level_performance.sbatch APL WGBS

sbatch  --job-name=meth.perf.paper.NA19240 \
	paper_read_level_performance.sbatch NA19240 RRBS_2Reps