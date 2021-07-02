#!/bin/bash

set -x

### METEORE include
sbatch  --job-name=meth.perf.HL60_METEORE \
	paper_read_level_performance.sbatch HL60 RRBS_2Reps_METEORE METEORE

sbatch  --job-name=meth.perf.K562_METEORE \
	paper_read_level_performance.sbatch K562 WGBS_2Reps_METEORE METEORE

sbatch  --job-name=meth.perf.APL_METEORE \
	paper_read_level_performance.sbatch APL WGBS_METEORE METEORE "  "

sbatch  --job-name=meth.perf.NA19240_METEORE \
	paper_read_level_performance.sbatch NA19240 RRBS_2Reps_METEORE METEORE "  "

exit 0

### six tool
sbatch  --job-name=meth.perf.HL60 \
	paper_read_level_performance.sbatch HL60 RRBS_2Reps six

sbatch  --job-name=meth.perf.K562 \
	paper_read_level_performance.sbatch K562 WGBS_2Reps six

sbatch  --job-name=meth.perf.APL \
	paper_read_level_performance.sbatch APL WGBS six "   "

sbatch  --job-name=meth.perf.NA19240 \
	paper_read_level_performance.sbatch NA19240 RRBS_2Reps six "   "

exit 0


