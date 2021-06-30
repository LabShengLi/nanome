#!/bin/bash
#SBATCH --job-name=meth-tool
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300G # memory pool for all cores
#SBATCH -t 72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

# /pod/2/li-lab/yang/workspace/nano-compare/src/nanocompare/utils/

#sbatch /pod/2/li-lab/yang/workspace/nano-compare/src/nanocompare/utils/meth_stats_tool.sh bismark-convert -i /pod/2/li-lab/Ziwei/Nanopore_methyl_compare/result/APL_BSseq/APL-bs_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bed -o .
set -x

NanoCompareDir=/projects/li-lab/yang/workspace/nano-compare
pythonFile=${NanoCompareDir}/src/nanocompare/utils/meth_stats_tool.py

PHTHONPATH=${NanoCompareDir}/src python ${pythonFile} $@
