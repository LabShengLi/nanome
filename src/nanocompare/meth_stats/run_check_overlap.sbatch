#!/bin/bash
#SBATCH --job-name=overlap
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 1-18:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

baseDir=/projects/liuya/workspace/tcgajax/nanocompare/meth_stats
pythonFn=check_overlap_four_tools.py

prjBaseDir=/projects/li-lab/nmf_epihet/tcgajax

PYTHONPATH=${prjBaseDir} python ${baseDir}/${pythonFn}
