#!/bin/bash
#SBATCH --job-name=report.nanocompare
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 2-18:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.figures # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

baseDir=/projects/liuya/workspace/tcgajax/nanocompare/meth_stats
pythonFn=report_generation.py

prjBaseDir=/projects/li-lab/nmf_epihet/tcgajax

PYTHONPATH=${prjBaseDir} python ${baseDir}/${pythonFn}
