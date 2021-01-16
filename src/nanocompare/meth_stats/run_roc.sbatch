#!/bin/bash
#SBATCH --job-name=run.roc
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 100g # memory pool for all cores
#SBATCH -t 8:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.figures # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#cd $PBS_O_WORKDIR # change the drectory for the one, from which the job was submitted

set -x

baseDir=/projects/liuya/workspace/tcgajax/nanocompare/meth_stats
pythonFn=roc_compute_and_plot.py

prjBaseDir=/projects/li-lab/nmf_epihet/tcgajax

PYTHONPATH=${prjBaseDir} python ${baseDir}/${pythonFn}


