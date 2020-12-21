#!/bin/bash
#SBATCH --job-name=meth-common
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem 100g # memory pool for all cores
#SBATCH -t 05:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
#SBATCH --array=1-200 # job array index

set -x

cmd=$1

prj_dir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${prj_dir}/src/nanocompare/meth_stats/meth_stats_common.py

mkdir -p log

#python ${pythonFile} $@

python ${pythonFile} ${cmd} -n ${SLURM_ARRAY_TASK_COUNT} -t ${SLURM_ARRAY_TASK_ID}
