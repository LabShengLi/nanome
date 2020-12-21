#!/bin/bash
#SBATCH --job-name=tsv.combine
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

cmd=$1
output_dir=$2

if [ $cmd = "tombo-add-seq" ]; then
   cat ${output_dir}/*tombo*.tsv > ${output_dir}/K562.tombo.combined.tsv
else
   cat ${output_dir}/*deepmod*.tsv > ${output_dir}/K562.deepmod.combined.tsv
fi
