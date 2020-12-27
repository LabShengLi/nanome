#!/bin/bash
#SBATCH --job-name=tsv.combine
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 13:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

cmd=$1
inputfn=$2
output_dir=$3
ntask=$4

dsname=$(basename ${inputfn})
dsname=$(echo "${dsname%%.*}")

if [ $cmd = "tombo-add-seq" ]; then
	rm -rf ${output_dir}/${dsname}.tombo.combined.tsv
    cat ${output_dir}/${dsname}*tombo_*-n${ntask}*.tsv > ${output_dir}/${dsname}.tombo.combined.onlycpg.tsv
    wc -l ${inputfn}
    wc -l ${output_dir}/${dsname}.tombo.combined.onlycpg.tsv

else
	rm -rf ${output_dir}/${dsname}.deepmod.combined.tsv
    cat ${output_dir}/${dsname}*deepmod_*-n${ntask}*.tsv > ${output_dir}/${dsname}.deepmod.combined.onlycpg.tsv
    wc -l ${inputfn}
    wc -l ${output_dir}/${dsname}.deepmod.combined.onlycpg.tsv
fi
