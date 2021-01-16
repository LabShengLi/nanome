#!/bin/bash
#SBATCH --job-name=tsv.combine
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 03:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.figures # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

cmd=$1
inputfn=$2
output_dir=$3
ntask=$4

dsname=$(basename ${inputfn})
dsname=$(echo "${dsname%%.*}")

if [ $cmd = "tombo-add-seq" ]; then
	# input file '/fastscratch/liuya/nanocompare/K562-Tombo-N50/K562-Tombo-N50-meth-call/K562.tombo.perReadsStats.combine.bed'
	# batch output file 'K562.tombo.perReadsStats.combine-with-seq-info-n300-t001-chr1.tsv'
	# outfn: ${dsname}.tombo.perReadsStatsOnlyCG.combine.tsv
	outfn=${output_dir}/${dsname}.tombo.perReadsStatsOnlyCG.combine.tsv
	rm -rf ${outfn}
    cat ${output_dir}/${dsname}.tombo.perReadsStats.combine-with-seq-info-n${ntask}*.tsv > ${outfn}
    wc -l ${inputfn}
    wc -l ${outfn}

else
	# input file: K562.deepmod.C.combined.bed
	# batch output file: 'K562.deepmod.C.combined-with-seq-info-n300-t001-chr1.tsv'
	# outfn: ${dsname}.deepmod.OnlyCG.combine.tsv
	outfn=${output_dir}/${dsname}.deepmod.OnlyCG.combine.tsv
	rm -rf ${outfn}
    cat ${output_dir}/${dsname}.deepmod.C.combine-with-seq-info-n${ntask}*.tsv > ${outfn}
    wc -l ${inputfn}
    wc -l ${outfn}
fi
