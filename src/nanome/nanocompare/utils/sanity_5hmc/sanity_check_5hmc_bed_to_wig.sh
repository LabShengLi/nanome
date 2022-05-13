#!/bin/bash
#SBATCH --job-name=santi.5hmc.bed2wig
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300G # memory pool for all cores
#SBATCH -t 72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

## Convert the 5hmC level in bed to BigWig format

set -x

baseDir=/projects/li-lab/yang/results/2021-07-04
flist=$(find $baseDir -name "APL.5hmc.tss.cov*.bed.gz")

bedGraphToBigWig=/projects/li-lab/yang/scripts/bedGraphToBigWig
hg38ChromSize=/projects/li-lab/yang/workspace/nano-compare/data/reference/hg38/hg38.chrom.sizes

for fn in ${flist}; do
	s1=${fn%.bed.gz}
	cutoff_num=${s1#*APL.5hmc.tss.cov}
	python ../../plot_figure.py bed-to-bedGraph -i $fn --cutoff ${cutoff_num} --cov-col -1  --signal-col 3
done

flist=$(find $baseDir -name "*.sorted.bedGraph")
for fn in ${flist}; do
	outfn=${fn/".sorted.bedGraph"/".bw"}
	echo fn=$fn, outfn=$outfn
	# Convert bedgraph file into bigwig file for deeptools plotting
	${bedGraphToBigWig} $fn ${hg38ChromSize} $outfn
done

rm -f ${baseDir}/*.bedGraph
