#!/bin/bash
#SBATCH --job-name=bedtools.closest
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

# Caculate the nearest distance of each tool with bg-truth bed files by bedtools's closest functions

set -x

RunPrefix=${1:-K562_WGBS_Joined}
baseFormat=0

data_base_dir=/projects/li-lab/yang/results/$(date +%F)/${RunPrefix}

bgtruth_fn=$(ls ${data_base_dir}/${RunPrefix}*-meth-cov-bgtruth*-baseCount${baseFormat}.bed)

mkdir -p results
mkdir -p results/${RunPrefix}

outdir=results/${RunPrefix}

echo "Processing file:${bgtruth_fn}"
bedtools sort -i ${bgtruth_fn} > ${outdir}/fb-sorted.bed

i=1
for fn in ${data_base_dir}/${RunPrefix}*-meth-cov-tool-*-baseCount${baseFormat}.bed;
do
	basefn=$(basename ${fn})
	echo "Processing file:$basefn"

	outfn=${basefn/-baseCount0.bed/}-bgtruth-closest.bed

	bedtools sort -i ${data_base_dir}/${basefn} > ${outdir}/f${i}-sorted.bed
	bedtools closest -a ${outdir}/f${i}-sorted.bed -b ${outdir}/fb-sorted.bed -d -s > ${outdir}/${outfn}

	python density_plot.py plot -i ${outdir}/${outfn}
	echo "Save results to ${outdir}/${outfn}"

	i=$((i+1))
done

echo "Closest task is finished"