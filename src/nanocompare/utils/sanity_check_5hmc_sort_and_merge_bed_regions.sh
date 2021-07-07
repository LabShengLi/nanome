#!/bin/bash
#SBATCH --job-name=santi.5hmc.sort.merge
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300G # memory pool for all cores
#SBATCH -t 72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

baseDir=/projects/li-lab/yang/results/2021-07-04

outDir=$baseDir/regions_sorted_merged
mkdir -p $outDir

regionFileList=$(find $baseDir -maxdepth 1 -name "APL.bgtruth.*.regions.*.*.*.bed.gz")

for regionFile in ${regionFileList} ; do
  basenameRegion=$(basename $regionFile)
  bedtools sort -i ${regionFile} | bedtools merge -s -i stdin | gzip > $outDir/${basenameRegion/.bed.gz}.sorted.merged.bed.gz &
done
wait
echo "### sort and merged DONE"

