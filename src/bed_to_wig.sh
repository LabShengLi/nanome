#!/bin/bash
#SBATCH --job-name=cov.bed.to.bw
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem 150G # memory pool for all cores
#SBATCH -t 20:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

bedGraphToBigWig=/projects/li-lab/yang/scripts/bedGraphToBigWig

hg38ChromSize=/projects/li-lab/yang/workspace/nano-compare/data/reference/hg38/hg38.chrom.sizes
basedir=/projects/li-lab/yang/results/2021-04-23

#dsnamelist=(NA19240 APL K562 HL60)
dsnamelist=(K562 HL60)

for dsname in "${dsnamelist[@]}"; do
  echo dsname=$dsname
  # NA19240_RRBS_2Reps.tss.bgtruth.cov5.sorted.bedGraph
  flist=$(find $basedir -name "${dsname}*.sorted.bedGraph")

  for fn in $flist; do
    outfn=${fn/".sorted.bedGraph"/".bw"}
    echo fn=$fn, outfn=$outfn
    # Convert bedgraph file into bigwig file for deeptools plotting
    ${bedGraphToBigWig} $fn ${hg38ChromSize} $outfn
  done

done
