#!/bin/bash
#SBATCH --job-name=cov.bed.to.bw
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem 150G # memory pool for all cores
#SBATCH -t 20:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#### Convert bed file into bigwig file format
set -x

bedGraphToBigWig=/projects/li-lab/yang/scripts/bedGraphToBigWig

hg38ChromSize=/projects/li-lab/yang/workspace/nano-compare/data/reference/hg38/hg38.chrom.sizes

# Input dir for bed files: such as *.cov1.bed
indir=/projects/li-lab/yang/results/2021-04-22

# Output dir for bw files
outdir=/projects/li-lab/yang/results/$(date +%F)/tss-data
mkdir -p $outdir

dsnamelist=(NA19240 APL K562 HL60)
#dsnamelist=(K562 HL60)

for dsname in "${dsnamelist[@]}"; do
    echo dsname=$dsname

    ### apply cov cutoff and convert BED to bedGraph format, and sorted BED
    flist=$(find ${indir} -name "*.cov1.bed")

    for fn in $flist; do # Not recommended, will break on whitespace
        echo process "$fn"
        bglabel=bgtruth
        if [[ "$fn" == *"$bglabel"* ]]; then
            CUTOFF=5
        else
            CUTOFF=3
        fi
        python plot_figure.py bed-to-bedGraph -i $fn --cutoff $CUTOFF -o $outdir
    done

    # NA19240_RRBS_2Reps.tss.bgtruth.cov5.sorted.bedGraph
    flist=$(find $outdir -name "${dsname}*.sorted.bedGraph")

    for fn in $flist; do
        outfn=${fn/".sorted.bedGraph"/".bw"}
        echo fn=$fn, outfn=$outfn
        # Convert bedgraph file into bigwig file for deeptools plotting
        ${bedGraphToBigWig} $fn ${hg38ChromSize} $outfn
    done

done

## Clean up
rm -rf $outdir/*.bedGraph