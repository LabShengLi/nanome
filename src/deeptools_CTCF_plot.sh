#!/bin/bash
#SBATCH --job-name=plot.CTCF
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem 150G # memory pool for all cores
#SBATCH -t 05:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x
set -e

# Input dir for bw files
basedir=/projects/li-lab/yang/results/2021-04-23

# Region input dir
regionDir=/fastscratch/c-panz/CTCFpeak

# Output dir
outdir=/projects/li-lab/yang/results/$(date +%F)/CTCF-plots
mkdir -p $outdir

# Get input params
dsname=${1:-NA19240}
actions=${2:-"compute plot"}
binSize=${3:-50}

# Plotting configuration
flank=2000

regionFile=$(find ${regionDir} -name "${dsname}*.bed")

# Get each bw file
DeepSignal=$(find $basedir -name "${dsname}*DeepSignal.cov3.bw")
Tombo=$(find $basedir -name "${dsname}*Tombo.cov3.bw")
Nanopolish=$(find $basedir -name "${dsname}*Nanopolish.cov3.bw")
DeepMod=$(find $basedir -name "${dsname}*DeepMod.cov3.bw")
Megalodon=$(find $basedir -name "${dsname}*Megalodon.cov3.bw")
BGTruth=$(find $basedir -name "${dsname}*bgtruth.cov5.bw")

if [ "$dsname" == "K562" ] || [ "$dsname" == "APL" ]; then
    BSLabel="WGBS"
else
    BSLabel="RRBS"
fi

# Compute matrix if needed
if [[ " ${actions} " =~ " compute " ]]; then
    ### Generate the matrix for heatmap (REFERENCE-Point):
    computeMatrix reference-point -p 16 \
        --referencePoint center \
        -S $BGTruth $DeepSignal $Tombo $Nanopolish $DeepMod $Megalodon \
        --samplesLabel ${BSLabel} DeepSignal Tombo Nanopolish DeepMod Megalodon \
        -R $regionFile \
        --beforeRegionStartLength $flank \
        --afterRegionStartLength $flank \
        --binSize $binSize \
        --skipZeros \
        -o $outdir/$dsname.bin${binSize}.flank${flank}.CTCF.mat.gz
fi

#### Profile plot:
plotProfile -m $outdir/$dsname.bin${binSize}.flank${flank}.CTCF.mat.gz \
    -out $outdir/$dsname.bin${binSize}.profileplot.CTCF.png \
    --plotType=lines \
    --perGroup \
    --colors black "#999999" "#E69F00" "#56B4E9" "#009E73" "#CC79A7" "#0072B2" \
    --plotTitle $dsname \
    --yMin 0 \
    --yMax 1

#### Generate the heatmap:
plotHeatmap --colorMap Reds \
    --missingDataColor '#ffffff' \
    -m $outdir/$dsname.bin${binSize}.flank${flank}.CTCF.mat.gz \
    --whatToShow 'plot, heatmap and colorbar' \
    --refPointLabel "CTCF" \
    --xAxisLabel "dist. to CTCF" \
    --outFileSortedRegions $outdir/$dsname.bin${binSize}.flank${flank}.CTCF.sortedRegions.txt \
    --outFileNameMatrix $outdir/$dsname.bin${binSize}.flank${flank}.CTCF.outMatrix.tsv.gz \
    --yMin 0 \
    --yMax 1 \
    -out $outdir/$dsname.bin${binSize}.plotheatmap.CTCF.png

#### Run as: bash script.sh NA19240 "compute plot"
