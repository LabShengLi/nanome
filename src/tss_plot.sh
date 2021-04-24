#!/bin/bash
#SBATCH --job-name=plot.tss
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem 150G # memory pool for all cores
#SBATCH -t 20:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

# Input dir for bw files
basedir=/projects/li-lab/yang/results/2021-04-23

# Output dir
outdir=/projects/li-lab/yang/results/$(date +%F)/tss-plots
mkdir -p $outdir

# Promoter +/- 2000 region file
regionFile=/pod/2/li-lab/Nanopore_compare/data/genome_annotation/hg38.promoter.2000.bed

# Plotting configuration
flank=2000
binSize=50

# Get input params
dsname=${1:-NA19240}
actions=${2:-"compute plot"}

# Get each bw file
DeepSignal=$(find $basedir -name "${dsname}*DeepSignal.cov3.bw")
Tombo=$(find $basedir -name "${dsname}*Tombo.cov3.bw")
Nanopolish=$(find $basedir -name "${dsname}*Nanopolish.cov3.bw")
DeepMod=$(find $basedir -name "${dsname}*DeepMod.cov3.bw")
Megalodon=$(find $basedir -name "${dsname}*Megalodon.cov3.bw")
BGTruth=$(find $basedir -name "${dsname}*bgtruth.cov5.bw")

# Compute matrix if needed
if [[ " ${actions} " =~ " compute " ]]; then
  ### Generate the matrix for heatmap (REFERENCE-Point):
  computeMatrix reference-point -p 16 \
    --referencePoint center \
    -S $BGTruth $DeepSignal $Tombo $Nanopolish $DeepMod $Megalodon \
    --samplesLabel BS-seq DeepSignal Tombo Nanopolish DeepMod Megalodon \
    -R $regionFile \
    --beforeRegionStartLength $flank \
    --afterRegionStartLength $flank \
    --binSize $binSize \
    --skipZeros \
    -o $outdir/$dsname.mat.gz
fi

#### Profile plot:
plotProfile -m $outdir/$dsname.mat.gz \
  -out $outdir/$dsname.profileplot.png \
  --plotType=lines \
  --perGroup \
  --colors "#D55E00" "#999999" "#E69F00" "#56B4E9" "#009E73" "#CC79A7" "#0072B2" \
  --plotTitle $dsname \
  --yMin 0 \
  --yMax 1

#### Generate the heatmap:
plotHeatmap --colorMap Reds \
  --missingDataColor '#ffffff' \
  -m $outdir/$dsname.mat.gz \
  --whatToShow 'plot, heatmap and colorbar' \
  --refPointLabel "TSS" \
  --xAxisLabel "dist. to TSS" \
  --outFileSortedRegions $outdir/$dsname.sortedRegions.txt \
  --outFileNameMatrix $outdir/$dsname.outMatrix.tsv \
  --yMin 0 \
  --yMax 1 \
  -out $outdir/$dsname.plotheatmap.png

#### Run as: bash script.sh NA19240 "compute plot"
