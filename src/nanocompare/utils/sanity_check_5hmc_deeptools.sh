#!/bin/bash
#SBATCH --job-name=santi.5hmc.deeptools
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300G # memory pool for all cores
#SBATCH -t 72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x
set -e

flank=100
binSize=10
processors=32

baseDir=/projects/li-lab/yang/results/2021-07-07
baseDir=/projects/li-lab/Nanopore_compare/suppdata/5hmc-sanity

# Output dir
outdir=/projects/li-lab/yang/results/$(date +%F)/sanity-check-5hmc
mkdir -p $outdir

bw_5hmc=$(find $baseDir -maxdepth 1 -name "APL.5hmc.tss.cov1.bw")
regionList=$(find $baseDir/regions_sorted_merged -maxdepth 1 -name "APL.bgtruth.*.regions.*.bed")

for regionFile in ${regionList}; do
  echo "Process regionFile=${regionFile}"
  basenameRegion=$(basename $regionFile)

  ### Generate the matrix for heatmap (REFERENCE-Point):
  computeMatrix reference-point -p ${processors} \
    --referencePoint center \
    -S $bw_5hmc \
    --samplesLabel 5hmC \
    -R $regionFile \
    --beforeRegionStartLength $flank \
    --afterRegionStartLength $flank \
    --binSize $binSize \
    --skipZeros \
    -o $outdir/${basenameRegion/.bed.gz/}.deeptools.flank${flank}.bin${binSize}.mat.gz

  #### Profile plot:
  plotProfile -m $outdir/${basenameRegion/.bed.gz/}.deeptools.flank${flank}.bin${binSize}.mat.gz \
    -out $outdir/${basenameRegion/.bed.gz/}.deeptools.flank${flank}.bin${binSize}.5hmc.profileplot.png \
    --plotType=lines \
    --perGroup \
    --colors red \
    --plotTitle 5hmC \
    --yMin 0 \
    --yMax 1

  #### Generate the heatmap:
  plotHeatmap --colorMap Reds \
    --missingDataColor '#ffffff' \
    -m $outdir/${basenameRegion/.bed.gz/}.deeptools.flank${flank}.bin${binSize}.mat.gz \
    --whatToShow 'plot, heatmap and colorbar' \
    --refPointLabel "5hmC" \
    --xAxisLabel "dist. to Region" \
    --outFileSortedRegions $outdir/${basenameRegion/.bed.gz/}.deeptools.flank${flank}.bin${binSize}.sortedRegions.txt \
    --outFileNameMatrix $outdir/${basenameRegion/.bed.gz/}.deeptools.flank${flank}.bin${binSize}.outMatrix.tsv.gz \
    --yMin 0 \
    --yMax 1 \
    -out $outdir/${basenameRegion/.bed.gz/}.deeptools.flank${flank}.bin${binSize}.5hmc.plotheatmap.TSS.png
done
