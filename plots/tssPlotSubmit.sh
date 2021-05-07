#!/bin/bash

set -x
set -e


baseDir=plots
inputDir=result

## Local Rscript running
PATH=/Library/Frameworks/R.framework/Resources/bin:$PATH
Rscript="Rscript"

## Server Rscript in singularity
Rscript="Rscript.sh"

###################################
###################################
## CTCF plot for call R script

## Plot CTCF for NA19240
${Rscript} $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/NA19240.bin100.flank2000.CTCF.outMatrix.tsv.gz \
    --bin-size 100 \
    --dsname NA19240 \
    --bslabel "RRBS" \
    --regionLabel "CTCF"

## Plot CTCF for W562
${Rscript} $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/K562.bin200.flank2000.CTCF.outMatrix.tsv.gz \
    --bin-size 200 \
    --dsname K562 \
    --bslabel "WGBS" \
    --regionLabel "CTCF"

## Plot CTCF for HL60
${Rscript} $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/HL60.bin200.flank2000.CTCF.outMatrix.tsv.gz \
    --bin-size 200 \
    --dsname HL60 \
    --bslabel "RRBS" \
    --regionLabel "CTCF"

####################################
####################################
### TSS plot for NA19240
${Rscript} $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/NA19240.bin50.TSS.outMatrix.tsv.gz \
    --dsname NA19240 \
    --bslabel "RRBS"

### TSS plot for APL
${Rscript} $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/APL.bin50.TSS.outMatrix.tsv.gz \
    --dsname APL \
    --bslabel "WGBS"

### TSS plot for K562
${Rscript} $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/K562.bin200.TSS.outMatrix.tsv.gz \
    --bin-size 200 \
    --dsname K562 \
    --bslabel "WGBS"

### TSS plot for HL60
${Rscript} $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/HL60.bin200.TSS.outMatrix.tsv.gz \
    --bin-size 200 \
    --dsname HL60 \
    --bslabel "RRBS"
