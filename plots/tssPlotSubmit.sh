#!/bin/bash

set -x
set -e

export PATH=/Library/Frameworks/R.framework/Resources/bin:$PATH

baseDir=plots
inputDir=result


###################################
###################################
## CTCF plot for call R script

Rscript $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/NA19240.bin100.flank2000.CTCF.outMatrix.tsv.gz \
    --bin-size 100 \
    --dsname NA19240 \
    --bslabel "RRBS" \
    --regionLabel "CTCF"

Rscript $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/K562.bin200.flank2000.CTCF.outMatrix.tsv.gz \
    --bin-size 200 \
    --dsname K562 \
    --bslabel "WGBS" \
    --regionLabel "CTCF"


Rscript $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/HL60.bin200.flank2000.CTCF.outMatrix.tsv.gz \
    --bin-size 200 \
    --dsname HL60 \
    --bslabel "RRBS" \
    --regionLabel "CTCF"

####################################
####################################
### TSS plot for call R script
Rscript $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/NA19240.bin50.TSS.outMatrix.tsv.gz \
    --dsname NA19240 \
    --bslabel "RRBS"
#
Rscript $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/APL.bin50.TSS.outMatrix.tsv.gz \
    --dsname APL \
    --bslabel "WGBS"
#
Rscript $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/K562.bin200.TSS.outMatrix.tsv.gz \
    --bin-size 200 \
    --dsname K562 \
    --bslabel "WGBS"



Rscript $baseDir/tssPlot.R \
    -i $inputDir/tss.CTCF.data/HL60.bin200.TSS.outMatrix.tsv.gz \
    --bin-size 200 \
    --dsname HL60 \
    --bslabel "RRBS"

