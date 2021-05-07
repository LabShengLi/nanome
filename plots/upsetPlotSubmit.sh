#!/bin/bash

set -x

export PATH=/Library/Frameworks/R.framework/Resources/bin:$PATH

baseDir=plots

Rscript $baseDir/upsetPlot.R --dsname APL

Rscript $baseDir/upsetPlot.R --dsname HL60

Rscript $baseDir/upsetPlot.R --dsname K562

Rscript $baseDir/upsetPlot.R --dsname NA19240
