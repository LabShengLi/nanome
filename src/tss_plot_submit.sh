#!/bin/bash

set -x

###################################
###################################
## CTCF plot for call R script



/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/tss.CTCF.data/NA19240.bin100.flank2000.CTCF.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --bin-size 100 \
    --dsname NA19240 \
    --bslabel "RRBS" \
    --regionLabel "CTCF"

/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/tss.CTCF.data/K562.bin200.flank2000.CTCF.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --bin-size 200 \
    --dsname K562 \
    --bslabel "WGBS" \
    --regionLabel "CTCF"

exit 0

#/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
#    -i ../result/tss.CTCF.data/HL60.bin100.flank2000.CTCF.outMatrix.tsv.gz \
#    -o ../figures/tss-plot \
#    --bin-size 100 \
#    --dsname HL60 \
#    --bslabel "RRBS" \
#    --regionLabel "CTCF"

###################################
###################################
## TSS plot for call R script
/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/tss.CTCF.data/NA19240.bin50.TSS.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --dsname NA19240 \
    --bslabel "RRBS"
#
/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/tss.CTCF.data/APL.bin50.TSS.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --dsname APL \
    --bslabel "WGBS"

/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/tss.CTCF.data/K562.bin200.TSS.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --bin-size 200 \
    --dsname K562 \
    --bslabel "WGBS"

#
#/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
#    -i ../result/tss.CTCF.data/HL60.bin100.TSS.outMatrix.tsv.gz \
#    -o ../figures/tss-plot \
#    --bin-size 100 \
#    --dsname HL60 \
#    --bslabel "RRBS"
#