#!/bin/bash

set -x

/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/NA19240.bin50.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --dsname NA19240

/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/APL.bin50.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --dsname APL

/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/HL60.bin200.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --bin-size 200 \
    --dsname HL60 \
    --bslabel "RRBS"

/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/K562.bin200.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --bin-size 200 \
    --dsname K562 \
    --bslabel "RRBS"

/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/HL60.bin100.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --bin-size 100 \
    --dsname HL60 \
    --bslabel "RRBS"

/Library/Frameworks/R.framework/Resources/bin/Rscript tss_plot.R \
    -i ../result/K562.bin100.outMatrix.tsv.gz \
    -o ../figures/tss-plot \
    --bin-size 100 \
    --dsname K562 \
    --bslabel "RRBS"
