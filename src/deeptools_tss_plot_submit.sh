#!/bin/bash

set -x

#sbatch deeptools_tss_plot.sh NA19240 "compute plot"
#sbatch deeptools_tss_plot.sh APL "compute plot"
sbatch deeptools_tss_plot.sh HL60 "compute plot" 200
sbatch deeptools_tss_plot.sh K562 "compute plot" 200


#sbatch deeptools_tss_plot.sh HL60 "compute plot" 100
#sbatch deeptools_tss_plot.sh K562 "compute plot" 100
