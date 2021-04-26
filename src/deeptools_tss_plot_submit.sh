#!/bin/bash

set -x

#sbatch tss_plot.sh NA19240
#sbatch tss_plot.sh APL
#sbatch tss_plot.sh HL60 "compute plot" 100
sbatch tss_plot.sh HL60 "compute plot" 150
sbatch tss_plot.sh HL60 "compute plot" 200
#sbatch tss_plot.sh K562 "compute plot" 100
sbatch tss_plot.sh K562 "compute plot" 150
sbatch tss_plot.sh K562 "compute plot" 200