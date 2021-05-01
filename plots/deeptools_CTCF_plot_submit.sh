#!/bin/bash

set -x

#sbatch deeptools_CTCF_plot.sh NA19240 "compute plot" 50
sbatch deeptools_CTCF_plot.sh NA19240 "compute plot" 100
#sbatch deeptools_CTCF_plot.sh NA19240 "compute plot" 200

#sbatch deeptools_CTCF_plot.sh K562 "compute plot" 50
#sbatch deeptools_CTCF_plot.sh K562 "compute plot" 200
#sbatch deeptools_CTCF_plot.sh K562 "compute plot" 200
#sbatch deeptools_CTCF_plot.sh HL60 "compute plot" 50
#sbatch deeptools_CTCF_plot.sh HL60 "compute plot" 200
#sbatch deeptools_CTCF_plot.sh HL60 "compute plot" 200
