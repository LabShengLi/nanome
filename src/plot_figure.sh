#!/bin/bash
#SBATCH --job-name=plot-figure
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem 250G # memory pool for all cores
#SBATCH -t 1-20:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

# revert master by ly1


set -x

pythonFile=plot_figure.py

# Step 1: Table S2,S3
#python plot_figure.py export-data -i /projects/li-lab/yang/results/2021-05-19/MethPerf-HL60_RRBS_2Reps /projects/li-lab/yang/results/2021-05-19/MethPerf-K562_WGBS_2Reps /projects/li-lab/yang/results/2021-05-19/MethPerf-APL_RRBS /projects/li-lab/yang/results/2021-05-19/MethPerf-NA19240_RRBS_2Reps

### Plot figure 5a
#find /projects/li-lab/Nanopore_compare/result/meth-exp/MethCorr-cut5-tool3 -name "Meth_corr_plot_data_joined-*.csv" -exec python ${pythonFile} fig5a -i {} -o . \;

# Step 4: Figure 5B data
#python ${pythonFile} export-corr-data  --beddir /projects/li-lab/yang/results/2021-04-12 \
#	-i /projects/li-lab/yang/results/2021-04-13/MethCorr-HL60_RRBS_2Reps \
#	 /projects/li-lab/yang/results/2021-04-13/MethCorr-K562_WGBS_2Reps \
#	 /projects/li-lab/yang/results/2021-04-13/MethCorr-APL_RRBS \
#	 /projects/li-lab/yang/results/2021-04-13/MethCorr-NA19240_RRBS_2Reps \

# Step 2: Export curve data for Figure 3B
#python plot_figure.py export-curve-data -i /projects/li-lab/yang/results/2021-04-12/MethPerf-HL60_RRBS_2Reps /projects/li-lab/yang/results/2021-04-12/MethPerf-K562_WGBS_2Reps /projects/li-lab/yang/results/2021-04-12/MethPerf-APL_RRBS /projects/li-lab/yang/results/2021-04-12/MethPerf-NA19240_RRBS_2Reps

# Step 3: Figure 3B Plot curves data
#find /projects/li-lab/yang/results/2021-04-13/plot-curve-data -name '*.pkl' -exec python plot_figure.py plot-curve-data -i {} \;

## Plot figure
# find /projects/li-lab/yang/results/2021-04-13/plot-curve-data -name 'MethPerf-NA19240*.pkl' -exec python plot_figure.py plot-curve-data -i {} \;