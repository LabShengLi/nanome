#!/bin/bash
#SBATCH --job-name=plot-figure
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem 250G # memory pool for all cores
#SBATCH -t 20:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

prj_dir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${prj_dir}/src/plot_figure.py

mkdir -p log


## Plot figure 5a
find /projects/li-lab/Nanopore_compare/result/meth-exp/MethCorr-cut5-tool3 -name "Meth_corr_plot_data_joined-*.csv" -exec python ${pythonFile} fig5a -i {} \;

### apply cov cutoff and convert to bedGraph format for TSS plot
#flist=$(find /projects/li-lab/yang/results/2021-04-22 -name "*.cov1.bed")
#
#for fn in $flist; do # Not recommended, will break on whitespace
#    echo process "$fn"
#    bglabel=bgtruth
#    if [[ "$fn" == *"$bglabel"* ]]; then
#        CUTOFF=5
#    else
#        CUTOFF=3
#    fi
#    python ${pythonFile} bed-to-bedGraph -i $fn --cutoff $CUTOFF
#done

# Step 4: Figure 5B data
#python ${pythonFile} export-corr-data  --beddir /projects/li-lab/yang/results/2021-04-12 \
#	-i /projects/li-lab/yang/results/2021-04-13/MethCorr-HL60_RRBS_2Reps \
#	 /projects/li-lab/yang/results/2021-04-13/MethCorr-K562_WGBS_2Reps \
#	 /projects/li-lab/yang/results/2021-04-13/MethCorr-APL_RRBS \
#	 /projects/li-lab/yang/results/2021-04-13/MethCorr-NA19240_RRBS_2Reps \

# Step 1: Table S2,S3
#python plot_figure.py export-data -i /projects/li-lab/yang/results/2021-04-12/MethPerf-HL60_RRBS_2Reps /projects/li-lab/yang/results/2021-04-12/MethPerf-K562_WGBS_2Reps /projects/li-lab/yang/results/2021-04-12/MethPerf-APL_RRBS /projects/li-lab/yang/results/2021-04-12/MethPerf-NA19240_RRBS_2Reps

# Step 2: Export curve data for Figure 3B
#python plot_figure.py export-curve-data -i /projects/li-lab/yang/results/2021-04-12/MethPerf-HL60_RRBS_2Reps /projects/li-lab/yang/results/2021-04-12/MethPerf-K562_WGBS_2Reps /projects/li-lab/yang/results/2021-04-12/MethPerf-APL_RRBS /projects/li-lab/yang/results/2021-04-12/MethPerf-NA19240_RRBS_2Reps

# Step 3: Figure 3B Plot curves data
#find /projects/li-lab/yang/results/2021-04-13/plot-curve-data -name '*.pkl' -exec python plot_figure.py plot-curve-data -i {} \;
