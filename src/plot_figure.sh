#!/bin/bash
#SBATCH --job-name=plot-figure
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 20:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


set -x

prj_dir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${prj_dir}/src/plot_figure.py

mkdir -p log

#python ${pythonFile} $@

python ${pythonFile} export-corr-data  --beddir /projects/li-lab/yang/results/2021-04-07/MethPerf-cut5 \
	-i /projects/li-lab/yang/results/2021-04-08/MethCorr-HL60_RRBS_2Reps \
	 /projects/li-lab/yang/results/2021-04-08/MethCorr-K562_WGBS_2Reps \
	 /projects/li-lab/yang/results/2021-04-08/MethCorr-APL_RRBS \
	 /projects/li-lab/yang/results/2021-04-08/MethCorr-NA19240_RRBS_2Reps \


#python plot_figure.py export-curve-data -i /projects/li-lab/yang/results/2021-04-07/MethPerf-cut5/MethPerf-HL60_RRBS_2Reps /projects/li-lab/yang/results/2021-04-07/MethPerf-cut5/MethPerf-K562_WGBS_2Reps /projects/li-lab/yang/results/2021-04-07/MethPerf-cut5/MethPerf-APL_RRBS /projects/li-lab/yang/results/2021-04-07/MethPerf-cut5/MethPerf-NA19240_RRBS_2Reps

find /projects/li-lab/yang/results/2021-04-11/plot-curve-data -name '*.pkl' -exec python plot_figure.py plot-curve-data -i {} \;
