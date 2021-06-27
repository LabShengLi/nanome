#!/bin/bash
#SBATCH --job-name=plot-figure
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300G # memory pool for all cores
#SBATCH -t 72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

## Script used for generate data and plot figures
## Parameters:
## 				command
## 				outputDir

#set -x

pythonFile=plot_figure.py

command=${1:-"Step1"}
outputDir=${2:-"/projects/li-lab/yang/results/2021-06-27"}

bedDir="/projects/li-lab/yang/results/2021-06-26"

if [ $command == "Step1" ]; then
	# Step 1: Table S2,S3
	python plot_figure.py export-data -i \
		${outputDir}/MethPerf-HL60_RRBS_2Reps \
		${outputDir}/MethPerf-K562_WGBS_2Reps \
		${outputDir}/MethPerf-APL_WGBS \
		${outputDir}/MethPerf-NA19240_RRBS_2Reps

elif [ $command == "Step1_METEORE" ]; then
	# bash plot_figure.sh Step1_METEORE
	python plot_figure.py export-data \
		-i 	${outputDir}/MethPerf-HL60_RRBS_2Reps_METEORE \
			${outputDir}/MethPerf-K562_WGBS_2Reps_METEORE \
			${outputDir}/MethPerf-APL_WGBS_METEORE \
			${outputDir}/MethPerf-NA19240_RRBS_2Reps_METEORE \
		--tagname METEORE

elif [ $command == "Step2" ]; then
	# Step 2.1: Export curve data for Figure 3B
	# sbatch plot_figure.sh Step2
	python plot_figure.py export-curve-data -i \
	 		${outputDir}/MethPerf-HL60_RRBS_2Reps \
	 		${outputDir}/MethPerf-K562_WGBS_2Reps \
	 		${outputDir}/MethPerf-APL_WGBS \
	 		${outputDir}/MethPerf-NA19240_RRBS_2Reps

	# Step 2.2: Figure 3B Plot curves data
	find /projects/li-lab/yang/results/$(date +%F)/plot-curve-data -name '*.pkl' -exec python plot_figure.py plot-curve-data -i {} \;

elif [ $command == "Fig5a" ]; then
	### Plot figure 5a
	# sbatch plot_figure.sh Fig5a
	fileList=$(find ${outputDir} -name "Meth_corr_plot_data_joined-*.csv")
	for fn in $fileList; do
	 	python ${pythonFile} fig5a -i $fn &
	done
	wait

elif [ $command == "Fig5b-data" ]; then
	### Figure 5B data: COE on each region data for bar plot
	# sbatch plot_figure.sh Fig5b-data
	python ${pythonFile} export-corr-data  --beddir ${bedDir} \
		-i	${outputDir}/MethCorr-HL60_RRBS_2Reps \
			${outputDir}/MethCorr-K562_WGBS_2Reps \
			${outputDir}/MethCorr-APL_WGBS \
			${outputDir}/MethCorr-NA19240_RRBS_2Reps

elif [ $command == "Fig5b-data-METEORE" ]; then
	### Figure 5B data: COE on each region data for bar plot
	# sbatch plot_figure.sh Fig5b-data-METEORE
	python ${pythonFile} export-corr-data  --beddir ${bedDir} \
		-i 	${outputDir}/MethCorr-HL60_RRBS_2Reps_METEORE \
			${outputDir}/MethCorr-K562_WGBS_2Reps_METEORE \
			${outputDir}/MethCorr-APL_WGBS_METEORE \
			${outputDir}/MethCorr-NA19240_RRBS_2Reps_METEORE \
		--tagname METEORE
else
	echo "### Unsupported command=${command}"
	exit 3
fi

echo "### DONE"
