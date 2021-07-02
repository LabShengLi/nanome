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
## 				resultsDir

pythonFile=plot_figure.py

command=${1:-"Step1"}
type=${2:-"METEORE"}
resultsDir=${3:-"/projects/li-lab/yang/results/2021-06-29"}

bedDir="/projects/li-lab/yang/results/2021-06-26"

if [ "$type" == "six" ]; then
	HL60_Result_dir=${resultsDir}/MethPerf-HL60_RRBS_2Reps
	K562_Result_dir=${resultsDir}/MethPerf-K562_WGBS_2Reps
	APL_Result_dir=${resultsDir}/MethPerf-APL_WGBS
	NA19240_Result_dir=${resultsDir}/MethPerf-NA19240_RRBS_2Reps
	tagOptions=""
	plotCurveDataDir="plot-curve-data"
elif [ "$type" == "METEORE" ]; then
	HL60_Result_dir=${resultsDir}/MethPerf-HL60_RRBS_2Reps_METEORE
	K562_Result_dir=${resultsDir}/MethPerf-K562_WGBS_2Reps_METEORE
	APL_Result_dir=${resultsDir}/MethPerf-APL_WGBS_METEORE
	NA19240_Result_dir=${resultsDir}/MethPerf-NA19240_RRBS_2Reps_METEORE
	tagOptions="--tagname METEORE"
	plotCurveDataDir="plot-curve-data-METEORE"
else
	echo "### Unsupported type=$type"
	exit 156
fi

if [ $command == "Step1" ]; then
	# Step 1: Table S2,S3
	# bash plot_figure.sh Step1
	# bash plot_figure.sh Step1 METEORE
	# bash plot_figure.sh Step1 METEORE /projects/li-lab/yang/results/2021-06-28
	python plot_figure.py export-data -i \
		${HL60_Result_dir} \
		${K562_Result_dir} \
		${APL_Result_dir} \
		${NA19240_Result_dir} ${tagOptions}

elif [ $command == "Step2" ]; then
	# Step 2.1: Export curve data for Figure 3B
	# sbatch plot_figure.sh Step2 METEORE
	python  plot_figure.py export-curve-data -i \
		${HL60_Result_dir} \
		${K562_Result_dir} \
		${APL_Result_dir} \
		${NA19240_Result_dir}    ${tagOptions}

	# Step 2.2: Figure 3B Plot curves data
	## python plot_figure.py plot-curve-data -i  /projects/li-lab/yang/results/2021-06-29/plot-curve-data-METEORE/MethPerf-NA19240_RRBS_2Reps_METEORE.plot.curve.data.ytrue.ypred.yscore.Singletons.pkl
	find  /projects/li-lab/yang/results/$(date +%F)/${plotCurveDataDir} -name '*.pkl' \
		-exec  python plot_figure.py plot-curve-data -i {} ${tagOptions} \;

#	find /projects/li-lab/yang/results/2021-06-29/${plotCurveDataDir} -name '*.pkl' \
#		-exec python plot_figure.py plot-curve-data -i {} ${tagOptions} \;

elif [ $command == "Fig5a" ]; then
	### Plot figure 5a
	# sbatch plot_figure.sh Fig5a
	## python plot_figure.py fig5a -i /projects/li-lab/yang/results/2021-06-28/MethCorr-HL60_RRBS_2Reps_METEORE/Meth_corr_plot_data_joined-HL60_RRBS_2Reps_METEORE-bsCov5-minToolCov3-baseFormat1.csv
	fileList=$(find ${resultsDir}/*_METEORE -name "Meth_corr_plot_data_joined-*.csv")
	for fn in $fileList; do
		python  ${pythonFile} fig5a -i $fn &
	done
	wait

elif [ $command == "Fig5b-data" ]; then
	### Figure 5B data: COE on each region data for bar plot
	# sbatch plot_figure.sh Fig5b-data six/METEORE
	python ${pythonFile} export-corr-data  --beddir ${bedDir} \
		-i ${HL60_Result_dir/MethPerf/MethCorr} \
		${K562_Result_dir/MethPerf/MethCorr} \
		${APL_Result_dir/MethPerf/MethCorr} \
		${NA19240_Result_dir/MethPerf/MethCorr}   ${tagOptions}
else
	echo "### Unsupported command=${command}"
	exit 3
fi

echo "### DONE"
