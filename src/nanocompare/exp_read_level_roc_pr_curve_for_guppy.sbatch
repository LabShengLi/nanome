#!/bin/bash
#SBATCH --job-name=roc.pr.guppy
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=180g
#SBATCH --time=72:00:00
#SBATCH -o log/%x.%j.out
#SBATCH -e log/%x.%j.err

set -x

pythonFn=exp_read_level_roc_pr_curve_for_guppy.py

dsname=${1:-HL60}
runid=${2:-"ExpPRData-HL60_RRBS_2Reps"}
joinedSetFn=${3:-"/projects/li-lab/yang/results/2021-05-19/MethPerf-HL60_RRBS_2Reps/HL60_RRBS_2Reps.Tools_BGTruth_cov5_Joined.bed"}
bgtruthParams=${4:-"bismark:/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDA.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDF.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz"}
callsDir=${5:-"/projects/li-lab/Nanopore_compare/data/HL60/"}

python ${pythonFn} \
	--calls-dir ${callsDir} \
	--joined-set-fn ${joinedSetFn} \
	--bgtruth ${bgtruthParams} \
	--runid ${runid} --dsname ${dsname} \
	--min-bgtruth-cov 5 --enable-cache --using-cache
