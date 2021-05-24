#!/bin/bash

dsname=(HL60 K562 APL NA19240)
runid=(ExpPRData-HL60_RRBS_2Reps ExpPRData-K562_WGBS_2Reps ExpPRData-APL_WGBS ExpPRData-NA19240_RRBS_2Reps)
joinedSetFn=("/projects/li-lab/yang/results/2021-05-19/MethPerf-HL60_RRBS_2Reps/HL60_RRBS_2Reps.Tools_BGTruth_cov5_Joined.bed" \
			"/projects/li-lab/yang/results/2021-05-19/MethPerf-K562_WGBS_2Reps/K562_WGBS_2Reps.Tools_BGTruth_cov5_Joined.bed" \
			"/projects/li-lab/yang/results/2021-05-19/MethPerf-APL_RRBS/APL_RRBS.Tools_BGTruth_cov5_Joined.bed" \
			"/projects/li-lab/yang/results/2021-05-19/MethPerf-NA19240_RRBS_2Reps/NA19240_RRBS_2Reps.Tools_BGTruth_cov5_Joined.bed" )
bgtruthParams=("bismark:/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDA.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDF.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz" \
			"encode:/projects/li-lab/Nanopore_compare/data/K562/ENCFF721JMB.bed;/projects/li-lab/Nanopore_compare/data/K562/ENCFF867JRG.bed" \
			"bismark:/projects/li-lab/Nanopore_compare/data/APL/APL-bs_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.convert.add.strand.tsv.gz" \
			"bismark:/projects/li-lab/Nanopore_compare/data/NA19240/NA19240_RRBS_ENCFF000LZS_rep1.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/NA19240/NA19240_RRBS_ENCFF000LZT_rep2.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz" )
callsDir=("/projects/li-lab/Nanopore_compare/data/HL60" \
			"/projects/li-lab/Nanopore_compare/data/K562" \
			"/projects/li-lab/Nanopore_compare/data/APL" \
			"/projects/li-lab/Nanopore_compare/data/NA19240" )

i=0
while [ $i -ne 1 ]
do
	echo sbatch --job-name roc.pr.data."${dsname[$i]}" exp_read_level_roc_pr_curve_for_guppy.sbatch \
		"${dsname[$i]}" "${runid[$i]}" "${joinedSetFn[$i]}" "${bgtruthParams[$i]}" "${callsDir[$i]}"
	sbatch --job-name roc.pr.data."${dsname[$i]}" exp_read_level_roc_pr_curve_for_guppy.sbatch \
		"${dsname[$i]}" "${runid[$i]}" "${joinedSetFn[$i]}" "${bgtruthParams[$i]}" "${callsDir[$i]}"
	i=$(($i+1))
done