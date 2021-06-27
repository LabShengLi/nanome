#!/bin/bash
#SBATCH --job-name=meth.corr
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=300g
#SBATCH --time=72:00:00
#SBATCH -o log/%x.%j.out
#SBATCH -e log/%x.%j.err

dsname=${1:-HL60}
runidSuffix=${2:-""} # "_RRBS_2Reps
command=${3:-"six"} # six - six tools, METEORE - include METEORE
bedDir=${4:-"/projects/li-lab/yang/results/2021-06-26"} # "_RRBS_2Reps

if [ ! -z "$runidSuffix" ]; then
	runidSuffix="_${runidSuffix}"
fi

baseDir="/projects/li-lab/Nanopore_compare/data"

## Detect calls by tools
NanopolishCalls=$(find $baseDir/$dsname -name "$dsname.nanopolish.methylation_calls.combine.*.gz")
MegalodonCalls=$(find $baseDir/$dsname -name "$dsname.megalodon.per_read.sorted.*.gz")
DeepSignalCalls=$(find $baseDir/$dsname -name "$dsname.deepsignal.call_mods.combine.*.gz")
GuppyCalls=$(find /projects/li-lab/Nanopore_compare/suppdata/guppy_methcall_yang/ -name "$dsname.guppy.fast5mod_site_level.combine.*.gz")
TomboCalls=$(find $baseDir/$dsname -name "$dsname.tombo.*.combine.*.gz")
DeepModCalls=$(find $baseDir/$dsname -name "$dsname.deepmod.C_clusterCpG.combine.*.gz")
METEORECalls=$(find /projects/li-lab/Nanopore_compare/suppdata/METEORE_results -name "$dsname.METEORE.megalodon_deepsignal-optimized-model-perRead.combine.tsv.gz")

## Automatically get the input files for each tool
declare -A bgtruthEncodeDict
bgtruthEncodeDict=(["HL60"]="bismark" ["K562"]="encode" ["APL"]="bismark" ["NA19240"]="bismark")

declare -A bgtruthFilenameDict
bgtruthFilenameDict=(
	["HL60"]="/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDA.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDF.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz"
	["K562"]="/projects/li-lab/Nanopore_compare/data/K562/ENCFF721JMB.bed.gz;/projects/li-lab/Nanopore_compare/data/K562/ENCFF867JRG.bed.gz"
	["APL"]="/projects/li-lab/Nanopore_compare/data/APL/APL-bs_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.convert.add.strand.tsv.gz"
	["NA19240"]="/projects/li-lab/Nanopore_compare/data/NA19240/NA19240_RRBS_ENCFF000LZS_rep1.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/NA19240/NA19240_RRBS_ENCFF000LZT_rep2.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz"
)

echo "####################"
echo "encode=${bgtruthEncodeDict[$dsname]}"
echo "bgtruthFilename=${bgtruthFilenameDict[$dsname]}"

echo "####################"
echo "NanopolishCalls=$NanopolishCalls"
echo "MegalodonCalls=$MegalodonCalls"
echo "DeepSignalCalls=$DeepSignalCalls"
echo "GuppyCalls=$GuppyCalls"
echo "TomboCalls=$TomboCalls"
echo "DeepModCalls=$DeepModCalls"
echo "METEORECalls=$METEORECalls"

echo "####################"
pythonFile=nanocompare/site_level_eval.py

set -x

if [ "${command}" == "six" ] ; then
	## six tool
	python ${pythonFile} \
		--calls \
			Nanopolish:${NanopolishCalls} \
			Megalodon.ZW:${MegalodonCalls} \
			DeepSignal:${DeepSignalCalls} \
			Guppy:${GuppyCalls} \
			Tombo:${TomboCalls} \
			DeepMod.Cluster:${DeepModCalls} \
		--bgtruth ${bgtruthEncodeDict[$dsname]}:${bgtruthFilenameDict[$dsname]} \
		--runid MethCorr-${dsname}${runidSuffix} \
		--dsname ${dsname} \
		--toolcov-cutoff 3 --min-bgtruth-cov 5  \
		--beddir ${bedDir} \
		--enable-cache --using-cache --gen-venn
elif [ "${command}" == "METEORE" ] ; then
	## include METEORE
	python ${pythonFile} \
		--calls \
			Nanopolish:${NanopolishCalls} \
			Megalodon.ZW:${MegalodonCalls} \
			DeepSignal:${DeepSignalCalls} \
			Guppy:${GuppyCalls} \
			Tombo:${TomboCalls} \
			DeepMod.Cluster:${DeepModCalls} \
			METEORE:${METEORECalls} \
		--bgtruth ${bgtruthEncodeDict[$dsname]}:${bgtruthFilenameDict[$dsname]} \
		--runid MethCorr-${dsname}${runidSuffix} \
		--dsname ${dsname} \
		--toolcov-cutoff 3 --min-bgtruth-cov 5  \
		--beddir ${bedDir} \
		--enable-cache --using-cache --gen-venn
else
	echo "Not support command"
	exit 19
fi
echo "### DONE"