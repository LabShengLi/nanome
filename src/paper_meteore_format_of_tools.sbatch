#!/bin/bash
#SBATCH --job-name=to.meteore.format
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=300g
#SBATCH --time=72:00:00
#SBATCH -o log/%x.%j.out
#SBATCH -e log/%x.%j.err

dsname=${1:-HL60}

baseDir="/projects/li-lab/Nanopore_compare/data"

if [ ${dsname} == "NA12878" ]; then
	## Detect calls by tools for NA12878
	NanopolishCalls=$(find /projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs -name "$dsname.nanopolish.methylation_calls.combine*.*.gz")
	MegalodonCalls=$(find /projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs -name "$dsname.megalodon.per_read.combine_allchrs.bed.gz")
	DeepSignalCalls=$(find /projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs -name "$dsname.deepsignal.call_mods.combine_allchrs.tsv.gz")
	GuppyCalls=$(find /projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs -name "$dsname.guppy.gcf52ref_read_level.combine_allchrs.tsv.gz")
	TomboCalls=$(find /projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs -name "$dsname.tombo.perReadsStats.combine_allchrs.bed.gz")
	MegalodonEncode="Megalodon"
else
	## Detect calls by tools for HL60, K562, APL, and NA19240
	NanopolishCalls=$(find $baseDir/$dsname -name "$dsname.nanopolish.methylation_calls.combine.*.gz")
	MegalodonCalls=$(find $baseDir/$dsname -name "$dsname.megalodon.per_read.sorted.*.gz")
	DeepSignalCalls=$(find $baseDir/$dsname -name "$dsname.deepsignal.call_mods.combine.*.gz")
	GuppyCalls=$(find /projects/li-lab/Nanopore_compare/suppdata/guppy_methcall_yang -name "$dsname.guppy.fast5mod_site_level.combine.*.gz")
	TomboCalls=$(find $baseDir/$dsname -name "$dsname.tombo.*.combine.*.gz")
	MegalodonEncode="Megalodon.ZW"

fi

echo "####################"
echo "NanopolishCalls=$NanopolishCalls"
echo "MegalodonCalls=$MegalodonCalls"
echo "DeepSignalCalls=$DeepSignalCalls"
echo "GuppyCalls=$GuppyCalls"
echo "TomboCalls=$TomboCalls"
echo "####################"

pythonFile=nanocompare/tss_eval.py

set -x

## convert tools format to format of METEORE
python ${pythonFile} \
	--calls \
	Nanopolish:${NanopolishCalls} \
	${MegalodonEncode}:${MegalodonCalls} \
	DeepSignal:${DeepSignalCalls} \
	Guppy.gcf52ref:${GuppyCalls} \
	Tombo:${TomboCalls} \
	--runid METEORE_FORMAT-${dsname} \
	--dsname ${dsname} \
	--output-meteore

echo "### DONE"
