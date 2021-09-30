#!/bin/bash
#SBATCH --job-name=to.unified.format
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=300g
#SBATCH --time=72:00:00
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR

## Read level unified bed format (same as the METEORE input) for METEORE usage
## Output:  ID	Chr	Pos	Strand	Score
#           14a6f90a-220d-45dd-b699-d7b8c51a41eb	chr8	1579304	+	-1.2320743055695798

dsname=${1:-HL60}

baseDir="/projects/li-lab/Nanopore_compare/data"
guppyDataDir="/projects/li-lab/Nanopore_compare/nanome_paper_result/supp_data/guppy_methcall_yang"

if [ ${dsname} == "NA12878" ]; then
    ## Detect calls by tools for NA12878
    baseDir="/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs"
    NanopolishCalls=$(find ${baseDir} -name "$dsname.nanopolish.methylation_calls.combine*.*.gz")
    MegalodonCalls=$(find ${baseDir} -name "$dsname.megalodon.per_read.combine_allchrs.bed.gz")
    DeepSignalCalls=$(find ${baseDir} -name "$dsname.deepsignal.call_mods.combine_allchrs.tsv.gz")
    GuppyCalls=$(find ${baseDir} -name "$dsname.guppy.gcf52ref_read_level.combine_allchrs.tsv.gz")
    TomboCalls=$(find ${baseDir} -name "$dsname.tombo.perReadsStats.combine_allchrs.bed.gz")
    MegalodonEncode="Megalodon"
else
    ## Detect calls by tools for HL60, K562, APL, and NA19240
    NanopolishCalls=$(find $baseDir/$dsname -name "$dsname.nanopolish.methylation_calls.combine.*.gz")
    MegalodonCalls=$(find $baseDir/$dsname -name "$dsname.megalodon.per_read.sorted.*.gz")
    DeepSignalCalls=$(find $baseDir/$dsname -name "$dsname.deepsignal.call_mods.combine.*.gz")
    # file like NA19240.guppy.gcf52ref_read_level.combine.tsv.gz
    GuppyCalls=$(find ${guppyDataDir} -name "$dsname.guppy.gcf52ref_read_level.combine.*.gz")
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

pythonFile=tss_eval.py

set -x

## convert tools format to unified format (such as METEORE)
${pythonFile} \
    --calls \
        Nanopolish:${NanopolishCalls} \
        ${MegalodonEncode}:${MegalodonCalls} \
        DeepSignal:${DeepSignalCalls} \
        Guppy.gcf52ref:${GuppyCalls} \
        Tombo:${TomboCalls} \
    --runid TSS_METEORE-${dsname} \
    --dsname ${dsname} --output-unified-format \
    --processors 8

echo "### DONE"
