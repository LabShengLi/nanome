#!/bin/bash
#SBATCH --job-name=meth.perf.paper
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=300g
#SBATCH --time=72:00:00
#SBATCH --output=log/%x.%j.log # STDOUT & STDERR
#SBATCH --mail-user=yang.liu@jax.org
#SBATCH --mail-type=END

dsname=${1:-HL60}      # dsname HL60, K562, APL, NA19240, NA12878
runidSuffix=${2:-"BS"} # "RRBS_2Reps"
command=${3:-"seven"}  # six - six tools, seven - all seven tools
otherOptions=${4:-""}  # if using mpi to run tasks for evaluation

genome_annotation_dir="${NanoCompareDir}/data/genome-annotation"
guppyDataDir="/projects/li-lab/Nanopore_compare/nanome_paper_result/supp_data/guppy_methcall_yang"
METEOREDataDir="/projects/li-lab/Nanopore_compare/nanome_paper_result/supp_data/METEORE_results"
baseDir="/projects/li-lab/Nanopore_compare/data"

if [ ! -z "$runidSuffix" ]; then
    runidSuffix="_${runidSuffix}"
fi

set -x
## Detect calls by tools
NanopolishCalls=$(find $baseDir/$dsname -name "$dsname.nanopolish.methylation_calls.combine.*.gz")
MegalodonCalls=$(find $baseDir/$dsname -name "$dsname.megalodon.per_read.sorted.*.gz")
DeepSignalCalls=$(find $baseDir/$dsname -name "$dsname.deepsignal.call_mods.combine.*.gz")
GuppyCalls=$(find ${guppyDataDir} -name "$dsname.guppy.fast5mod_site_level.combine.*.gz")
TomboCalls=$(find $baseDir/$dsname -name "$dsname.tombo.*.combine.*.gz")
DeepModCalls=$(find $baseDir/$dsname -name "$dsname.deepmod.C.combine.*.gz")
METEORECalls=$(find ${METEOREDataDir} -name "$dsname.METEORE.megalodon_deepsignal-optimized-model-perRead.combine.tsv.gz")
MegalodonEncode="Megalodon.ZW"

if [ "$dsname" == "NA12878" ]; then
    ## Detect tool results for NA12878 at specific folders
    baseDir="/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs"
    NanopolishCalls=$(find $baseDir -name "$dsname.nanopolish.methylation_calls.combine_allchrs.*.gz")
    MegalodonCalls=$(find $baseDir -name "$dsname.megalodon.per_read.combine_allchrs.*.gz")
    DeepSignalCalls=$(find $baseDir -name "$dsname.deepsignal.call_mods.combine_allchrs.*.gz")
    GuppyCalls=$(find $baseDir -name "$dsname.guppy.fast5mod_site_level.combine_allchrs.*.gz")
    TomboCalls=$(find $baseDir -name "$dsname.tombo.perReadsStats.combine_allchrs.*.gz")
    DeepModCalls=$(find $baseDir -name "$dsname.deepmod.C.combine_allchrs.*.gz")
    METEORECalls=$(find ${METEOREDataDir} -name "$dsname.METEORE.megalodon_deepsignal-optimized-model-perRead.combine.tsv.gz")
    MegalodonEncode="Megalodon"
fi

## Automatically get the input files for each tool
declare -A bgtruthEncodeDict
bgtruthEncodeDict=(["HL60"]="bismark" ["K562"]="encode" ["APL"]="bismark" ["NA19240"]="bismark" ["NA12878"]="encode")

## BS seq data for each dataset
declare -A bgtruthFilenameDict
bgtruthFilenameDict=(
    ["HL60"]="/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDA.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDF.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz"
    ["K562"]="/projects/li-lab/Nanopore_compare/data/K562/ENCFF721JMB.bed.gz;/projects/li-lab/Nanopore_compare/data/K562/ENCFF867JRG.bed.gz"
    ["APL"]="/projects/li-lab/Nanopore_compare/data/APL/APL-bs_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.convert.add.strand.tsv.gz"
    ["NA19240"]="/projects/li-lab/Nanopore_compare/data/NA19240/NA19240_RRBS_ENCFF000LZS_rep1.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/NA19240/NA19240_RRBS_ENCFF000LZT_rep2.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz"
    ["NA12878"]="/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz;/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz"
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
echo "METEORECalls=$METEORECalls"
echo "DeepModCalls=$DeepModCalls"
echo "####################"

## Read level evaluation script
pythonFile=read_level_eval.py

set -x

if [ "${command}" == "six" ]; then
    ## six tool
    ${pythonFile} \
        --calls \
            Nanopolish:${NanopolishCalls} \
            ${MegalodonEncode}:${MegalodonCalls} \
            DeepSignal:${DeepSignalCalls} \
            Guppy:${GuppyCalls} \
            Tombo:${TomboCalls} \
            DeepMod.C:${DeepModCalls} \
        --bgtruth ${bgtruthEncodeDict[$dsname]}:${bgtruthFilenameDict[$dsname]} \
        --runid MethPerf-${dsname}${runidSuffix} \
        --dsname ${dsname} \
        --min-bgtruth-cov 5 \
        --distribution \
        --report-joined \
        --processors 30 \
        --enable-cache --using-cache\
        --genome-annotation ${genome_annotation_dir} ${otherOptions} ### --mpi
elif [ "${command}" == "seven" ]; then
    ## --distribution: generate sites distribution in singletons/non-singletons for each genomic regions
    ${pythonFile} \
        --calls \
            Nanopolish:${NanopolishCalls} \
            ${MegalodonEncode}:${MegalodonCalls} \
            DeepSignal:${DeepSignalCalls} \
            Guppy:${GuppyCalls} \
            Tombo:${TomboCalls} \
            METEORE:${METEORECalls} \
            DeepMod.C:${DeepModCalls} \
        --bgtruth ${bgtruthEncodeDict[$dsname]}:${bgtruthFilenameDict[$dsname]} \
        --runid MethPerf-${dsname}${runidSuffix} \
        --dsname ${dsname} \
        --min-bgtruth-cov 5 \
        --report-joined \
        --enable-cache --using-cache \
        --distribution --processors 30 \
        --genome-annotation ${genome_annotation_dir} ${otherOptions} ### --mpi
else
    echo "Not support command"
    exit 255
fi
echo "### DONE"
