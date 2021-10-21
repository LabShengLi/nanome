#!/bin/bash
#SBATCH --job-name=to.meteore.format.ecoli
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=300g
#SBATCH --time=72:00:00
#SBATCH -o log/%x.%j.out
#SBATCH -e log/%x.%j.err

baseDir="/projects/li-lab/Nanopore_compare/suppdata/ecoli-sanity-check/ECOLI_METROPAPER-methylation-callings"

## Detect calls by tools for Ecoli
NanopolishCalls=$(find $baseDir -name "ECOLI_METROPAPER.Nanopolish.combine.tsv.gz")
MegalodonCalls=$(find $baseDir -name "ECOLI_METROPAPER.Megalodon.combine.tsv.gz")
DeepSignalCalls=$(find $baseDir -name "ECOLI_METROPAPER.DeepSignal.combine.tsv.gz")
GuppyCalls=$(find $baseDir -name "EcoliDemo.guppy.gcf52ref_read_level.combine.tsv.gz")
TomboCalls=$(find $baseDir -name "ECOLI_METROPAPER.Tombo.combine.tsv.gz")
MegalodonEncode="Megalodon"

echo "####################"
echo "NanopolishCalls=$NanopolishCalls"
echo "MegalodonCalls=$MegalodonCalls"
echo "DeepSignalCalls=$DeepSignalCalls"
echo "GuppyCalls=$GuppyCalls"
echo "TomboCalls=$TomboCalls"
echo "MegalodonEncode=$MegalodonEncode"
echo "####################"

pythonFile=${NanoCompareDir}/src/nanocompare/tss_eval.py

set -x

## convert tools format to format of METEORE
python ${pythonFile} \
    --calls \
    Nanopolish:${NanopolishCalls} \
    ${MegalodonEncode}:${MegalodonCalls} \
    DeepSignal:${DeepSignalCalls} \
    Guppy.gcf52ref:${GuppyCalls} \
    Tombo:${TomboCalls} \
    --runid METEORE_FORMAT-Ecoli \
    --dsname Ecoli \
    --chrs NC_000913.3 \
    --output-meteore

echo "### DONE"
