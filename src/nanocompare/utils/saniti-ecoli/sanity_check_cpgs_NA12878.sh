#!/bin/bash
#SBATCH --job-name=santi.cpgs
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300G # memory pool for all cores
#SBATCH -t 72:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

cpgs=(159199943 159200094 159200216 159200386)

infn_meteore='/projects/li-lab/Nanopore_compare/suppdata/METEORE_results/NA19240.METEORE.megalodon_deepsignal-optimized-model-perRead.combine.tsv.gz'
infn_megalodon='/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.megalodon.per_read.combine_allchrs.bed.gz'
infn_deepsignal='/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.deepsignal.call_mods.combine_allchrs.tsv.gz'

set -x

zcat $infn_meteore | egrep '159199943' >>meteore_out1.txt &
zcat $infn_megalodon | egrep '159199942' >>megalodon_out1.txt &
zcat $infn_deepsignal | egrep '159199942' >>deepsignal_out1.txt &

#for cpg in "${cpgs[@]}"; do
#    echo $cpg
#    #    zcat $infn_meteore | grep $cpg >>meteore_out.txt &
#    #    zcat $infn_megalodon | grep $cpg >>megalodon_out.txt &
#    #    zcat $infn_deepsignal | grep $cpg >>deepsignal_out.txt &
#
#done

wait
echo "DONE"
