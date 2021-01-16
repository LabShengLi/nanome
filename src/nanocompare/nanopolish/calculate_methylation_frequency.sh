#!/bin/bash
#SBATCH --job-name=cal.meth.nanop
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 5-14:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.figures # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


prjBaseDir=/projects/li-lab/nmf_epihet/tcgajax
calMethFreq=${prjBaseDir}/nanocompare/nanopolish/calculate_methylation_frequency.py

methCallFileK562=/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.tsv

methFreqFileK562=K562.methylation_frequency.tsv

python ${calMethFreq} -n num_cpgs --split-groups ${methCallFileK562} > ${methFreqFileK562}
