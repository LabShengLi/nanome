#!/bin/bash
#SBATCH --job-name=deepmod
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 2-14:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

## optional params:
processors=8
part=0
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
deepModModel="/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/DeepMod/DeepMod/train_mod/rnn_conmodC_P100wd21_f7ne1u0_4/mod_train_conmodC_P100wd21_f3ne1u0"


basecallOutputDir=/fastscratch/liuya/nanocompare/NA19240_basecalled

methCallsDir=/fastscratch/liuya/nanocompare/NA19240_methcalled
mkdir -p ${methCallsDir}

processedFast5DIR=${basecallOutputDir}/$part

## Modify directory for processed files after basecalling:
processedFast5DIR=${processedFast5DIR}/workspace/pass/0

source /projects/liuya/workspace/tcgajax/nanocompare/methcall/conda_setup.sh

conda activate nanoai

## Call methylation from processed fast5 files:
python /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/DeepMod/DeepMod/bin/DeepMod.py detect --wrkBase $processedFast5DIR --Ref $refGenome --FileID batch_$part --modfile $deepModModel --threads $processors --outFolder $methCallsDir


conda deactivate
