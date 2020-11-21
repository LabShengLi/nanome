#!/bin/bash
#SBATCH --job-name=methcall.DeepMod
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 5-14:00 # time (D-HH:MM)
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

bash
#source /projects/liuya/workspace/tcgajax/nanocompare/methcall/conda_setup.sh

# default DeepMod model location
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
deepModModel="/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/DeepMod/DeepMod/train_mod/rnn_conmodC_P100wd21_f7ne1u0_4/mod_train_conmodC_P100wd21_f3ne1u0"

processors=8

jobindex=$((SLURM_ARRAY_TASK_ID-1))

set -u
basecalledDir=$1
methCalledDir=$2
set +u

## Modify directory for processed files after basecalling:
basecallInput=${basecalledDir}/${jobindex}/workspace/pass/0

# output
methCalledDir=${methCalledDir}/DeepMod/${jobindex}

mkdir -p ${methCalledDir}

conda activate nanoai

## Call methylation from processed fast5 files:
python /projects/liuya/tools/DeepMod/bin/DeepMod.py detect \
	--wrkBase ${basecallInput} --Ref $refGenome --FileID batch_${jobindex} --modfile $deepModModel \
	--threads $processors --outFolder ${methCalledDir}

conda deactivate

echo "###   DeepMod methylcation calling DONE    ###"
