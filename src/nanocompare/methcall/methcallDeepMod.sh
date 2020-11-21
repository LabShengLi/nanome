#!/bin/bash
#SBATCH --job-name=methcall.DeepMod
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 5-14:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/liuya/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/liuya/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/liuya/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/liuya/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


set -x

echo JOBINDEX=$((SLURM_ARRAY_TASK_ID-1))

jobindex=$((SLURM_ARRAY_TASK_ID-1))

targetNum=3

processors=8

# input
basecalledFast5DIR=/projects/liuya/results/nanocompare/demo/basecalled

basecalledFast5Input=${basecalledFast5DIR}/${jobindex}

# output
methCallsDir=/projects/liuya/results/nanocompare/demo/methcalledDeepMod


# default DeepMod model location
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
deepModModel="/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/DeepMod/DeepMod/train_mod/rnn_conmodC_P100wd21_f7ne1u0_4/mod_train_conmodC_P100wd21_f3ne1u0"

mkdir -p ${methCallsDir}

set +x
conda activate nanoai

conda info

set -x

## Modify directory for processed files after basecalling:
basecallInput=${basecalledFast5DIR}/${jobindex}/workspace/pass/0

## Call methylation from processed fast5 files:
python /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/DeepMod/DeepMod/bin/DeepMod.py detect --wrkBase ${basecallInput} --Ref $refGenome --FileID batch_${jobindex} --modfile $deepModModel --threads $processors --outFolder $methCallsDir


set +x
conda deactivate
set -x



