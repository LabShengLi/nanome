#!/bin/bash
#SBATCH --job-name=NA19240.failed
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem=270g # memory pool for all cores
#SBATCH --time=2-03:50:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

##SBATCH -p gpu
##SBATCH --gres=gpu:1
##SBATCH -q inference

source /home/liuya/.bash_profile

conda env list

conda activate nanoai

conda env list

deepsignal call_mods --input_path /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-N300-resquiggle/47/workspace --model_path /projects/li-lab/yang/workspace/nano-compare/data/dl-model/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt --result_file /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-DeepSignal-N300/NA19240-DeepSignal-N300-methcall/NA19240-DeepSignal-N300.batch_47.CpG.deepsignal.call_mods.tsv --reference_path /projects/li-lab/reference/hg38/hg38.fasta --corrected_group RawGenomeCorrected_000 --nproc 16 --is_gpu no

conda deactivate