#!/bin/bash
#SBATCH --job-name=NA19240.tombo.call
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem=170g # memory pool for all cores
#SBATCH --time=1-13:50:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

##SBATCH -p gpu
##SBATCH --gres=gpu:1
##SBATCH -q inference

source /home/liuya/.bash_profile

conda env list

conda activate nanoai

conda env list


time tombo detect_modifications alternative_model --fast5-basedirs /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-N300-resquiggle/43/workspace --dna --standard-log-likelihood-ratio --statistics-file-basename /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-Tombo-N300/NA19240-Tombo-N300-methcall/NA19240-Tombo-N300.batch_43 --per-read-statistics-basename /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-Tombo-N300/NA19240-Tombo-N300-methcall/NA19240-Tombo-N300.batch_43 --alternate-bases CpG --processes 16 --corrected-group RawGenomeCorrected_000

time tombo detect_modifications alternative_model --fast5-basedirs /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-N300-resquiggle/174/workspace --dna --standard-log-likelihood-ratio --statistics-file-basename /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-Tombo-N300/NA19240-Tombo-N300-methcall/NA19240-Tombo-N300.batch_174 --per-read-statistics-basename /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-Tombo-N300/NA19240-Tombo-N300-methcall/NA19240-Tombo-N300.batch_174 --alternate-bases CpG --processes 16 --corrected-group RawGenomeCorrected_000

tombo detect_modifications alternative_model --fast5-basedirs /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-N300-resquiggle/264/workspace --dna --standard-log-likelihood-ratio --statistics-file-basename /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-Tombo-N300/NA19240-Tombo-N300-methcall/NA19240-Tombo-N300.batch_264 --per-read-statistics-basename /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-Tombo-N300/NA19240-Tombo-N300-methcall/NA19240-Tombo-N300.batch_264 --alternate-bases CpG --processes 16 --corrected-group RawGenomeCorrected_000

#cd /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-DeepSignal-N300/NA19240-DeepSignal-N300-methcall
#
#> NA19240.deepsignal.call_mods.combine.tsv
#
#cat *.deepsignal.call_mods.tsv > NA19240.deepsignal.call_mods.combine.tsv
#
#wc -l NA19240.deepsignal.call_mods.combine.tsv

#deepsignal call_mods --input_path /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-N300-resquiggle/47/workspace --model_path /projects/li-lab/yang/workspace/nano-compare/data/dl-model/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt --result_file /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-DeepSignal-N300/NA19240-DeepSignal-N300-methcall/NA19240-DeepSignal-N300.batch_47.CpG.deepsignal.call_mods.tsv --reference_path /projects/li-lab/reference/hg38/hg38.fasta --corrected_group RawGenomeCorrected_000 --nproc 16 --is_gpu no

conda deactivate