#!/bin/bash
#SBATCH --job-name=albcore.debug
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem=170g # memory pool for all cores
#SBATCH --time=10:50:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

##SBATCH --array=297-300      # job array index
##SBATCH -o log/%x.batch%a.%j.out # STDOUT
##SBATCH -e log/%x.batch%a.%j.err # STDERR

##SBATCH -p gpu
##SBATCH --gres=gpu:1
##SBATCH -q inference


source /home/liuya/.bash_profile

conda activate nanoai

mkdir -p /fastscratch/liuya/nanocompare/HL60-Runs/HL60-N50-basecall-Albacore-1/1

read_fast5_basecaller.py -o fastq,fast5 -t 16 -s /fastscratch/liuya/nanocompare/HL60-Runs/HL60-N50-basecall-Albacore-1/1 -i /fastscratch/liuya/nanocompare/HL60-Runs/HL60-N50-sept/1 -c r94_450bps_linear.cfg -n 20000000000 --debug

conda deactivate