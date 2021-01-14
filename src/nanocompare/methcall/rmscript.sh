#!/bin/bash
#SBATCH --job-name=rm
##SBATCH -q batch
##SBATCH --partition=compute
#SBATCH -p gpu
#SBATCH -q inference
#SBATCH --gres=gpu:1
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 50g # memory pool for all cores
#SBATCH -t 02:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

# Run script on Sumner: sbatch --partition=compute -q batch --gres=  rmscript.sh

#set -x

#cd /projects/li-lab/ctimage/yang/results-ctimage-early-phase
#
##rm -rf 08-* 09-* 10-* 11-*
##rm -rf 07-2*
##rm -rf 07-15 07-16 07-18 07-19
#
#
#rm -rf 07-22  09-22  10-28

source /home/liuya/.bash_profile


conda activate nanoai

conda env list

conda deactivate

conda env list