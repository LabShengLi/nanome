#!/bin/bash
#SBATCH --job-name=cp
##SBATCH -q batch
#SBATCH --partition=compute
##SBATCH -p gpu
##SBATCH -q inference
##SBATCH --gres=gpu:1
#SBATCH -N 1 # number of nodes
#SBATCH -n 3 # number of cores
#SBATCH --mem 50g # memory pool for all cores
#SBATCH -t 2-12:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

# Run script on Sumner: sbatch --partition=compute -q batch --gres=  rmscript.sh

set -x
#mv /fastscratch/liuya/nanocompare/HL60-Runs/HL60-N50-sept /projects/li-lab/Nanopore_compare/nanopore_fast5
cp -rf /fastscratch/liuya/nanocompare/APL-Runs/APL-N50-sept /projects/li-lab/Nanopore_compare/nanopore_fast5
cp -rf /fastscratch/liuya/nanocompare/K562-Runs/K562-N50-sept /projects/li-lab/Nanopore_compare/nanopore_fast5
echo "done"




#
#source /home/liuya/.bash_profile
#
#
#conda activate nanoai
#
#conda env list
#
#conda deactivate
#
#conda env list