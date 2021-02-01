#!/bin/bash
#SBATCH --job-name=cp.na.bascal
##SBATCH -q batch
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem 100g # memory pool for all cores
#SBATCH -t 3-23:59:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

##SBATCH -p gpu
##SBATCH -q inference
##SBATCH --gres=gpu:1

# Run script on Sumner: sbatch --partition=compute -q batch --gres=  rmscript.sh

set -x
#mv /fastscratch/liuya/nanocompare/HL60-Runs/HL60-N50-sept /projects/li-lab/Nanopore_compare/nanopore_fast5
#cp -rf /fastscratch/liuya/nanocompare/APL-Runs/APL-N50-sept /projects/li-lab/Nanopore_compare/nanopore_fast5
#cp -rf /fastscratch/liuya/nanocompare/K562-Runs/K562-N50-sept /projects/li-lab/Nanopore_compare/nanopore_fast5
#cp -rf /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-N300-sept /projects/li-lab/Nanopore_compare/nanopore_fast5

#cp -rf /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-N300-basecall /projects/li-lab/Nanopore_compare/nanopore_fast5

#cp -rf /fastscratch/liuya/nanocompare/NA19240-Runs/NA19240-DeepMod-N300 /projects/li-lab/Nanopore_compare/nanopore_fast5

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