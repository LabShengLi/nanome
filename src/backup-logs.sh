#!/bin/bash
#SBATCH --job-name=tar
#SBATCH -q batch
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem 50g # memory pool for all cores
#SBATCH -t 02:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


set -x

cd /fastscratch/liuya/nanocompare

#tar czvf HL60-Runs.log.tar.gz HL60-Runs/*/log HL60-Runs/*/*/log
#tar czvf K562-Runs.log.tar.gz K562-Runs/*/log K562-Runs/*/*/log

tar czvf APL-Runs.log.tar.gz APL-Runs/*/log K562-Runs/*/*/log

mv *.log.tar.gz /projects/li-lab/Nanopore_compare/result/
