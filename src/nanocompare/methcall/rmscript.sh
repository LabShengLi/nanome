#!/bin/bash
#SBATCH --job-name=rm
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 250g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#
#cd /fastscratch/liuya/nanocompare/
#
#ls -d *-N70 *-N60 *-N51
#
#rm -rf *-N70 *-N60 *-N51

find . -name "*.sh" -print0 | while read -d $'\0' file
do
    echo "### processing $file"
done