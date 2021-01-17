#!/bin/bash
#SBATCH --job-name=rm
##SBATCH -q batch
#SBATCH --partition=compute
##SBATCH -p gpu
##SBATCH -q inference
##SBATCH --gres=gpu:1
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 50g # memory pool for all cores
#SBATCH -t 02:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

# Run script on Sumner: sbatch --partition=compute -q batch --gres=  rmscript.sh

set -x

cd  /projects/li-lab/yang/results/results-1st-year
rm -rf 07-02 &

cd /projects/li-lab/ctimage/yang/workspace/covid19-ctimage
rm -rf data results &

cd /projects/li-lab/yang/prev-workspace

rm -rf long_read nanocompa nanocompa_eval nanopolish_test testremote tmp &

wait

echo "RM done"




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