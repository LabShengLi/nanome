#!/bin/bash
#SBATCH --job-name=rm
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 250g # memory pool for all cores
#SBATCH -t 23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
set -x

#SCRIPT_PATH=${0%/*}
#if [ "$0" != "$SCRIPT_PATH" ] && [ "$SCRIPT_PATH" != "" ]; then
#    cd $SCRIPT_PATH
#fi

# change working path to script path

#cd "$(dirname "$0")"

pwd

ls

source utils.common.sh
