#!/bin/bash
#SBATCH --job-name=tar
#SBATCH -q batch
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of cores
#SBATCH --mem 50g # memory pool for all cores
#SBATCH -t 02:00:00 # time (D-HH:MM)
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR


set -x

NanoComp --summary APL-N50-basecall.sequencing_summary.txt K562-N50-basecall.sequencing_summary.txt HL60-N50-basecall.sequencing_summary.txt NA19240-N300-basecall.sequencing_summary.txt --names APL K562 HL60 NA19240 --format jpg