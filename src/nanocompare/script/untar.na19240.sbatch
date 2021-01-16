#!/bin/bash
#SBATCH --job-name=untar.na19240
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 2-14:00 # time (D-HH:MM)
#SBATCH -o %x.%j.figures # STDOUT
#SBATCH -e %x.%j.err # STDERR


NA19240DIR=/projects/li-lab/NanoporeData
outdir=/fastscratch/liuya/nanocompare
mkdir -p ${outdir}

tar -zxf ${NA19240DIR}/NA19240.tar.gz -C ${outdir}
