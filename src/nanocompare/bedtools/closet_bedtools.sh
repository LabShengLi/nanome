#!/bin/bash
#SBATCH --job-name=bedtools
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM:SS)
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR


#bedtools intersect -a K562_WGBS_Joined-meth-cov-bgtruth.bed -b K562_WGBS_Joined-meth-cov-Tombo.bed > K562_WGBS_Joined-bgtruth-Tombo.bed
set -x

fn1=/projects/li-lab/yang/results/11-24/K562_WGBS_Joined/K562_WGBS_Joined-meth-cov-Tombo.bed
fn2=/projects/li-lab/yang/results/11-24/K562_WGBS_Joined/K562_WGBS_Joined-meth-cov-bgtruth.bed

outfn=K562-WGBS-Tombo-BGTruth-closest.bed

bedtools sort -i ${fn1} > f1-sorted.bed
bedtools sort -i ${fn2} > f2-sorted.bed

bedtools closest -a f1-sorted.bed -b f2-sorted.bed -d > ${outfn}

rm -f f1-sorted.bed f2-sorted.bed
