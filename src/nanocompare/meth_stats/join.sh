#!/bin/bash
#PBS -q batch # batch - up to 72h, long - up to 504h. 
#PBS -l nodes=1:ppn=1   # if I want 32 cores, using machines that has 16 cores, I should specify nodes=2:ppn=16 and be sure 
##PBS -l mem=1000mb # this will ask for 1000mb of RAM; 300gb
#PBS -l walltime=72:00:00
#PBS -M wojciech.rosikiewicz@jax.org
##PBS -m abe # be emailed at the beginning, end and error
##PBS -o /home/rosikw/qsubReports
##PBS -e /home/rosikw/qsubReports
##PBS -t 1-24 # array 
cd $PBS_O_WORKDIR # change the directory for the one, from which the job was submitted


zcat Rep1_R1.ENCFF413KHN.fastq.gz > K562_R1_joined.WGBS.fastq
zcat Rep2_R1.ENCFF336KJH.fastq.gz >> K562_R1_joined.WGBS.fastq
gzip K562_R1_joined.WGBS.fastq

zcat Rep1_R2.ENCFF567DAI.fastq.gz > K562_R2_joined.WGBS.fastq
zcat Rep2_R2.ENCFF585HYM.fastq.gz >> K562_R2_joined.WGBS.fastq
gzip K562_R2_joined.WGBS.fastq


