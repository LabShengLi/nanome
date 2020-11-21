#!/bin/bash
#SBATCH --job-name=hl60
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 100g # memory pool for all cores
#SBATCH -t 18:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#cd $PBS_O_WORKDIR # change the drectory for the one, from which the job was submitted

set -x

which samtools

HL60_BAM=/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/joined_reps/RRBS/mappingBismark/HL60_joined_RRBS.Read_R1.Rep_3_trimmed_bismark_bt2.bam

HL60_BAM_SORTED=/projects/liuya/results/pkl/nanocompare/HL60_joined_RRBS.Read_R1.Rep_3_trimmed_bismark_bt2.sort.bam

HL60_BAM_BAI=/projects/liuya/results/pkl/nanocompare/HL60_joined_RRBS.Read_R1.Rep_3_trimmed_bismark_bt2.bam.bai

#samtools sort ${HL60_BAM} -o ${HL60_BAM_SORTED}

samtools index ${HL60_BAM_SORTED} ${HL60_BAM_BAI}
