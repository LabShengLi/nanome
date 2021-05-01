#!/bin/bash
#SBATCH --job-name=guppy.basecall2
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=170G # memory pool for all cores
#SBATCH --time=23:59:00 # time
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

inputdir=/projects/li-lab/Nanopore_compare/nanopore_fast5/APL-N50-sept

outdir=/fastscratch/liuya/nanocompare/APL_basecalled_long_time
rm -rf $outdir
mkdir -p $outdir

guppy_basecaller --input_path $inputdir --save_path ${outdir} \
    --config dna_r9.4.1_450bps_hac.cfg --num_callers 8 \
    --fast5_out --recursive --verbose_logs --gpu_runners_per_device 8 --device auto
