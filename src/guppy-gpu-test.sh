#!/bin/bash
#SBATCH --job-name=basecall.gpu
#SBATCH -p gpu
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem=70g # memory pool for all cores
#SBATCH --time=00:10:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

set -x

testout=/fastscratch/liuya/nanocompare/k562-test-our-guppy-out
rm -rf ${testout}
mkdir -p ${testout}

/projects/li-lab/software/ont-guppy-gpu_4.2.2/bin/guppy_basecaller \
	--input_path /fastscratch/liuya/nanocompare/K562-Runs/K562-N50-sept/1 \
	--save_path ${testout} --config dna_r9.4.1_450bps_hac.cfg --gpu_runners_per_device 16 --num_callers 3 \
	--fast5_out --verbose_logs --device auto
