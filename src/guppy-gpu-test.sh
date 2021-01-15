#!/bin/bash
#SBATCH --job-name=deepmod.gpu
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=70g # memory pool for all cores
#SBATCH --time=00:50:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

#set -x

testout=/fastscratch/liuya/nanocompare/deepmod-test-out

rm -rf ${testout}
mkdir -p ${testout}

python /projects/li-lab/yang/tools/latest-version/DeepMod/bin/DeepMod.py detect --wrkBase /fastscratch/liuya/nanocompare/K562-Runs/K562-N50-basecall/7/workspace --Ref /projects/li-lab/reference/hg38/hg38.fasta --outFolder ${testout} --Base C --modfile /projects/li-lab/yang/tools/latest-version/DeepMod/train_deepmod/rnn_conmodC_P100wd21_f7ne1u0_4/mod_train_conmodC_P100wd21_f3ne1u0 --FileID batch_7 --threads 16 --move