#!/bin/bash
#SBATCH --job-name=dpmod.motif
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=170g # memory pool for all cores
#SBATCH --time=10:50:00 # time
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR


set -x

DeepModDir=/projects/li-lab/yang/tools/latest-version/DeepMod
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"


outdir=/projects/li-lab/yang/workspace/nano-compare/src/motif

rm -rf ${outdir}/*
mkdir -p ${outdir}


# Ref: https://github.com/WGLab/DeepMod/blob/d822ebaaab3e673cfdb6017637e7603f07cadccd/docs/Usage.md#output-c-in-cpg-motifs-in-a-genome
python ${DeepModDir}/DeepMod_tools/generate_motif_pos.py ${refGenome} ${outdir} C CG 0



#
#echo SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}
#echo SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}
#echo SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}
#echo SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX}
#echo SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN}



#source /home/liuya/.bash_profile
#
#conda env list
#
#conda activate nanoai
#
#conda env list
#
#
#
#conda deactivate