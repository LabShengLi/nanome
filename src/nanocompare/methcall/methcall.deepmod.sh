#!/bin/bash
#SBATCH --job-name=deepsignal.methcall
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 250g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# DeepMod methylation call workflow
# Need to populate the parameters into this script
################################################################################
set +x
source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh

set -x

processors=8

numk=$((SLURM_ARRAY_TASK_ID-1))

numKBasecallOutputDir=${basecallOutputDir}/${numk}
echo NumKBasecallOutputDir=${numKBasecallOutputDir}

## Modify directory for processed files after basecalling:
processedFast5DIR=$numKBasecallOutputDir/workspace/pass/0
echo ${processedFast5DIR}

set +x
conda activate nanoai
set -x

## Call methylation from processed fast5 files:
date
python /projects/li-lab/yang/tools/DeepMod/bin/DeepMod.py detect --wrkBase ${processedFast5DIR} --Ref $refGenome --FileID batch_${numk} --modfile $deepModModel --threads $processors --outFolder $methCallsDir

echo "python /projects/li-lab/yang/tools/DeepMod/bin/DeepMod.py detect --wrkBase ${processedFast5DIR} --Ref $refGenome --FileID batch_${numk} --modfile $deepModModel --threads $processors --outFolder $methCallsDir"

echo "###   DeepMod methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   DeepMod Meth-call task is DONE    ###"
