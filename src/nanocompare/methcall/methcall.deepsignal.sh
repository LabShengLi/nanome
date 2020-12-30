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
# DeepSignal methylation call workflow
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

set +x
conda activate nanoai
set -x

if [ "${run_resquiggling}" = true ] ; then
	## Re-squiggle the data:
	tombo resquiggle $processedFast5DIR $refGenome --processes $processors --corrected-group $correctedGroup --basecall-group Basecall_1D_000 --overwrite

	echo "tombo resquiggle $processedFast5DIR $refGenome --processes $processors --corrected-group $correctedGroup --basecall-group Basecall_1D_000 --overwrite"
	echo "###   Re-squiggling DONE"
fi

## Call methylation from processed fast5 files:
date
deepsignal call_mods --input_path $processedFast5DIR --model_path $deepsignalModel --result_file $methCallsDir/$analysisPrefix.batch_${numk}.CpG.deepsignal.call_mods.tsv --reference_path $refGenome --corrected_group $correctedGroup --nproc $processors --is_gpu $isGPU

echo "deepsignal call_mods --input_path $processedFast5DIR --model_path $deepsignalModel --result_file $methCallsDir/$analysisPrefix.batch_${numk}.CpG.deepsignal.call_mods.tsv --reference_path $refGenome --corrected_group $correctedGroup --nproc $processors --is_gpu $isGPU"

echo "###   DeepSignal methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   DeepSignal Meth-call task is DONE    ###"
