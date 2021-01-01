#!/bin/bash
#SBATCH --job-name=tombo.methcall
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=250g # memory pool for all cores
#SBATCH --time=19:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# Tombo methylation call workflow
# Need to populate the parameters into this script
################################################################################
set -e
set +x

#source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh
source /home/liuya/.bash_profile
set -x

processors=8

job_index=$((SLURM_ARRAY_TASK_ID-1))

jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

## Modify directory for processed files after basecalling:
processedFast5DIR=${jobkBasecallOutputDir}/workspace/pass/0


set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "targetNum: ${targetNum}"
echo "analysisPrefix: ${analysisPrefix}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "jobkBasecallOutputDir: ${jobkBasecallOutputDir}"
echo "processedFast5DIR: ${processedFast5DIR}"
echo "methCallsDir: ${methCallsDir}"
echo "refGenome: ${refGenome}"
echo "run_resquiggling: ${run_resquiggling}"
echo "correctedGroup: ${correctedGroup}"
echo "chromSizesFile: ${chromSizesFile}"
echo "##################"
set +u

set +x
conda activate nanoai
set -x

if [ "${run_resquiggling}" = true ] ; then
	## Re-squiggle the data:
	tombo resquiggle $processedFast5DIR $refGenome --processes $processors --corrected-group $correctedGroup --basecall-group Basecall_1D_000 --overwrite

	echo "###   Re-squiggling DONE"
fi

## Call methylation from processed fast5 files:
date
tombo detect_modifications alternative_model --fast5-basedirs $processedFast5DIR --statistics-file-basename $methCallsDir/$analysisPrefix.batch_${job_index} --per-read-statistics-basename $methCallsDir/$analysisPrefix.batch_${job_index} --alternate-bases 5mC --processes $processors --corrected-group $correctedGroup --multiprocess-region-size 10000
date

echo "###   Tombo methylation calling DONE"

## Postprocess per read methylation calls for each run complete
python /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/Tombo_extract_per_read_stats.py $chromSizesFile $methCallsDir/$analysisPrefix.batch_${job_index}.5mC.tombo.per_read_stats $methCallsDir/$analysisPrefix.batch_${job_index}.perReadsStats.bed

wc -l $methCallsDir/$analysisPrefix.batch_${job_index}.perReadsStats.bed

echo "###   Postprocessing of per read stats files DONE"

set +x
conda deactivate
set -x

echo "###   Tombo Meth-call task is DONE    ###"
