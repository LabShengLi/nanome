#!/bin/bash
#SBATCH --job-name=tombo.methcall
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 250g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

set +x
source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh

set -x

processors=8

numk=$((SLURM_ARRAY_TASK_ID-1))

part=${numk}

numKSeptInputDir=${septInputDir}/${numk}
numKBasecallOutputDir=${basecallOutputDir}/${numk}

echo NumKSeptInputDir=${numKSeptInputDir}
echo NumKBasecallOutputDir=${numKBasecallOutputDir}


## Modify directory for processed files after basecalling:
processedFast5DIR=$numKBasecallOutputDir/workspace/pass/0

set +x
conda activate nanoai
set -x

## Re-squiggle the data:
tombo resquiggle $processedFast5DIR $refGenome --processes $processors --corrected-group $correctedGroup --basecall-group Basecall_1D_000 --overwrite
echo "###   Re-squiggling DONE"

## Call methylation from processed fast5 files:
date
tombo detect_modifications alternative_model --fast5-basedirs $processedFast5DIR --statistics-file-basename $methCallsDir/$analysisPrefix.batch_$part --per-read-statistics-basename $methCallsDir/$analysisPrefix.batch_$part --alternate-bases 5mC --processes $processors --corrected-group $correctedGroup --multiprocess-region-size 10000
date


## postprocess per read methylation calls
python /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/Tombo_extract_per_read_stats.py $chromSizesFile $methCallsDir/$analysisPrefix.batch_$part.5mC.tombo.per_read_stats $methCallsDir/$analysisPrefix.batch_$part.perReadsStats.bed
echo "###   Postprocessing of per read stats files DONE"

set +x
conda deactivate
set -x

echo "tombo detect_modifications alternative_model --fast5-basedirs $processedFast5DIR --statistics-file-basename $methCallsDir/$analysisPrefix.batch_$part --per-read-statistics-basename $methCallsDir/$analysisPrefix.batch_$part --alternate-bases 5mC --processes $processors --corrected-group $correctedGroup --multiprocess-region-size 10000"
echo "###   Tombo methylation calling DONE"


echo "###   Tombo Meth-call DONE    ###"