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

## Export nanopolish path
export PATH=/projects/li-lab/yang/tools/nanopolish:${PATH}
date

if [ "${run_resquiggling}" = true ] ; then
	fastqFile=${numKBasecallOutputDir}/workspace/pass/reads.fq
	rm -rf ${fastqFile}
	touch ${fastqFile}
	for f in $(ls -1 ${numKBasecallOutputDir}/workspace/pass/*.fastq)
	do
		cat $f >> $fastqFile
		echo "cat $f >> $fastqFile - COMPLETED"
	done

	bamFileName="${analysisPrefix}.batch_${numk}.sorted.bam"
	fastqNoDupFile="${fastqFile}.noDups.fq"

	rm -rf ${bamFileName} ${fastqNoDupFile}

	python /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanopore_nanopolish.NA19240_pipeline.step_02.preindexing_checkDups.py ${fastqFile} ${fastqNoDupFile}

	# indexing fastq with nanopolish (one needs to have all fast5 reads for that first):
	nanopolish index -d ${processedFast5DIR} ${fastqNoDupFile}

	# Reads align to the reference genome:
	minimap2 -t ${processors} -a -x map-ont ${refGenome} ${fastqNoDupFile} | samtools sort -T tmp -o ${numKBasecallOutputDir}/${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ threads ${numKBasecallOutputDir}/$bamFileName
	echo "### samtools finished"
fi

# calling methylation:
nanopolish call-methylation -t ${processors} -r ${fastqNoDupFile} -b ${numKBasecallOutputDir}/${bamFileName} -g ${refGenome} > ${methCallsDir}/${analysisPrefix}.batch_${numk}.nanopolish.methylation_calls.tsv

echo "nanopolish call-methylation -t ${processors} -r ${fastqNoDupFile} -b ${numKBasecallOutputDir}/${bamFileName} -g ${refGenome} > ${methCallsDir}/${analysisPrefix}.batch_${numk}.nanopolish.methylation_calls.tsv"

echo "### Nanopolish methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   Nanopolish Meth-call task is DONE    ###"
