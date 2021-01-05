#!/bin/bash
#SBATCH --job-name=deepsignal.methcall
#SBATCH --partition=compute
# SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=250g # memory pool for all cores
#SBATCH --time=17:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# DeepMod methylation call workflow
# Need to populate the parameters into this script
################################################################################
set -e
set +x
#source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh
source /home/liuya/.bash_profile
set -x

#processors=8

job_index=$((SLURM_ARRAY_TASK_ID))
jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

## Modify directory for processed files after basecalling:
processedFast5DIR=${jobkBasecallOutputDir}/workspace/pass/0
echo ${processedFast5DIR}

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
echo "processors: ${processors}"
echo "##################"
set +u


set +x
conda activate nanoai
set -x

## Export nanopolish path
export PATH=/projects/li-lab/yang/tools/nanopolish:${PATH}
date

fastqFile=${jobkBasecallOutputDir}/workspace/pass/reads.fq
fastqNoDupFile="${fastqFile}.noDups.fq"

bamFileName="${analysisPrefix}.batch_${job_index}.sorted.bam"


## Do alignment firstly
if [ "${run_resquiggling}" = true ] ; then

	rm -rf ${fastqFile}
	touch ${fastqFile}
	for f in $(ls -1 ${jobkBasecallOutputDir}/workspace/pass/*.fastq)
	do
		cat $f >> $fastqFile
		echo "cat $f >> $fastqFile - COMPLETED"
	done

	rm -rf ${bamFileName} ${fastqNoDupFile}

	python /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanopore_nanopolish.NA19240_pipeline.step_02.preindexing_checkDups.py ${fastqFile} ${fastqNoDupFile}

	# indexing fastq with nanopolish (one needs to have all fast5 reads for that first):
	nanopolish index -d ${processedFast5DIR} ${fastqNoDupFile}

	# Reads align to the reference genome:
	minimap2 -t ${processors} -a -x map-ont ${refGenome} ${fastqNoDupFile} | samtools sort -T tmp -o ${jobkBasecallOutputDir}/${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ threads ${jobkBasecallOutputDir}/${bamFileName}
	echo "### samtools finished"

	echo "### Alignment step DONE"
fi

# calling methylation:
time nanopolish call-methylation -t ${processors} -r ${fastqNoDupFile} -b ${jobkBasecallOutputDir}/${bamFileName} -g ${refGenome} > ${methCallsDir}/${analysisPrefix}.batch_${job_index}.nanopolish.methylation_calls.tsv

echo "### Nanopolish methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   Nanopolish Meth-call task is DONE    ###"
