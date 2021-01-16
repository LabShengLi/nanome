#!/bin/bash
#SBATCH --job-name=deepsignal.methcall
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=250g # memory pool for all cores
#SBATCH --time=1-17:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.figures # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# Nanopolish methylation call workflow
# Need to populate the parameters into this script
################################################################################

set -e
set +x
source ../../utils.common.sh

set -x

job_index=$((SLURM_ARRAY_TASK_ID))
jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

## Modify directory for processed files after basecalling:
processedFast5DIR=${jobkBasecallOutputDir}/workspace
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

date

fastqFile=${jobkBasecallOutputDir}/reads.fq
fastqNoDupFile="${fastqFile}.noDups.fq"

bamFileName="${analysisPrefix}.batch_${job_index}.sorted.bam"


## Do alignment firstly
#if [ "${run_resquiggling}" = true ] ; then

rm -rf ${fastqFile}
touch ${fastqFile}
for f in $(ls -1 ${jobkBasecallOutputDir}/*.fastq)
do
	cat $f >> $fastqFile
	echo "cat $f >> $fastqFile - COMPLETED"
done

rm -rf ${bamFileName} ${fastqNoDupFile}

python nanopore_nanopolish.NA19240_pipeline.step_02.preindexing_checkDups.py ${fastqFile} ${fastqNoDupFile}

# indexing fastq with nanopolish (one needs to have all fast5 reads for that first):
${NanopolishDir}/nanopolish index -d ${processedFast5DIR} ${fastqNoDupFile}

# Reads align to the reference genome:
minimap2 -t ${processors} -a -x map-ont ${refGenome} ${fastqNoDupFile} | samtools sort -T tmp -o ${jobkBasecallOutputDir}/${bamFileName}
echo "### minimap2 finished"

samtools index -@ threads ${jobkBasecallOutputDir}/${bamFileName}
echo "### samtools finished"

echo "### Alignment step DONE"
#fi

# calling methylation:
time ${NanopolishDir}/nanopolish call-methylation -t ${processors} -r ${fastqNoDupFile} -b ${jobkBasecallOutputDir}/${bamFileName} -g ${refGenome} > ${methCallsDir}/${analysisPrefix}.batch_${job_index}.nanopolish.methylation_calls.tsv

echo "### Nanopolish methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   Nanopolish Meth-call task is DONE    ###"
