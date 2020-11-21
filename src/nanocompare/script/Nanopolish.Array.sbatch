#!/bin/bash
#SBATCH --job-name=methcall.DeepMod
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 5-14:00 # time (D-HH:MM)
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

bash
#source /projects/liuya/workspace/tcgajax/nanocompare/methcall/conda_setup.sh

# default DeepMod model location
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
processors=8

jobindex=$((SLURM_ARRAY_TASK_ID-1))

set -u
basecalledDir=$1
methCalledDir=$2
dsname=$3
set +u

bamFileName="${dsname}.sorted.bam"

## Modify directory for processed files after basecalling:
basecallInput=${basecalledDir}/${jobindex}/workspace/pass/0

# output
methCalledDir=${methCalledDir}/Nanopolish/${jobindex}

mkdir -p ${methCalledDir}

## Join fastq files after basecalling
fastqFile=${basecalledDir}/${jobindex}/workspace/pass/reads.fq

for f in $(ls -1 ${basecalledDir}/${jobindex}/workspace/pass/*.fastq)
do
  cat $f >> $fastqFile
  echo "cat $f >> $fastqFile - COMPLETED"
done

conda activate nanoai

# Indexing
/projects/liuya/tools/nanopolish/nanopolish index -d ${basecallInput} $fastqFile
echo "### index finished"

#### MAPPING:
minimap2 -t $processors -a -x map-ont ${refGenome} $fastqFile > ${methCalledDir}/${bamFileName/sorted.bam/sam}
echo "### minimap2 finished"

samtools view -@ $processors -T ${refGenome} -C -o ${methCalledDir}/${bamFileName/sorted.bam/cram} ${methCalledDir}/${bamFileName/sorted.bam/sam}
echo "### CRAM file created"

samtools sort -@ $processors -O CRAM -o ${methCalledDir}/${bamFileName/bam/cram} ${methCalledDir}/${bamFileName/sorted.bam/cram}
echo "### CRAM file sorted"

samtools index -@ processors ${methCalledDir}/${bamFileName/bam/cram}
echo "### CRAM file indexed"

#### METHYLATION CALLING:
/projects/liuya/tools/nanopolish/nanopolish call-methylation -t $processors -r $fastqFile -b ${methCalledDir}/${bamFileName/bam/cram} -g ${refGenome} > ${methCalledDir}/${dsname}.nanopolish.methylation_calls.tsv
echo "### nanopolish call-methylation finished"
date

conda deactivate

echo "###   nanopolish methylcation calling DONE    ###"
