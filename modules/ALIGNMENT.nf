/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : ALIGNMENT.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// Align each basecalled outputs
process ALIGNMENT {
	tag "${basecallDir.baseName}"

	input:
	path 	basecallDir
	each 	path(reference_genome)

	output:
	path "${basecallDir.baseName}.alignment", 		optional:true,	emit: alignment
	tuple val(basecallDir.baseName), path ("${basecallDir.baseName}.alignment"),	optional:true,		emit: alignment_tuple

	when:
	params.runAlignment

	shell:
	cores = task.cpus * params.mediumProcTimes
	'''
	mkdir -p !{basecallDir.baseName}.alignment

	## After basecall, we align results to merged, sorted bam, can be for ONT coverage analyses/output bam
	# align FASTQ files to reference genome, write sorted alignments to a BAM file
	minimap2 -t !{cores} -a  -x map-ont \
		!{params.referenceGenome} \
		!{basecallDir}/batch_basecall_combine_fq_*.fq.gz | \
		samtools sort -@ !{cores} -T tmp -o \
			!{basecallDir.baseName}.alignment/!{basecallDir.baseName}_bam.bam &&\
		samtools index -@ !{cores}  !{basecallDir.baseName}.alignment/!{basecallDir.baseName}_bam.bam
	echo "### Samtools alignment DONE"
	'''
}
