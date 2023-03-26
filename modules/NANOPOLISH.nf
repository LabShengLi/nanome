/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : NANOPOLISH.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process NANOPOLISH {
	tag "${basecallDir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/nanopolish",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	// path basecallDir
	tuple val(id), path (untarDir), path (basecallDir), path(alignmentDir)
	each path(reference_genome)

	output:
	path "${params.dsname}_nanopolish_batch_${basecallDir.baseName}.*.gz", 	emit: nanopolish_tsv

	when:
	params.runBasecall && params.runMethcall && params.runNanopolish

	script:
	samtools_cores = task.cpus * params.mediumProcTimes
	nanopolish_cores = (task.cpus*params.reduceProcTimes).intValue()

	"""
	## Put all fq and bam files into working dir, DO NOT affect the basecall dir
	## bamFileName="${params.dsname}.batch_${basecallDir.baseName}.sorted.bam"
	bamFileName=\$(find ${alignmentDir}/ -name "*_bam.bam")

	## Do alignment firstly, find the combined fastq file
	fastqFile=\$(find ${basecallDir}/ -name 'batch_basecall_combine_fq_*.fq.gz' -type f)

	## Index, ref: https://github.com/jts/nanopolish#data-preprocessing
	## Index the raw read with fastq, we do not index in basecalled dir, in case of cache can be work
	ln -s \${fastqFile}  \${fastqFile##*/}
	nanopolish index -d ${untarDir}/ \
		-s ${basecallDir}/${basecallDir.baseName}-sequencing_summary.txt \
		\${fastqFile##*/}

	## Calling methylation, ref: https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html#calling-methylation
	## there are segment fault issues, if set -t to a large number or use low memory,
	## ref: https://github.com/jts/nanopolish/issues/872
	## ref: https://github.com/jts/nanopolish/issues/683, https://github.com/jts/nanopolish/issues/580
	nanopolish call-methylation \
		-t ${nanopolish_cores}\
	 	-r \${fastqFile##*/} \
		-b \${bamFileName} -g ${params.referenceGenome} -q cpg | \
		awk 'NR>1' | \
		gzip -f > ${params.dsname}_nanopolish_batch_${basecallDir.baseName}.tsv.gz

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f *.sorted.bam *.sorted.bam.bai
		rm -f *.fq.gz.index*
	fi
	echo "### Nanopolish methylation calling DONE"
	"""
}


// Combine Nanopolish runs' all results together
process NPLSHCOMB {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_nanopolish_per_read_combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov*.gz"

	input:
	path x
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_nanopolish_per_read_combine.tsv.gz",	emit: nanopolish_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify

	//tuple val("Nanopolish"), path ("Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"),
	//	optional:true,		emit: read_unify_map
	//tuple val("Nanopolish"), path ("Site_Level-${params.dsname}/*-perSite-cov*.gz"),
	//	optional:true,		emit: site_unify_map

	when:
	x.size() >= 1 && params.runCombine

	"""
	> ${params.dsname}_nanopolish_per_read_combine.tsv.gz
	cat ${x} > ${params.dsname}_nanopolish_per_read_combine.tsv.gz

	if [[ ${params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		zcat ${params.dsname}_nanopolish_per_read_combine.tsv.gz |\
			sort -V -u -k1,1 -k3,3n -k4,4n -k5,5 -k2,2 |\
			gzip -f > ${params.dsname}_nanopolish_per_read_combine.sort.tsv.gz
		rm ${params.dsname}_nanopolish_per_read_combine.tsv.gz &&\
			mv ${params.dsname}_nanopolish_per_read_combine.sort.tsv.gz \
				${params.dsname}_nanopolish_per_read_combine.tsv.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  Nanopolish Nanopolish \
		${params.dsname}_nanopolish_per_read_combine.tsv.gz \
		.  ${task.cpus}  12 ${params.sort  ? true : false}   "${params.chrSet1.replaceAll(',', ' ')}"

	echo "### Nanopolish combine DONE"
	"""
}
