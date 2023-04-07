/**
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : GUPPY6.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
**/

// methylation calling for Guppy6
process Guppy6 {
	tag "${fast5Untar.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/guppy6",
		mode: "copy",
		pattern: "${fast5Untar.baseName}_batch_merge_bam_out.bam*",
		enabled: params.outputIntermediate

	publishDir "${params.outdir}/${params.dsname}-run-log",
		mode: "copy", pattern: "*.Guppy6.run.log"

	input:
	path fast5Untar
	each path(reference_genome)
	each path(utils)

	output:
	path "${fast5Untar.baseName}_batch_merge_bam_out.bam*",	emit: guppy_batch_bam_out, optional: true
	path "${fast5Untar.baseName}_guppy6_per_read_batch.tsv.gz", emit: guppy_batch_per_read, optional: true
	path "*.Guppy6.run.log", optional:true,	emit: runlog

	when:
	params.runMethcall && params.runGuppy

	script:
	cores = task.cpus * params.highProcTimes
	samtools_cores = task.cpus * params.mediumProcTimes
	"""
	date; hostname; pwd

	echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"
	if [[ "\${CUDA_VISIBLE_DEVICES:-}" == "" ]] ; then
		echo "Detect no GPU, using CPU commandType"
		commandType='cpu'
		gpuOptions=" "
	else
		echo "Detect GPU, using GPU commandType"
		commandType='gpu'
		gpuOptions="-x auto"
	fi

	if [[ ${params.skipBasecall} == false ]]; then
		indir=${fast5Untar}
	else
		indir=${fast5Untar}/workspace
	fi

	mkdir -p ${fast5Untar.baseName}.methcalled

	## CPU/GPU version command
	guppy_basecaller \
		--config ${params.GUPPY_METHCALL_MODEL} \
		--align_ref ${params.referenceGenome} \
		-i \${indir}  \
		-s ${fast5Untar.baseName}.methcalled \
		--num_callers  ${params.guppy_num_callers} \
		--gpu_runners_per_device  ${params.guppy_gpu_runners_per_device} \
		--cpu_threads_per_caller  ${params.guppy_cpu_threads_per_caller} \
		--bam_out --recursive --compress_fastq \
		--verbose_logs \
		\${gpuOptions} &>> ${params.dsname}.${fast5Untar.baseName}.Guppy6.run.log

	echo "### Guppy6 methylation calling DONE"

	samtools cat -@ ${samtools_cores} \
		-o ${fast5Untar.baseName}_batch_merge_bam_out.bam \
    	\$(find "${fast5Untar.baseName}.methcalled/" "${fast5Untar.baseName}.methcalled/pass/" \
    		-maxdepth 1 -name '*.bam' -type f)

	samtools sort -@ ${samtools_cores}  ${fast5Untar.baseName}_batch_merge_bam_out.bam \
		-o ${fast5Untar.baseName}_batch_merge_bam_out.sort.bam

	rm -f ${fast5Untar.baseName}_batch_merge_bam_out.bam
	mv ${fast5Untar.baseName}_batch_merge_bam_out.sort.bam \
		${fast5Untar.baseName}_batch_merge_bam_out.bam
	samtools index -@ ${samtools_cores}  ${fast5Untar.baseName}_batch_merge_bam_out.bam

	python utils/modbam2bed_extract_read_cpg.py \
		-r ${params.referenceGenome} \
		-i ${fast5Untar.baseName}_batch_merge_bam_out.bam \
    	-o ${fast5Untar.baseName}_guppy6_per_read_batch.tsv \
    	-a ${params.guppy_canon_threshold} -b ${params.guppy_mod_threshold}

    gzip -f  ${fast5Untar.baseName}_guppy6_per_read_batch.tsv

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -rf ${fast5Untar.baseName}.methcalled/{pass,fail}
	fi

	echo "### Guppy6 merged bam_out DONE"
	"""
}


// Combine Guppy6 runs' all results together
process Guppy6Comb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}_guppy6_bam_out_combine.bam*",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}_guppy6_per_read_combine.tsv.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Site_Level-${params.dsname}/*-perSite-cov*.gz"

	input:
	path batch_bam_out_collect
	path batch_per_read_collect
	path reference_genome
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_guppy6_bam_out_combine.bam*", emit: guppy6_combine_bam_out, optional: true
	path "${params.dsname}_guppy6_per_read_combine.tsv.gz", emit: guppy6_combine_tsv, optional: true
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify, optional: true
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify, optional: true

	when:
	batch_bam_out_collect.size() >= 1 && params.runCombine

	script:
	cores = task.cpus * params.highProcTimes
	samtools_cores = task.cpus * params.mediumProcTimes
	"""
	samtools --version

	## combine bam output for Guppy6
	## using cat instead of merge, due to issues for merge large NA12878 chr22 data
	samtools cat -@ ${samtools_cores} -o ${params.dsname}_guppy6_bam_out_combine.bam \
    	\$(find . -maxdepth 1 -name '*_batch_merge_bam_out.bam')

	samtools sort -@ ${samtools_cores}  ${params.dsname}_guppy6_bam_out_combine.bam  \
		-o ${params.dsname}_guppy6_bam_out_combine.sort.bam

	rm -f ${params.dsname}_guppy6_bam_out_combine.bam
	mv ${params.dsname}_guppy6_bam_out_combine.sort.bam \
		${params.dsname}_guppy6_bam_out_combine.bam
	samtools index -@ ${samtools_cores}  ${params.dsname}_guppy6_bam_out_combine.bam

	## combine per read results

	headerFile=\$(ls *_guppy6_per_read_batch.tsv.gz|sort| head -n 1)
	zcat \${headerFile} | head -n 1 | gzip -f > \
		${params.dsname}_guppy6_per_read_combine.tsv.gz

	ls *_guppy6_per_read_batch.tsv.gz | sort | \
		while read infn ; do
			echo \$infn
			zcat \$infn | awk 'NR>1' | gzip -f >> \
				${params.dsname}_guppy6_per_read_combine.tsv.gz
		done

	if [[ ${params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		## sort order: Chr, Start, (End), ID, Strand
		zcat ${params.dsname}_guppy6_per_read_combine.tsv.gz |\
			sort -V -u -k2,2 -k3,3n -k1,1 -k4,4 |\
			gzip -f > ${params.dsname}_guppy6_per_read_combine.sort.tsv.gz
		rm ${params.dsname}_guppy6_per_read_combine.tsv.gz &&\
			mv ${params.dsname}_guppy6_per_read_combine.sort.tsv.gz\
				${params.dsname}_guppy6_per_read_combine.tsv.gz
	fi

	## Unify format output
	echo "### generate read/site level results"
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  Guppy  NANOME\
		${params.dsname}_guppy6_per_read_combine.tsv.gz \
		.  $task.cpus  12  ${params.sort ? true : false}  "${params.chrSet1.replaceAll(',', ' ')}"

	echo "### Guppy6Comb DONE"
	"""
}
