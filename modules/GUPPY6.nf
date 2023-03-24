/*
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
*/

// methylation calling for Guppy6
process Guppy6 {
	tag "${fast5Untar.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/guppy6",
		mode: "copy",
		pattern: "${fast5Untar.baseName}_batch_merge_bam_out.bam*",
		enabled: params.outputIntermediate

	input:
	path fast5Untar
	each path(reference_genome)
	each path(utils)

	output:
	path "${fast5Untar.baseName}_batch_merge_bam_out.bam*",	emit: guppy_batch_bam_out, optional: true

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

	## sort and combine bam output
	mkdir -p  sort_bam

	find "${fast5Untar.baseName}.methcalled/" "${fast5Untar.baseName}.methcalled/pass/" \
		-maxdepth 1 -name '*.bam' -type f\
		-print0 2>/dev/null | \
		while read -d \$'\0' file ; do
			basefn=\$(basename \$file)
			samtools sort -@ ${samtools_cores} \$file \
				-o  sort_bam/sort_\${basefn}
		done

	samtools merge -@ !{samtools_cores}  ${fast5Untar.baseName}_batch_merge_bam_out.bam  \
		sort_bam/sort*.bam  &&\
		samtools index -@ !{samtools_cores}   ${fast5Untar.baseName}_batch_merge_bam_out.bam

	echo "### Guppy6 merged bam_out DONE"
	"""
}


// Combine Guppy6 runs' all results together
process Guppy6Comb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}_guppy6_merge_bam_out.bam*",
		enabled: params.outputRaw

	input:
	path batch_bam_out
	path reference_genome
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_guppy6_merge_bam_out.bam*", emit: guppy6_combine_bam_out, optional: true

	when:
	batch_bam_out.size() >= 1 && params.runCombine

	script:
	cores = task.cpus * params.highProcTimes
	samtools_cores = task.cpus * params.mediumProcTimes
	"""
	samtools merge -@ !{samtools_cores}  ${params.dsname}_guppy6_merge_bam_out.bam  \
		*_batch_merge_bam_out.bam  &&\
		samtools index -@ !{samtools_cores}  ${params.dsname}_guppy6_merge_bam_out.bam

	echo "### Guppy6Comb DONE"
	"""
}
