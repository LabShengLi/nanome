/**
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : DEEPSIGNAL.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
**/
// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
// will be deprecated, use DeepSignal v2 instead
process DEEPSIGNAL {
	tag "${indir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/deepsignal",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path indir
	each path(reference_genome)
	each path(deepsignal_model_dir)

	output:
	path "${params.dsname}_deepsignal_batch_${indir.baseName}.*.gz",	emit: deepsignal_tsv, optional: true

	when:
	params.runBasecall && params.runMethcall && params.runDeepSignal1 && ! params.stopDeepSignal

	script:
	cores = task.cpus * params.highProcTimes
	"""
	DeepSignalModelBaseDir="."
	commandType='gpu'
	outFile="${params.dsname}_deepsignal_batch_${indir.baseName}.tsv"

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path "\${DeepSignalModelBaseDir}/${params.DEEPSIGNAL_MODEL_DIR}/${params.DEEPSIGNAL_MODEL}" \
			--result_file \${outFile} \
			--reference_path ${params.referenceGenome} \
			--corrected_group ${params.ResquiggleCorrectedGroup} \
			--nproc $cores \
			--is_gpu no   &>> ${params.dsname}.${indir.baseName}.DeepSignal.run.log
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path "\${DeepSignalModelBaseDir}/${params.DEEPSIGNAL_MODEL_DIR}/${params.DEEPSIGNAL_MODEL}" \
			--result_file \${outFile} \
			--reference_path ${params.referenceGenome} \
			--corrected_group ${params.ResquiggleCorrectedGroup} \
			--nproc $cores \
			--is_gpu yes  &>> ${params.dsname}.${indir.baseName}.DeepSignal.run.log
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi

	gzip -f \${outFile}
	echo "### DeepSignal methylation DONE"
	"""
}


// Combine DeepSignal runs' all results together
process DPSIGCOMB {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_deepsignal_per_read_combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1*.gz"

	input:
	path x
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_deepsignal_per_read_combine.*.gz",	emit: deepsignal_combine_out, optional: true
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify, optional: true
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify, optional: true

	when:
	x.size() >= 1 && params.runCombine

	"""
	touch ${params.dsname}_deepsignal_per_read_combine.tsv.gz
	cat ${x} > ${params.dsname}_deepsignal_per_read_combine.tsv.gz

	if [[ ${params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		## sort order: Chr, Start, (End), ID, Strand
		zcat ${params.dsname}_deepsignal_per_read_combine.tsv.gz |\
			sort -V -u -k1,1 -k2,2n -k5,5 -k3,3 |\
			gzip -f > ${params.dsname}_deepsignal_per_read_combine.sort.tsv.gz
		rm ${params.dsname}_deepsignal_per_read_combine.tsv.gz &&\
			mv ${params.dsname}_deepsignal_per_read_combine.sort.tsv.gz  \
				${params.dsname}_deepsignal_per_read_combine.tsv.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  DeepSignal DeepSignal\
		${params.dsname}_deepsignal_per_read_combine.tsv.gz \
		.  $task.cpus  12 ${params.sort  ? true : false}  "${params.chrSet1.replaceAll(',', ' ')}"
	echo "### DeepSignal combine DONE"
	"""
}
