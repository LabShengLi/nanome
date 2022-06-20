/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : NEWTOOLS.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// NewTool runs
process NewTool {
	tag "${module.name}(${input})"
	container  "${workflow.containerEngine == 'singularity' ? module.container_singularity : workflow.containerEngine == 'docker'? module.container_docker: null}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/${module.name}",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	tuple val (module), path (input)
	each path (genome_path)
	val genome

	output:
	path "${params.dsname}_${module.name}_batch_${input.baseName}.*.gz",  emit: batch_out

	script:
	"""
	echo "### Check input"
	echo ${module.name}
	input=${input}
	genome=${genome}

	echo "### Perform calling"
	${module['cmd']}

	echo "### Rename output"
	echo output=${module['output']}

	if [[ ${module['output']} == *.gz ]] ; then
		mv ${module['output']}  ${params.dsname}_${module.name}_batch_${input.baseName}.tsv.gz
	else
		gzip  -cvf ${module['output']} > ${params.dsname}_${module.name}_batch_${input.baseName}.tsv.gz
	fi
	echo "### NewTool=${module.name} DONE"
	"""
}


// Combine NewTool runs' all results together
process NewToolComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_${module.name}_per_read_combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov*.gz"

	input:
	path batch_in
	val module
	each path(src)

	output:
	path "${params.dsname}_${module.name}_per_read_combine.*.gz",	emit: combine_out
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify

	when:
	params.runCombine

	"""
	touch ${params.dsname}_${module.name}_per_read_combine.tsv.gz

	if [[ ${module.outputHeader} == true ]] ; then
		## remove header before combine
		find . -maxdepth 1 -name '${params.dsname}_${module.name}_batch*.gz' \
	 		-print0 2>/dev/null | \
	 		while read -d \$'\0' file ; do
	 			zcat \$file | \
	 				awk 'NR>1' | \
	 				gzip -f \
	 				>> ${params.dsname}_${module.name}_per_read_combine.tsv.gz
	 		done
	else
		cat ${params.dsname}_${module.name}_batch*.gz\
			> ${params.dsname}_${module.name}_per_read_combine.tsv.gz
	fi

	mkdir -p Read_Level-${params.dsname}
	mkdir -p Site_Level-${params.dsname}

	python src/nanome/nanocompare/newtool_parser.py\
	 	-i  ${params.dsname}_${module.name}_per_read_combine.tsv.gz\
	 	--read-out Read_Level-${params.dsname}/${params.dsname}_${module.name}-perRead-score.tsv.gz \
	 	--site-out Site_Level-${params.dsname}/${params.dsname}_${module.name}-perSite-cov1.sort.bed.gz\
	 	--column-order ${module.outputOrder.join(' ')} \
	 	--score-cols ${module.outputScoreCols.join(' ')}  ${module.logScore ? '--log-score': ' '}\
	 	--chrSet ${chrSet} ${params.deduplicate ? '--deduplicate': ' '} ${params.sort ? '--sort': ' '}

	echo "### NewTool Combine DONE"
	"""
}
