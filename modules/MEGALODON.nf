// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process MEGALODON {
	tag "${fast5Untar.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/megalodon",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path fast5Untar
	each path(reference_genome)
	each path(rerio_dir)

	output:
	path "${params.dsname}_megalodon_batch_${fast5Untar.baseName}.*.gz", emit: megalodon_tsv
	path "${params.dsname}_megalodon_batch_${fast5Untar.baseName}_mod_mappings.bam*", emit: megalodon_mod_mappings

	when:
	params.runMethcall && params.runMegalodon

	script:
	cores = task.cpus * params.mediumProcTimes
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
		gpuOptions="--devices 0"
	fi

	echo "### Guppy dir:"
	which guppy_basecall_server

	if [[ ${params.rerio} == true ]] ; then
		## Rerio model running
		megalodon \
				${fast5Untar}   --overwrite \
				--mod-motif m CG 0 \
				--outputs per_read_mods mods mod_mappings per_read_refs \
				--mod-output-formats bedmethyl wiggle \
				--write-mods-text --write-mod-log-probs \
				--guppy-server-path \$(which guppy_basecall_server) \
				--guppy-config ${params.MEGALODON_MODEL} \
				--guppy-params "-d ./${rerio_dir}/basecall_models/" \
				--guppy-timeout ${params.GUPPY_TIMEOUT} \
				--reference ${params.referenceGenome} \
				--processes $cores \${gpuOptions} \
				&>> ${params.dsname}.${fast5Untar.baseName}.Megalodon.run.log
	else
		## Run Remora model 5mc or 5hmc_5mc
		megalodon ${fast5Untar} --overwrite\
				--guppy-config ${params.GUPPY_BASECALL_MODEL}\
				--remora-modified-bases ${params.remoraModel} fast 0.0.0 ${params.hmc ? "5hmc_5mc" : "5mc"} CG 0\
				--outputs mod_mappings mods per_read_mods \
				--guppy-server-path \$(which guppy_basecall_server) \
				--guppy-timeout ${params.GUPPY_TIMEOUT} \
				--mod-output-formats bedmethyl wiggle \
				--write-mods-text --write-mod-log-probs\
				--reference ${params.referenceGenome}\
				--processes $cores \${gpuOptions} \
				&>> ${params.dsname}.${fast5Untar.baseName}.Megalodon.run.log
	fi

	awk 'NR>1' megalodon_results/per_read_modified_base_calls.txt | gzip -f > \
		${params.dsname}_megalodon_batch_${fast5Untar.baseName}.tsv.gz

	samtools sort -@ ${cores} megalodon_results/mod_mappings.bam -o \
		${params.dsname}_megalodon_batch_${fast5Untar.baseName}_mod_mappings.bam
	samtools index -@ ${cores} ${params.dsname}_megalodon_batch_${fast5Untar.baseName}_mod_mappings.bam

	### Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		### keep guppy server log, due to it may fail when remove that folder, rm -rf megalodon_results
		find megalodon_results/  -maxdepth 1 -type f |\
		 	parallel -j$cores 'rm {}'
	fi
	echo "### Megalodon DONE"
	"""
}


// Combine Megalodon runs' all results together
process MGLDNCOMB {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_megalodon_merge_mod_mappings.bam*",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_megalodon_per_read_combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov*.gz"

	input:
	path x // read-level tsv.gz files
	path y // mod_mappings.bam files
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_megalodon_per_read_combine.*.gz",	emit: megalodon_combine
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify
	path "${params.dsname}_megalodon_merge_mod_mappings.bam*", emit: megalodon_merge_mod_mappings

	when:
	x.size() >= 1  && params.runCombine

	script:
	cores = task.cpus * params.mediumProcTimes
	"""
	> ${params.dsname}_megalodon_per_read_combine.bed.gz
	cat ${x} > ${params.dsname}_megalodon_per_read_combine.bed.gz

	samtools merge -@ ${cores} ${params.dsname}_megalodon_merge_mod_mappings.bam \
		\$(find . -name "${params.dsname}_megalodon_batch_*_mod_mappings.bam")
	samtools index -@ ${cores} ${params.dsname}_megalodon_merge_mod_mappings.bam

	if [[ ${params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		## sort order: Chr, Start, (End), ID, Strand
		zcat ${params.dsname}_megalodon_per_read_combine.bed.gz |\
			sort -V -u -k2,2 -k4,4n -k1,1 -k3,3 |\
			gzip -f > ${params.dsname}_megalodon_per_read_combine.sort.bed.gz
		rm ${params.dsname}_megalodon_per_read_combine.bed.gz &&\
			mv ${params.dsname}_megalodon_per_read_combine.sort.bed.gz \
			 	${params.dsname}_megalodon_per_read_combine.bed.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  Megalodon Megalodon \
		${params.dsname}_megalodon_per_read_combine.bed.gz \
		.  ${task.cpus}  12  ${params.sort  ? true : false}  "${params.chrSet1}"

	echo "### Megalodon combine DONE"
	"""
}
