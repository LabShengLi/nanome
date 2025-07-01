/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : ENVCHECK.nf
 @Software : NANOME project
 @Organization : Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// Check all tools and conditions work well
process ENVCHECK {
	tag "${params.dsname}"
	errorStrategy 'terminate'

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "tools_version_table.tsv", overwrite: true

	input:
	path reference_genome
	path utils
	path rerioDir
	path deepsignalDir

	output:
	path "reference_genome",				emit: reference_genome, optional: true
	path "rerio", 							emit: rerio, optional: true  // used by Megalodon
	path "${params.DEEPSIGNAL_MODEL_DIR}",	emit: deepsignal_model, optional: true // used by DeepSignal v1, will be deprecated
	path "tools_version_table.tsv",			emit: tools_version_tsv, optional: true
	path "basecall_version.txt",			emit: basecall_version_txt, optional: true

	shell:
	'''
	date; hostname; pwd
	echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"

	## Validate nanome container/environment is correct
	bash utils/validate_nanome_container.sh  tools_version_table.tsv

	guppy_basecaller -v |  head -n 1 | python utils/getGuppyVersion.py > basecall_version.txt

	if [[ !{params.runNewTool} == true ]] ; then
		newTools=(!{params.newModuleConfigs.collect{it.name}.join(' ')})
		newToolsVersion=(!{params.newModuleConfigs.collect{it.version}.join(' ')})

		for i in "${!newTools[@]}"; do
			printf "%s\t%s\n" "${newTools[$i]}" "${newToolsVersion[$i]}" >> tools_version_table.tsv
		done
	fi

	## Untar and prepare megalodon model
	if [[ !{params.runMegalodon} == true && !{params.runMethcall} == true ]]; then
		if [[ !{rerioDir} == null* && !{params.rerio} == true ]] ; then
			# Obtain and run R9.4.1, MinION, 5mC CpG model from Rerio
			git clone !{params.rerioGithub}
			rerio/download_model.py rerio/basecall_models/!{params.MEGALODON_MODEL.replace('.cfg', '')}
		elif [[ !{rerioDir} != rerio && -d !{rerioDir} && !{params.rerio} == true ]] ; then
			## rename it to rerio for output channel
			cp  -a !{rerioDir}  rerio
		else
			mkdir -p rerio
			touch rerio/test.txt
		fi
		## Check Rerio model
		ls -lh rerio/
	fi

	if [[ !{params.runBasecall} == true || !{params.runMethcall} == true ]]; then
		## Build dir for reference_genome
		mkdir -p !{params.GENOME_DIR}
		find_dir="$PWD/!{params.GENOME_DIR}"
		if [[ !{reference_genome} == *.tar.gz && -f !{reference_genome}  ]] ; then
			tar -xzf !{reference_genome} -C !{params.GENOME_DIR}
		elif [[ !{reference_genome} == *.tar && -f !{reference_genome} ]] ; then
			tar -xf !{reference_genome} -C !{params.GENOME_DIR}
		elif [[ -d !{reference_genome} ]] ; then
			## for folder, use ln, note this is a link to a folder
			## find_dir=$( readlink -f !{reference_genome} )
			## Copy reference genome, avoid singularity/docker access out data problem
			## cp -f -L !{reference_genome}/*   !{params.GENOME_DIR}/
			find !{reference_genome}/ -maxdepth 1 \\( -type f -o -type l \\) | \
				parallel -j0 'cp -f -L {} !{params.GENOME_DIR}/ || echo "File {} not copied!!!"'
		else
			echo "### ERROR: not recognized reference_genome=!{reference_genome}"
			exit -1
		fi

		# Rename and link reference file
		if [[ ! -z $(find ${find_dir}/ \\( -name '*.fasta' -o -name '*.fasta.gz' \\)  ) ]] ; then
			[[ ! -f !{params.GENOME_DIR}/!{params.GENOME_FN} ]] && \
			 	find ${find_dir} -name '*.fasta*' | \
				 	parallel -j1 -v  'fn={/} ; ln  -s {}   !{params.GENOME_DIR}/${fn/*.fasta/!{params.GENOME_FN}}'
		elif [[ ! -z $(find ${find_dir}/ \\( -name '*.fa' -o -name '*.fa.gz' \\)  ) ]] ; then
			## note here, do not replace .fa, due to example.fa.fai will match all fa pattern
			[[ ! -f !{params.GENOME_DIR}/!{params.GENOME_FN} ]] && \
				find ${find_dir} -name '*.fa*' | \
				 	parallel -j1 -v  'fn={/} ; ln -s {}   !{params.GENOME_DIR}/!{params.GENOME_FN}${fn#*.fa}'
		fi

		## Chrom size file if exists, relink it if name is not same
		[[ ! -f !{params.GENOME_DIR}/!{params.CHROM_SIZE_FN} ]] && \
			find ${find_dir} -name '*.sizes' | head -n 1 |\
				parallel -j1 -v ln  -s -f {} !{params.GENOME_DIR}/!{params.CHROM_SIZE_FN}

		ls -lhiR reference_genome/
	fi
	echo "### Check env DONE"
	'''
}
