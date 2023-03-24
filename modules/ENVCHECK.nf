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
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// Check all tools work well
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
	path "${params.DEEPSIGNAL_MODEL_DIR}",	emit: deepsignal_model, optional: true
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

	## Untar and prepare deepsignal model
	if [[ !{params.runDeepSignal} == true && !{params.runMethcall} == true ]]; then
		if [[ !{deepsignalDir} == *.tar.gz ]] ; then
			## Get DeepSignal Model online
			tar -xzf !{deepsignalDir}
		elif [[ !{deepsignalDir} != !{params.DEEPSIGNAL_MODEL_DIR} && -d !{deepsignalDir}  ]] ; then
			## rename it to deepsignal default dir name
			cp  -a !{deepsignalDir}  !{params.DEEPSIGNAL_MODEL_DIR}
		fi
		## Check DeepSignal model
		ls -lh !{params.DEEPSIGNAL_MODEL_DIR}/
	fi

	if [[ !{params.runBasecall} == true || !{params.runMethcall} == true ]]; then
		## Build dir for reference_genome
		mkdir -p reference_genome
		find_dir="$PWD/reference_genome"
		if [[ !{reference_genome} == *.tar.gz && -f !{reference_genome}  ]] ; then
			tar -xzf !{reference_genome} -C reference_genome
		elif [[ !{reference_genome} == *.tar && -f !{reference_genome} ]] ; then
			tar -xf !{reference_genome} -C reference_genome
		elif [[ -d !{reference_genome} ]] ; then
			## for folder, use ln, note this is a symbolic link to a folder
			## find_dir=$( readlink -f !{reference_genome} )
			## Copy reference genome, avoid singularity/docker access out data problem
			cp !{reference_genome}/*   reference_genome/ -f
		else
			echo "### ERROR: not recognized reference_genome=!{reference_genome}"
			exit -1
		fi

		# Rename reference file
		if [[ ! -z $(find ${find_dir}/ \\( -name '*.fasta' -o -name '*.fasta.gz' \\)  ) ]] ; then
			find ${find_dir} -name '*.fasta*' | \
				 parallel -j0 -v  'fn={/} ; ln -s -f  {}   reference_genome/${fn/*.fasta/ref.fasta}'
		elif [[ ! -z $(find ${find_dir}/ \\( -name '*.fa' -o -name '*.fa.gz' \\)  ) ]] ; then
			find ${find_dir} -name '*.fa*' | \
				 parallel -j0 -v  'fn={/} ; ln -s -f  {}   reference_genome/${fn/*.fa/ref.fasta}'
		fi

		## Chrom size file if exists
		find ${find_dir} -name '*.sizes' | \
				parallel -j1 -v ln -s -f {} reference_genome/chrom.sizes

		ls -lh reference_genome/
	fi
	echo "### Check env DONE"
	'''
}

/**
	echo "### Check reference genome and chrSet"
	echo "params.referenceGenome=!{params.referenceGenome}"
	echo "params.chromSizesFile=!{params.chromSizesFile}"
	echo "params.chrSet1=[!{params.chrSet1}]"
	echo "params.dataType1=!{params.dataType1}"
	echo "cpus=!{task.cpus}"
**/