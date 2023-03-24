/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : DEEPSIGNAL2.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
process DEEPSIGNAL2 {
	tag "${resquiggle.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/deepsignal2",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path resquiggle
	path reference_genome
	path ch_src
	path ch_utils
	path deepsignal2_model_file // online model file

	output:
	path "batch_${resquiggle.baseName}_deepsignal2_per_read.tsv.gz",	emit: deepsignal2_batch_per_read, optional: true
	path "batch_${resquiggle.baseName}_deepsignal2_feature.tsv.gz",		emit: deepsignal2_batch_feature, optional: true

	when:
	params.runMethcall && params.runDeepSignal2

	script:
	cores = task.cpus * params.highProcTimes

	shell:
	'''
	set +xu
	. /opt/conda/etc/profile.d/conda.sh
	conda activate /opt/conda/envs/deepsignal2

	set -x
	export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
	export HDF5_PLUGIN_PATH="$CONDA_PREFIX/hdf5/lib/plugin"
	which deepsignal2

	## wget !{params.DEEPSIGNAL2_MODEL_FILE}
	tar -xzf !{deepsignal2_model_file}

	echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"
	if [[ "${CUDA_VISIBLE_DEVICES:-}" == "" ]] ; then
		echo "Detect no GPU, using CPU commandType"
		commandType='cpu'
		gpuOptions=" "
	else
		echo "Detect GPU, using GPU commandType"
		commandType='gpu'
		gpuOptions="--nproc_gpu 1"
	fi

	> batch_!{resquiggle.baseName}_deepsignal2_feature.tsv.gz
	> batch_!{resquiggle.baseName}_deepsignal2_per_read.tsv.gz

	deepsignal2  extract  \
		-i !{resquiggle}/workspace/ \
		-o batch_!{resquiggle.baseName}_deepsignal2_feature.tsv  \
		--corrected_group !{params.ResquiggleCorrectedGroup} \
		--nproc !{cores}  --motifs CG \
		&>> !{params.dsname}.DeepSignal2.run.log

	deepsignal2 call_mods \
		--model_path !{params.DEEPSIGNAL2_MODEL_NAME} \
		--input_path batch_!{resquiggle.baseName}_deepsignal2_feature.tsv \
		--result_file batch_!{resquiggle.baseName}_deepsignal2_per_read.tsv \
		--nproc !{cores}  ${gpuOptions} \
		&>> !{params.dsname}.DeepSignal2.run.log

	cat batch_!{resquiggle.baseName}_deepsignal2_feature.tsv | gzip -f >> \
		batch_!{resquiggle.baseName}_deepsignal2_feature.tsv.gz

	cat batch_!{resquiggle.baseName}_deepsignal2_per_read.tsv| gzip -f >> \
		batch_!{resquiggle.baseName}_deepsignal2_per_read.tsv.gz
	echo "### DeepSignal2 methylation DONE"
	echo "### DeepSignal2 batch DONE"
	'''
}

process DEEPSIGNAL2COMB {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_deepsignal2_per_read_combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Features-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_deepsignal2_feature_combine.*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Site_Level-${params.dsname}/*-perSite-cov1*.gz"

	input:
	path deepsignal2_batch_per_read_collect
	path deepsignal2_batch_feature_collect
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_deepsignal2_per_read_combine.tsv.gz",	emit: deepsignal2_per_read_combine
	path "${params.dsname}_deepsignal2_feature_combine.tsv.gz",	emit: deepsignal2_feature_combine
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify

	when:
	params.runCombine

	script:
	cores = task.cpus * params.highProcTimes

	shell:
	'''
	## combine batches
	> !{params.dsname}_deepsignal2_feature_combine.tsv.gz
	> !{params.dsname}_deepsignal2_per_read_combine.tsv.gz

	cat batch_*deepsignal2_feature.tsv.gz > !{params.dsname}_deepsignal2_feature_combine.tsv.gz
	cat batch_*deepsignal2_per_read.tsv.gz > !{params.dsname}_deepsignal2_per_read_combine.tsv.gz

	if [[ !{params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		## sort order: Chr, Start, (End), ID, Strand
		zcat !{params.dsname}_deepsignal2_per_read_combine.tsv.gz |\
			sort -V -u -k1,1 -k2,2n -k5,5 -k3,3 |\
			gzip -f > !{params.dsname}_deepsignal2_per_read_combine.sort.tsv.gz
		rm !{params.dsname}_deepsignal2_per_read_combine.tsv.gz &&\
			mv !{params.dsname}_deepsignal2_per_read_combine.sort.tsv.gz  \
				!{params.dsname}_deepsignal2_per_read_combine.tsv.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		!{params.dsname}  DeepSignal2 DeepSignal\
		!{params.dsname}_deepsignal2_per_read_combine.tsv.gz \
		.  !{cores}  12 !{params.sort  ? true : false} \
		"!{params.chrSet1.replaceAll(',', ' ')}"
	'''
}
