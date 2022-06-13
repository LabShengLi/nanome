// These tools will be deprecated soon

// methylation calling for Guppy
process Guppy {
	tag "${fast5Untar.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/guppy",
		mode: "copy",
		pattern: "bamfile_${params.dsname}_guppy_fast5mod_batch_${fast5Untar.baseName}_guppy2sam.bam.tar.gz",
		enabled: params.outputIntermediate

	publishDir "${params.outdir}/${params.dsname}_intermediate/guppy",
		mode: "copy",
		pattern: "${params.dsname}_guppy_gcf52ref_batch_${fast5Untar.baseName}.*.gz",
		enabled: params.outputIntermediate && params.runGuppyGcf52ref

	input:
	path fast5Untar
	each path(reference_genome)
	each path(utils)

	output:
	path "${params.dsname}_guppy_fast5mod_batch_${fast5Untar.baseName}_guppy2sam.bam*",	emit: guppy_fast5mod_bam
	path "${params.dsname}_guppy_gcf52ref_batch_${fast5Untar.baseName}.*.gz",	optional: true, emit: guppy_gcf52ref_tsv

	path "bamfile_${params.dsname}_guppy_fast5mod_batch_${fast5Untar.baseName}_guppy2sam.bam.tar.gz",	optional: true,	emit: guppy_fast5mod_bam_gz

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
	guppy_basecaller --input_path \${indir} --recursive\
		--save_path ${fast5Untar.baseName}.methcalled \
		--config ${params.GUPPY_METHCALL_MODEL} \
		--num_callers $task.cpus \
		--fast5_out --compress_fastq\
		--verbose_logs  \${gpuOptions} &>> ${params.dsname}.${fast5Untar.baseName}.Guppy.run.log

	echo "### Guppy methylation calling DONE"

	## Extract guppy methylation-callings
	## Combine fastq
	touch batch_combine_fq.fq.gz

	find "${fast5Untar.baseName}.methcalled/" "${fast5Untar.baseName}.methcalled/pass/"\
		"${fast5Untar.baseName}.methcalled/fail/" -maxdepth 1 -name '*.fastq.gz' -type f\
		-print0 2>/dev/null | \
		while read -d \$'\0' file ; do
			cat \$file >> batch_combine_fq.fq.gz
		done

	## Remove fastq.gz
	find "${fast5Untar.baseName}.methcalled/"   "${fast5Untar.baseName}.methcalled/pass/"\
		"${fast5Untar.baseName}.methcalled/fail/" -maxdepth 1 -name '*.fastq.gz' \
		-type f 2>/dev/null |\
		parallel -j${task.cpus * params.mediumProcTimes} 'rm -f {}'

	if [[ ${params.runGuppyGcf52ref} == true ]]; then
		## gcf52ref ways
		minimap2 -t ${task.cpus * params.mediumProcTimes} -a -x map-ont ${params.referenceGenome} \
			batch_combine_fq.fq.gz | \
			samtools sort -@ ${samtools_cores} \
				-T tmp -o gcf52ref.batch.${fast5Untar.baseName}.bam &&\
			samtools index -@ ${samtools_cores} \
				gcf52ref.batch.${fast5Untar.baseName}.bam
		echo "### gcf52ref minimap2 alignment is done"

		## Modified version, support dir input, not all fast5 files (too long arguments)
		extract_methylation_fast5_support_dir.py \
			-p ${samtools_cores}  ${fast5Untar.baseName}.methcalled/workspace
		echo "### gcf52ref extract to db done"

		## gcf52ref files preparation
		### git clone https://github.com/kpalin/gcf52ref.git
		tar -xzf utils/gcf52ref.tar.gz -C .
		patch gcf52ref/scripts/extract_methylation_from_rocks.py < \
			utils/gcf52ref.patch

		python gcf52ref/scripts/extract_methylation_from_rocks.py \
			-d base_mods.rocksdb \
			-a gcf52ref.batch.${fast5Untar.baseName}.bam \
			-r ${params.referenceGenome} \
			-o tmp.batch_${fast5Untar.baseName}.guppy.gcf52ref_per_read.tsv
		echo "### gcf52ref extract to tsv done"

		awk 'NR>1' tmp.batch_${fast5Untar.baseName}.guppy.gcf52ref_per_read.tsv | gzip -f > \
			${params.dsname}_guppy_gcf52ref_batch_${fast5Untar.baseName}.tsv.gz &&\
			rm -f tmp.batch_${fast5Untar.baseName}.guppy.gcf52ref_per_read.tsv
		echo "### gcf52ref extraction DONE"
	fi

	## fast5mod ways
	FAST5PATH=${fast5Untar.baseName}.methcalled/workspace
	OUTBAM=${params.dsname}_guppy_fast5mod_batch_${fast5Untar.baseName}_guppy2sam.bam

	fast5mod guppy2sam \${FAST5PATH} --reference ${params.referenceGenome} \
		--workers ${samtools_cores}  --recursive --quiet \
		| samtools sort -@   ${samtools_cores}  | \
		samtools view -b -@ ${samtools_cores}  > \${OUTBAM} &&\
		samtools index -@ ${samtools_cores}   \${OUTBAM}

	if [[ "${params.outputIntermediate}" == true ]] ; then
		tar -czf bamfile_\${OUTBAM}.tar.gz \
			\${OUTBAM}*
	fi
	echo "### fast5mod extraction DONE"

	## Clean
	## methcalled folder is no need, keep only gcf52ref's tsv and fast5mod's bam for combine step
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -rf ${fast5Untar.baseName}.methcalled
		rm -f gcf52ref.*.bam gcf52ref.*.bam.bai tmp*.tsv batch_combine_fq.fq.gz
		rm -rf gcf52ref/
		rm -rf base_mods.rocksdb/
		echo "### Clean DONE"
	fi
	echo "### Guppy methcall, extracted by fast5mod and gcf52ref DONE"
	"""
}


// Combine Guppy runs' all results together
process GuppyComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}_guppy_*_per_*_combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}_intermediate/guppy",
		mode: "copy", pattern: "bamfile_${params.dsname}_guppy_fast5mod_combine.bam.tar.gz",
		enabled: params.outputIntermediate

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov*.gz"

	input:
	path x // fast5mod bam inputs
	path y // gcf52ref tsv.gz inputs
	path reference_genome
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_guppy_fast5mod_per_site_combine.*.gz", emit: guppy_fast5mod_combine_out
	path "${params.dsname}_guppy_gcf52ref_per_read_combine.*.gz", optional: true, emit: guppy_gcf52ref_combine_out
	path "bamfile_${params.dsname}_guppy_fast5mod_combine.bam.tar.gz", emit: guppy_combine_raw_out_ch,  optional: true

	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz", optional: true,  emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz", emit: site_unify

	when:
	x.size() >= 1 && params.runCombine

	script:
	cores = task.cpus * params.highProcTimes
	samtools_cores = task.cpus * params.mediumProcTimes
	"""
	if [[ ${params.runGuppyGcf52ref} == true ]] ; then
		## gcf52ref ways
		cat ${params.dsname}_guppy_gcf52ref_batch_*.gz > ${params.dsname}_guppy_gcf52ref_per_read_combine.tsv.gz
		echo "### gcf52ref combine DONE"
	fi

	## fast5mod ways combine
	## find name like batch_*.guppy.fast5mod_guppy2sam.bam*
	find . -maxdepth 1  -name '${params.dsname}_guppy_fast5mod_batch_*_guppy2sam.bam' |
		parallel -j$task.cpus --xargs -v \
		samtools merge -@${samtools_cores} \
		 	${params.dsname}_guppy_fast5mod_combine.bam {}

	### sort is not needed due to merge the sorted bam, ref: http://www.htslib.org/doc/samtools-merge.html
	### samtools sort -@ ${samtools_cores} total.meth.bam
	samtools index -@ ${samtools_cores}  ${params.dsname}_guppy_fast5mod_combine.bam
	echo "Samtool merge and index for fast5mod DONE"

	if [[ ${params.outputIntermediate} == true ]] ; then
		tar -czf bamfile_${params.dsname}_guppy_fast5mod_combine.bam.tar.gz \
			${params.dsname}_guppy_fast5mod_combine.bam*
	fi

	awk '/^>/' ${params.referenceGenome} | awk '{print \$1}' \
		> rf_chr_all_list.txt
	if [[ "${params.dataType1}" == "human" ]] ; then
		echo "### For human, extract chr1-22, X and Y"
		> chr_all_list.txt
		for chrname in ${params.chrSet1.replaceAll(',', ' ')}
		do
			if cat rf_chr_all_list.txt | grep -w ">\${chrname}" -q ; then
				echo \${chrname} >> chr_all_list.txt
			fi
		done
		rm  -f rf_chr_all_list.txt
		echo "### Chromosome list"
		cat chr_all_list.txt

		## Ref: https://github.com/nanoporetech/medaka/issues/177
		parallel -j$cores  -v \
			"fast5mod call  ${params.dsname}_guppy_fast5mod_combine.bam  ${params.referenceGenome} \
				meth.chr_{}.tsv  --meth cpg --quiet --regions {} ; \
				gzip -f meth.chr_{}.tsv" :::: chr_all_list.txt

		touch ${params.dsname}_guppy_fast5mod_per_site_combine.tsv.gz
		for chrname in ${params.chrSet1.replaceAll(',', ' ')}
		do
			if [ -f "meth.chr_\${chrname}.tsv.gz" ]; then
				cat  meth.chr_\${chrname}.tsv.gz >> \
					${params.dsname}_guppy_fast5mod_per_site_combine.tsv.gz
			fi
		done
	else
		echo "### For other genome, params.chrSet1=[${params.chrSet1}]"
		for chrname in ${params.chrSet1.replaceAll(',', ' ')} ;
		do
			fast5mod call ${params.dsname}_guppy_fast5mod_combine.bam  ${params.referenceGenome} \
				meth.chr_\${chrname}.tsv \
				--meth cpg --quiet \
				--regions \${chrname}
			gzip -f  meth.chr_\${chrname}.tsv
		done
		cat meth.chr_*.tsv.gz > ${params.dsname}_guppy_fast5mod_per_site_combine.tsv.gz &&\
			rm -f meth.chr_*.tsv.gz
	fi

	if [[ ${params.deduplicate} == true && ${params.runGuppyGcf52ref} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		## sort order: Chr, Start, (End), ID, Strand
		zcat ${params.dsname}_guppy_gcf52ref_per_read_combine.tsv.gz |\
			sort -V -u -k1,1 -k3,3n -k5,5 -k2,2 |\
			gzip -f > ${params.dsname}_guppy_gcf52ref_per_read_combine.sort.tsv.gz
		rm ${params.dsname}_guppy_gcf52ref_per_read_combine.tsv.gz &&\
			mv ${params.dsname}_guppy_gcf52ref_per_read_combine.sort.tsv.gz \
				${params.dsname}_guppy_gcf52ref_per_read_combine.tsv.gz
	fi

	if [[ ${params.runGuppyGcf52ref} == true ]] ; then
		## Unify format output for read level
		bash utils/unify_format_for_calls.sh \
			${params.dsname}  Guppy Guppy.gcf52ref\
			 ${params.dsname}_guppy_gcf52ref_per_read_combine.tsv.gz \
			.  ${task.cpus}  1  ${params.sort  ? true : false}  "${params.chrSet1.replaceAll(',', ' ')}"
	fi

	## Unify format output for site level
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  Guppy Guppy\
		${params.dsname}_guppy_fast5mod_per_site_combine.tsv.gz \
		.  ${task.cpus}  2  ${params.sort  ? true : false}  "${params.chrSet1.replaceAll(',', ' ')}"

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f meth.chr*.tsv.gz
		rm -f ${params.dsname}_guppy_fast5mod_combine.bam*
	fi
	echo "### Guppy combine DONE"
	"""
}


// Tombo runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Tombo {
	tag "${resquiggleDir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/tombo",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path resquiggleDir
	each path(reference_genome)

	output:
	path "${params.dsname}_tombo_batch_${resquiggleDir.baseName}.*.gz",	emit: tombo_tsv, optional: true

	when:
	params.runMethcall && params.runTombo

	"""
	## Check if there is a BrokenPipeError: [Errno 32] Broken pipe
	## Ref: https://github.com/nanoporetech/tombo/issues/183
	## Note 1 is still fast for tombo
	tombo detect_modifications alternative_model \
		--fast5-basedirs ${resquiggleDir}/workspace \
		--dna\
		--statistics-file-basename ${params.dsname}_tombo_batch_${resquiggleDir.baseName} \
		--per-read-statistics-basename ${params.dsname}_tombo_batch_${resquiggleDir.baseName} \
		--alternate-bases CpG \
		--processes $task.cpus \
		--corrected-group ${params.ResquiggleCorrectedGroup} \
		--multiprocess-region-size ${params.tomboMultiprocessRegionSize} &> \
		${params.dsname}.${resquiggleDir.baseName}.Tombo.run.log

	retry=1
	## while grep -q "BrokenPipeError:" ${resquiggleDir.baseName}.Tombo.run.log
	while ! tail -n 1 ${params.dsname}.${resquiggleDir.baseName}.Tombo.run.log |  grep -q "100%"
	do
		echo "### Found error in tombo detect_modifications, repeat tombo running again!!!"
		tombo detect_modifications alternative_model \
			--fast5-basedirs ${resquiggleDir}/workspace \
			--dna\
			--statistics-file-basename ${params.dsname}_tombo_batch_${resquiggleDir.baseName} \
			--per-read-statistics-basename ${params.dsname}_tombo_batch_${resquiggleDir.baseName} \
			--alternate-bases CpG \
			--processes $task.cpus \
			--corrected-group ${params.ResquiggleCorrectedGroup} \
			--multiprocess-region-size ${params.tomboMultiprocessRegionSize} &> \
			${params.dsname}.${resquiggleDir.baseName}.Tombo.run.log
		retry=\$(( retry+1 ))
		if (( retry >= 5 )); then
			break
		fi
	done

	## if grep -q "BrokenPipeError: \\[Errno 32\\] Broken pipe" ${resquiggleDir.baseName}.Tombo.run.log; then
	if ! tail -n 1 ${params.dsname}.${resquiggleDir.baseName}.Tombo.run.log |  grep -q "100%" ; then
		## Grep the broken pipeline bug for Tombo
		echo "### Tombo seems not finish 100% after retry reached at \${retry} times, please check by yourself, it may be software or genome reference problem."
	else
		## Tombo was ok
		echo "### Tombo log passed, OK"
	fi

	## Tombo lib need h5py lower than 3.0
	## Error may occur with higher version of h5py: AttributeError: 'Dataset' object has no attribute 'value'
	## ref: https://github.com/nanoporetech/tombo/issues/325

	if [ -f "${params.chromSizesFile}" ]; then
		ln -s  ${params.chromSizesFile}  genome.chome.sizes
	else
		cut -f1,2 ${params.referenceGenome}.fai > genome.chome.sizes
	fi

	tombo_extract_per_read_stats.py \
		genome.chome.sizes \
		"${params.dsname}_tombo_batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats" \
		"${params.dsname}_tombo_batch_${resquiggleDir.baseName}.bed" \
		&>> ${params.dsname}.${resquiggleDir.baseName}.Tombo.run.log

	gzip -f ${params.dsname}_tombo_batch_${resquiggleDir.baseName}.bed

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f *.tombo.per_read_stats   *.tombo.stats
	fi
	echo "### Tombo methylation calling DONE"
	"""
}


// Combine Tombo runs' all results together
process TomboComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_tombo_per_read_combine.*.gz",
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
	path "${params.dsname}_tombo_per_read_combine.*.gz",	emit: tombo_combine_out
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify

	when:
	x.size() >= 1 && params.runCombine

	"""
	touch ${params.dsname}_tombo_per_read_combine.bed.gz
	cat ${x} > ${params.dsname}_tombo_per_read_combine.bed.gz

	if [[ ${params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		## sort order: Chr, Start, (End), ID, Strand
		zcat ${params.dsname}_tombo_per_read_combine.bed.gz |\
			sort -V -u -k1,1 -k3,3n -k4,4 -k6,6 |\
			gzip -f > ${params.dsname}_tombo_per_read_combine.sort.bed.gz
		rm ${params.dsname}_tombo_per_read_combine.bed.gz &&\
			mv ${params.dsname}_tombo_per_read_combine.sort.bed.gz  \
				${params.dsname}_tombo_per_read_combine.bed.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  Tombo Tombo\
		${params.dsname}_tombo_per_read_combine.bed.gz \
		.  $task.cpus  12  ${params.sort  ? true : false}  "${params.chrSet1.replaceAll(',', ' ')}"
	echo "### Tombo combine DONE"
	"""
}


// DeepMod runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepMod {
	tag "${basecallDir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/deepmod",
		mode: "copy", pattern: "batch_${basecallDir.baseName}_num.tar.gz",
		enabled: params.outputIntermediate

	input:
	path basecallDir
	each path(reference_genome)

	output:
	path "mod_output/batch_${basecallDir.baseName}_num", emit: deepmod_out
	path "batch_${basecallDir.baseName}_num.tar.gz", emit: deepmod_out_gz,	 optional: true

	when:
	params.runMethcall && params.runDeepMod

	"""
	## Activate nanome conda env if possible
	if [ -d ${params.conda_base_dir} ]; then
		set +u
		source ${params.conda_base_dir}/etc/profile.d/conda.sh
		conda activate ${params.conda_name}
		set -u
	fi

	echo "### Env set ok"
	## Find the model dir
	DeepModTrainModelDir=\$(find \$CONDA_PREFIX -name 'train_deepmod' -type d)
	if [[ \${DeepModTrainModelDir:-} == "" ]]; then
		wget ${params.DeepModGithub} --no-verbose
		tar -xzf v0.1.3.tar.gz
		DeepModTrainModelDir="DeepMod-0.1.3/train_deepmod"
	fi

	if [[ "${params.dataType1}" == "human" ]] ; then
		mod_cluster=1 ## Human will use cluster model
	else
		mod_cluster=0 ## Not human will skip cluser model
	fi

	## Usage ref: https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md#1-how-to-detect-modifications-from-fast5-files
	## DeepMod depends on h5py < 3.0,
	## issue may occur if use a greater version: AttributeError: 'Dataset' object has no attribute 'value'
	rm -rf mod_output
	DeepMod.py detect \
			--wrkBase ${basecallDir}/workspace \
			--outLevel 0\
			--Ref ${params.referenceGenome} \
			--Base C \
			--modfile \${DeepModTrainModelDir}/${params.DEEPMOD_RNN_MODEL} \
			--FileID batch_${basecallDir.baseName}_num \
			--threads ${task.cpus * params.mediumProcTimes} \
			${params.moveOption ? '--move' : ' '} \
			&>> ${params.dsname}.${basecallDir.baseName}.DeepMod.run.log

	if [[ "${params.outputIntermediate}" == true ]] ; then
		tar -czf batch_${basecallDir.baseName}_num.tar.gz mod_output/batch_${basecallDir.baseName}_num/
	fi

	## Clean unused files
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f mod_output/batch_${basecallDir.baseName}_num/rnn.pred.ind.*
		rm -rf mod_output/batch_${basecallDir.baseName}_num/0
	fi
	echo "### DeepMod methylation DONE"
	"""
}


// Combine DeepMod runs' all results together
process DpmodComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}_deepmod_*_per_site_combine.*.gz",
		enabled: params.outputRaw
	publishDir "${params.outdir}/${params.dsname}_intermediate/deepmod",
		mode: "copy", pattern: "${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz",
		enabled: params.outputIntermediate
	publishDir "${params.outdir}/${params.dsname}_intermediate/deepmod",
		mode: "copy", pattern: "${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz",
		enabled: params.outputIntermediate
	publishDir "${params.outdir}/${params.dsname}_intermediate/deepmod",
		mode: "copy", pattern: "${params.dsname}.deepmod.all_batch.C.bed.tar.gz",
		enabled: params.outputIntermediate

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov*.gz"

	input:
	path x
	path deepmod_cfile
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_deepmod_*_per_site_combine.*.gz",	emit: deepmod_combine_out,	 optional: true
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify,	 optional: true

	path "${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz",	emit: deepmod_combine_sum_chrs_mod_ch,	optional: true
	path "${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz",	emit: deepmod_combine_c_cluster_all_chrs_ch,	optional: true
	path "${params.dsname}.deepmod.all_batch.C.bed.tar.gz",	emit: deepmod_combine_all_batch_c_ch,	optional: true

	when:
	x.size() >= 1 && params.runCombine

	"""
	## Activate nanome conda env if possible
	if [ -d ${params.conda_base_dir} ]; then
		set +u
		source ${params.conda_base_dir}/etc/profile.d/conda.sh
		conda activate ${params.conda_name}
		set -u
	fi
	echo "### Env set ok"
	## Find the model dir
	DeepModTrainModelDir=\$(find \$CONDA_PREFIX -name 'train_deepmod' -type d)
	if [[ \${DeepModTrainModelDir:-} == "" ]]; then
		wget ${params.DeepModGithub} --no-verbose
		tar -xzf v0.1.3.tar.gz
		DeepModTrainModelDir="DeepMod-0.1.3/train_deepmod"
	fi

	## Copy all batch results, then summarize site level outputs by chromosome
	## This way will no contamine DeepMod process folders for cache use of nextflow
	mkdir -p indir
	for dx in $x
	do
		mkdir -p indir/\$dx
		cp -rf \$dx/*.C.bed indir/\$dx/
	done

	## merge different runs of modification detection
	## ref: https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md#2-how-to-merge-different-runs-of-modification-detection
	python utils/sum_chr_mod.py \
		indir/ C ${params.dsname}.deepmod ${params.chrSet1.split(' ').join(',')} \
		&>> ${params.dsname}.DpmodComb.run.log

	> ${params.dsname}_deepmod_c_per_site_combine.bed

	## Note: for ecoli data, no pattern for chr*, but N*
	for f in \$(ls -1 indir/${params.dsname}.deepmod.*.C.bed)
	do
	  cat \$f >> ${params.dsname}_deepmod_c_per_site_combine.bed
	done
	gzip -f ${params.dsname}_deepmod_c_per_site_combine.bed

	if [[ "${params.dataType1}" == "human" && "${params.isDeepModCluster}" == "true" ]] ; then
		## Only apply to human genome
		echo "### For human, apply cluster model of DeepMod"

		## Get dir for deepmod cluster-model inputs
		if [[ "${deepmod_cfile}" == *.tar.gz ]] ; then
			tar -xzf ${deepmod_cfile}
		else
			if [[ "${deepmod_cfile}" != "C" ]] ; then
				cp -a  ${deepmod_cfile}   C
			fi
		fi

		## consider modification cluster effect.
		## ref: https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md#3-how-to-consider-modification-cluster-effect
		hm_cluster_predict.py \
			indir/${params.dsname}.deepmod \
			./C \
			\${DeepModTrainModelDir}/${params.DEEPMOD_CLUSTER_MODEL} \
			&>> ${params.dsname}.DpmodComb.run.log

		> ${params.dsname}_deepmod_clusterCpG_per_site_combine.bed
		for f in \$(ls -1 indir/${params.dsname}.deepmod_clusterCpG.*.C.bed)
		do
		  cat \$f >> ${params.dsname}_deepmod_clusterCpG_per_site_combine.bed
		done

		gzip -f ${params.dsname}_deepmod_clusterCpG_per_site_combine.bed

		if [[ "${params.outputIntermediate}" == true ]] ; then
			tar -czf ${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz \
				indir/${params.dsname}.deepmod_clusterCpG.chr*.C.bed
		fi
	fi

	if [[ "${params.outputIntermediate}" == true ]] ; then
		tar -czf ${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz \
			indir/${params.dsname}.deepmod.*.C.bed
		tar -czf ${params.dsname}.deepmod.all_batch.C.bed.tar.gz \
			indir/batch_*_num/mod_pos.*.C.bed
	fi

	if [[ "${params.isDeepModCluster}" == "true" ]] ; then
		callfn=${params.dsname}_deepmod_clusterCpG_per_site_combine.bed.gz
	else
		callfn=${params.dsname}_deepmod_c_per_site_combine.bed.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  DeepMod DeepMod\
		\${callfn} \
		.  $task.cpus  2  ${params.sort ? true : false} "${params.chrSet1.replaceAll(',', ' ')}"\
		&>> ${params.dsname}.DpmodComb.run.log

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -rf indir/
		echo "### Clean DONE"
	fi
	echo "### DeepMod combine DONE"
	"""
}


// Read level unified output, and get METEORE output
process METEORE {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_meteore_*_model_per_read_combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov*.gz"

	input:
	path naonopolish
	path megalodon
	path deepsignal
	path ch_src
	path ch_utils
	path METEOREDir

	output:
	path "${params.dsname}_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.*.gz",	emit: meteore_combine_out, optional: true
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify, 	optional: true
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify,	optional: true

	when:
	params.runMethcall && params.runMETEORE

	"""
	if [[ ${METEOREDir} == *.tar.gz ]] ; then
		## Get METEORE model online
		tar -xzf ${METEOREDir}
	else
		if [[ ${METEOREDir} != ${params.METEOREDirName} ]] ; then
			## rename link folder
			cp -a  ${METEOREDir}  ${params.METEOREDirName}
		fi
	fi

	if [[ -e "${params.METEOREDirName}" ]] ; then
		METEOREDIR="${params.METEOREDirName}"
	else
		echo "### METEORE folder is not correct"
		exit -1
	fi

	## METEORE outputs by combining other tools
	modelContentFileName=${params.dsname}_METEORE_RF_DeepSignal_Megalodon.model_content.tsv
	> \$modelContentFileName
	printf '%s\t%s\n' deepsignal ${deepsignal} >> \$modelContentFileName
	printf '%s\t%s\n' megalodon ${megalodon} >> \$modelContentFileName

	## Degrade sk-learn for METEORE program if needed, it's model load need lower version
	## 0.23.2 version work both for NANOME>=0.23.2 and METEORE<=0.23.2
	## pip install -U scikit-learn==0.23.2
	combineScript="python utils/combination_model_prediction.py"
	## combineScript="combination_model_prediction.py"

	# Use the optimized model
	# Please note this optimized model is reported in METEORE paper, ref: https://github.com/comprna/METEORE#command
	# paper: https://doi.org/10.1101/2020.10.14.340315, text:  random forest (RF) (parameters: max_depth=3 and n_estimator=10)
	# therefore, we output optimized model results
	\${combineScript}\
		-i \${modelContentFileName} \
		-m optimized -b \${METEOREDIR} \
		-o ${params.dsname}_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.tsv.gz\
		&>> ${params.dsname}.METEORE.run.log

	# Use default parameters from the sklearn library (n_estimator = 100 and max_dep = None)
	##\${combineScript}\
	##	-i \${modelContentFileName} -m default -b \${METEOREDIR} \
	##	-o ${params.dsname}.meteore.megalodon_deepsignal_default_rf_model_per_read.combine.tsv.gz\
	##	&>> ${params.dsname}.METEORE.run.log

	if [[ ${params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		## sort order: Chr, Start, (End), ID, Strand
		zcat ${params.dsname}_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.tsv.gz |\
			sort -V -u -k2,2 -k3,3n -k1,1 -k6,6 |\
			gzip -f > ${params.dsname}_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.sort.tsv.gz
		rm ${params.dsname}_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.tsv.gz &&\
			mv ${params.dsname}_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.sort.tsv.gz\
			  	${params.dsname}_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.tsv.gz
	fi

	# Read level and site level output
	if [ -f ${params.dsname}_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.tsv.gz ] ; then
		## Unify format output for read/site level
		bash utils/unify_format_for_calls.sh \
			${params.dsname}  METEORE METEORE\
			${params.dsname}_meteore_deepsignal_megalodon_optimized_rf_model_per_read_combine.tsv.gz \
			.  $task.cpus  12   ${params.sort ? true : false}  "${params.chrSet1.replaceAll(',', ' ')}"\
			&>> ${params.dsname}.METEORE.run.log
	fi
	echo "### METEORE consensus DONE"
	"""
}
