#!/usr/bin/env nextflow

log.info """\
NANOME - NF PIPELINE (v1.0)
by Li Lab at The Jackson Laboratory
http://nanome.jax.org
=================================
dsname			:${params.dsname}
input			:${params.input}
reference_genome	:${params.referenceGenome}
chromSizesFile		:${params.chromSizesFile}
runBasecall		:${params.runBasecall}
runMethcall		:${params.runMethcall}
evaluation		:${params.eval}
=================================
"""
.stripIndent()

projectDir = workflow.projectDir
ch_utils = Channel.fromPath("${projectDir}/utils",  type: 'dir', followLinks: false)
ch_src   = Channel.fromPath("${projectDir}/src",  type: 'dir', followLinks: false)

ch_utils.into{ch_utils1; ch_utils2; ch_utils3; ch_utils4; ch_utils5; ch_utils6; ch_utils7}
ch_src.into{ch_src1; ch_src2}


// We collect all folders of fast5 files, and send into Channels for pipelines
if (params.input.endsWith(".filelist.txt")) { // filelist
	Channel.fromPath( params.input )
		.splitCsv(header: false)
		.map { file(it[0]) }
		.toList()
		.set{ fast5_tar_ch }
	// emit one time for each fast5.tar file
	fast5_tar_ch.flatten().into{fast5_tar_ch1; fast5_tar_ch4}
} else { // single file
	Channel.fromPath( params.input ).into {fast5_tar_ch1; fast5_tar_ch4}
}


// Check all tools work well on the platform
process EnvCheck {
	tag 'EnvCheck'
	errorStrategy 'terminate'
	label 'EnvCheck'

	input:
	file reference_genome_tar from Channel.fromPath(params.reference_genome_tar)

	output:
	file "reference_genome" into reference_genome_ch

	"""
	which tombo
	tombo -v

	which nanopolish
	nanopolish --version

	which megalodon
	megalodon -v

	which deepsignal
	deepsignal

	which DeepMod.py
	DeepMod.py

	which guppy_basecaller
	guppy_basecaller -v

	## Untar to dir reference_genome
	tar -xzf ${reference_genome_tar}

	minimap2 --version
	samtools --version

	echo "### Check env DONE"
	"""
}


//duplicate reference_genome dir to all other processes' input
reference_genome_ch.into{reference_genome_ch1; reference_genome_ch2;
	reference_genome_ch3; reference_genome_ch4; reference_genome_ch5;
	reference_genome_ch6; reference_genome_ch7; reference_genome_ch8}


// Untar of subfolders named 'M1', ..., 'M10', etc.
process Untar {
	tag "${fast5_tar.baseName}"

	input:
	file fast5_tar from fast5_tar_ch4
	each file("*") from ch_utils5

	output:
	file "${fast5_tar.baseName}.untar" into untar_out_ch

	"""
	infn=${fast5_tar}

	mkdir -p untarDir
	if [ "\${infn##*.}" == "tar" ]; then ### deal with tar
		tar -xf ${fast5_tar} -C untarDir
		## move fast5 files in tree folders into a single folder
		mkdir -p untarDir1
		find untarDir -name "*.fast5" -type f -exec mv {} untarDir1/ \\;
	elif [ "\${infn##*.}" == "gz" ]; then ### deal with tar.gz
		tar -xzf ${fast5_tar} -C untarDir
		## move fast5 files in tree folders into a single folder
		mkdir -p untarDir1
		find untarDir -name "*.fast5" -type f -exec mv {} untarDir1/ \\;
	else ### deal with ready folder
		mv  ${fast5_tar} untarDir1
	fi

	## Clean old analyses in input fast5 files
	if [[ "${params.cleanAnalyses}" = true ]] ; then
		echo "### Start cleaning old analysis"
		python utils/clean_old_basecall_in_fast5.py -i untarDir1 --is-indir
	fi

	mv untarDir1 ${fast5_tar.baseName}.untar

	## Clean unused files
	rm -rf untarDir
	"""
}


// Untar output will be used by basecall, megalodon and guppy
untar_out_ch.into{ untar_out_ch1; untar_out_ch2; untar_out_ch3 }


// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
	tag "${fast5_dir.baseName}"

	label 'with_gpus'

	input:
	file fast5_dir from untar_out_ch1
	each file("*") from ch_utils1

	output:
	file "${fast5_dir.baseName}.basecalled" into basecall_out_ch  // try to fix the christina proposed problems
	file "${fast5_dir.baseName}.basecalled/${fast5_dir.baseName}-sequencing_summary.txt" into qc_ch

	when:
	params.runBasecall

	"""
	mkdir -p ${fast5_dir.baseName}.basecalled
	guppy_basecaller --input_path ${fast5_dir} \
		--save_path "${fast5_dir.baseName}.basecalled" \
		--config ${params.GUPPY_BASECALL_MODEL} \
		--num_callers ${params.GuppyNumCallers} \
		--fast5_out \
		--recursive \
		--verbose_logs ${params.GuppyGPUOptions}

	## After basecall, rename and publishe summary file names
	mv ${fast5_dir.baseName}.basecalled/sequencing_summary.txt ${fast5_dir.baseName}.basecalled/${fast5_dir.baseName}-sequencing_summary.txt

	echo "### Basecalled by Guppy DONE"
	"""
}


// Collect and output QC results for basecall
process QCExport {
	publishDir "${params.outputDir}/${params.dsname}-qc-report" , mode: "copy"

	input:
	file flist from qc_ch.collect()

	output:
	file "${params.dsname}-qc-report.tar.gz" into qc_out_ch

	"""
	mkdir -p ${params.dsname}-qc-report
	cp -L -f *-sequencing_summary.txt ${params.dsname}-qc-report/
	tar -czf ${params.dsname}-qc-report.tar.gz ${params.dsname}-qc-report/
	"""
}


// Duplicates basecall outputs
basecall_out_ch
	.into { resquiggle_in_ch; nanopolish_in_ch; deepmod_in_ch }


// methylation calling for Guppy
process Guppy {
	tag "${fast5_dir.baseName}"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/guppy" , mode: "copy"
	label 'with_gpus'

	input:
	file fast5_dir from untar_out_ch2
	each file("*") from ch_utils4
	each file(reference_genome) from reference_genome_ch1

	output:
	file "batch_${fast5_dir.baseName}.bam*" into guppy_methcall_out_ch
	file "batch_${fast5_dir.baseName}.bam.tar.gz" into guppy_methcall_gz_out_ch

	when:
	params.runMethcall && params.runGuppy

	"""
	refGenome=${params.referenceGenome}

	mkdir -p ${fast5_dir.baseName}_methcalled

	guppy_basecaller --input_path ${fast5_dir} \
		--save_path ${fast5_dir.baseName}_methcalled \
		--config dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg \
		--gpu_runners_per_device 32 \
		--num_callers 16 --fast5_out \
		--verbose_logs --device auto

	## install fast5mod
	## pip install fast5mod

	FAST5PATH=${fast5_dir.baseName}_methcalled/workspace
	OUTBAM=batch_${fast5_dir.baseName}.bam
	fast5mod guppy2sam \${FAST5PATH} --reference \${refGenome} \
		--workers 74 --recursive --quiet\
		| samtools sort -@ 8 | samtools view -b -@ 8 > \${OUTBAM}

	samtools sort \${OUTBAM}
	samtools index \${OUTBAM}

	tar -czf batch_${fast5_dir.baseName}.bam.tar.gz batch_${fast5_dir.baseName}.bam*

	## Clean unused files
	rm -rf ${fast5_dir.baseName}_methcalled
	echo "###   Guppy methylation calling DONE"
	"""
}


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${fast5_dir.baseName}"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/megalodon" , mode: "copy"
	label 'with_gpus'

	input:
	file fast5_dir from untar_out_ch3
	each file(reference_genome) from reference_genome_ch2
	each file (megalodonModelTar) from Channel.fromPath(params.megalodon_model_tar)

	output:
	file "batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt.gz" into megalodon_out_ch

	when:
	params.runMethcall && params.runMegalodon

	"""
	refGenome=${params.referenceGenome}

	## Get megalodon model dir
	tar -xzf ${megalodonModelTar}

	megalodon \
		${fast5_dir} \
		--overwrite \
		--outputs per_read_mods mods per_read_refs \
		--guppy-server-path guppy_basecall_server \
		--guppy-config ${params.MEGALODON_MODEL_FOR_GUPPY_CONFIG} \
		--guppy-params "-d ./megalodon_model/ --num_callers 16 --ipc_threads 80" \
		--guppy-timeout ${params.GUPPY_TIMEOUT} \
		--samtools-executable ${params.SAMTOOLS_PATH} \
		--sort-mappings \
		--mappings-format bam \
		--reference \${refGenome} \
		--mod-motif m CG 0 \
		--mod-output-formats bedmethyl wiggle \
		--write-mods-text \
		--write-mod-log-probs \
		--processes ${params.processors} ${params.megalodonGPUOptions}

	### mv megalodon_results/per_read_modified_base_calls.txt batch_${fast5_dir.baseName}.per_read_modified_base_calls.txt
	sed '1d' megalodon_results/per_read_modified_base_calls.txt > batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt
	gzip batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt
	"""
}


// Resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
	tag "${basecallIndir.baseName}"

	input:
	file basecallIndir from resquiggle_in_ch
	//each file(reference_genome_tar) from Channel.fromPath(params.reference_genome_tar)
	each file(reference_genome) from reference_genome_ch3

	output:
	file "${basecallIndir.baseName}_resquiggle_dir" into resquiggle_out_ch
	file "${basecallIndir.baseName}.resquiggle.run.log" into resquiggle_logs

	when:
	params.runMethcall && params.runResquiggle

	"""
	refGenome=${params.referenceGenome}

	### copy basecall workspace files, due to tombo resquiggle modify base folder
	mkdir -p ${basecallIndir.baseName}_resquiggle_dir

	### original basecalled results will be parrallelly used by other processes
	cp -rf ${basecallIndir}/* ${basecallIndir.baseName}_resquiggle_dir/

	### Now set processes=1, to avoid Tombo bug of BrokenPipeError, it is very fast even set to 1.
	### ref: https://github.com/nanoporetech/tombo/issues/139
	tombo resquiggle --dna --processes 1 \
		--corrected-group ${params.resquiggleCorrectedGroup} \
		--basecall-group ${params.BasecallGroupName} --overwrite \
		${basecallIndir.baseName}_resquiggle_dir/workspace \${refGenome} &>  ${basecallIndir.baseName}.resquiggle.run.log

	echo "### Resquiggle DONE"
	"""
}


// duplicate resquiggle results for DeepSignal and Tombo
resquiggle_out_ch.into { deepsignal_in_ch; tombo_in_ch }


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${indir.baseName}"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepsignal" , mode: "copy"
	label 'with_gpus'

	input:
	file indir from deepsignal_in_ch // ttt is [basecallDir, reference_genome]
	each file(deepsignal_model_tar) from Channel.fromPath(params.deepsignel_model_tar)
	each file(reference_genome) from reference_genome_ch4

	output:
	file "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv.gz" into deepsignal_out_ch

	when:
	params.runMethcall && params.runDeepSignal

	"""
	refGenome=${params.referenceGenome}

	tar -xzf ${deepsignal_model_tar}

	deepsignal call_mods --input_path ${indir}/workspace \
		--model_path ./model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt \
		--result_file "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv" \
		--reference_path \${refGenome} \
		--corrected_group ${params.resquiggleCorrectedGroup} \
		--nproc ${params.processors} \
		--is_gpu ${params.DeepSignal_isgpu}

	gzip batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv
	"""
}


// Tombo runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Tombo {
	tag "${resquiggleDir.baseName}"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/tombo" , mode: "copy"

	input:
	each file(reference_genome) from reference_genome_ch5
	each file("*") from ch_utils2
	file resquiggleDir from tombo_in_ch

	output:
	file "batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats.bed.gz" into tombo_out_ch
	file "${resquiggleDir.baseName}.tombo.run.log" into tombo_log_ch

	when:
	params.runMethcall && params.runTombo

	"""
	## using --processes 1 due to tombo bug for BrokenPipeError: [Errno 32] Broken pipe
	## Note 1 is still fast for tombo
	tombo detect_modifications alternative_model \
		--fast5-basedirs ${resquiggleDir}/workspace \
		--dna --standard-log-likelihood-ratio \
		--statistics-file-basename \
		batch_${resquiggleDir.baseName} \
		--per-read-statistics-basename batch_${resquiggleDir.baseName} \
		--alternate-bases CpG \
		--processes 1 \
		--corrected-group ${params.resquiggleCorrectedGroup}  &> ${resquiggleDir.baseName}.tombo.run.log

	python utils/tombo_extract_per_read_stats.py \
		${params.chromSizesFile} \
		"batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats" \
		"batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats.bed"

	gzip batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats.bed
	echo "### Tombo methylation calling DONE"
	"""
}


// DeepMod runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepMod {
	tag "${basecallDir.baseName}"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepmod" , mode: "copy"
	errorStrategy 'ignore'

	// using cpu save costs
	//label 'with_gpus'

	input:
	each file(reference_genome) from reference_genome_ch6
	file basecallDir from deepmod_in_ch

	output:
	file "mod_output/batch_${basecallDir.baseName}_num" into deepmod_out_ch
	file "batch_${basecallDir.baseName}_num.tar.gz" into deepmod_gz_out_ch

	when:
	params.runMethcall && params.runDeepMod

	"""
	wget ${params.DeepModGithub} --no-verbose
	tar -xzf v0.1.3.tar.gz
	DeepModProjectDir="DeepMod-0.1.3"

	refGenome=${params.referenceGenome}

	if [[ "${params.dataType}" = "human" ]] ; then
		mod_cluster=1 ## Human will use cluster model
	else
		mod_cluster=0 ## Not human will skip cluser model
	fi

	DeepMod.py detect \
			--wrkBase ${basecallDir}/workspace --Ref \${refGenome} \
			--Base C --modfile \${DeepModProjectDir}/train_deepmod/${params.deepModModel} \
			--FileID batch_${basecallDir.baseName}_num \
			--threads 16 ${params.DeepModMoveOptions}  \
			--basecall_1d ${params.BasecallGroupName} ###	--mod_cluster 0 is not work

	tar -czf batch_${basecallDir.baseName}_num.tar.gz mod_output/batch_${basecallDir.baseName}_num/
	"""
}


// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Nanopolish {
	tag "${basecallDir.baseName}"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/nanopolish" , mode: "copy"

	input:
	file basecallDir from nanopolish_in_ch
	each file(reference_genome) from reference_genome_ch7
	each file("*") from ch_utils3

	output:
	file "batch_${basecallDir.baseName}.nanopolish.methylation_calls.tsv.gz" into nanopolish_out_ch

	when:
	params.runMethcall && params.runNanopolish

	"""
	refGenome=${params.referenceGenome}

	### We put all fq and bam files into working dir, DO NOT affect the basecall dir
	fastqFile=reads.fq
	fastqNoDupFile="\${fastqFile}.noDups.fq"
	bamFileName="${params.dsname}.batch_${basecallDir.baseName}.sorted.bam"

	echo \${fastqFile}
	echo \${fastqNoDupFile}

	## Do alignment firstly
	touch \${fastqFile}
	for f in \$(ls -1 ${basecallDir}/*.fastq)
	do
		cat \$f >> \$fastqFile
	done

	python utils/nanopore_nanopolish_preindexing_checkDups.py \${fastqFile} \${fastqNoDupFile}

	# Index the raw read with fastq
	nanopolish index -d ${basecallDir}/workspace \${fastqNoDupFile}

	minimap2 -t ${params.processors} -a -x map-ont \${refGenome} \${fastqNoDupFile} | samtools sort -T tmp -o \${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ threads \${bamFileName}
	echo "### samtools finished"
	echo "### Alignment step DONE"

	nanopolish call-methylation -t ${params.processors} -r \${fastqNoDupFile} \
		-b \${bamFileName} -g \${refGenome} > tmp.tsv

	sed '1d' tmp.tsv > batch_${basecallDir.baseName}.nanopolish.methylation_calls.tsv
	gzip batch_${basecallDir.baseName}.nanopolish.methylation_calls.tsv
	echo "### Nanopolish methylation calling DONE"
	"""
}


// prepare combining results
nanopolish_combine_in_ch = nanopolish_out_ch.toList()
deepmod_combine_in_ch=deepmod_out_ch.toList()
megalodon_combine_in_ch = megalodon_out_ch.toList()
tombo_combine_in_ch = tombo_out_ch.toList()
deepsignal_combine_in_ch = deepsignal_out_ch.toList()
guppy_combine_in_ch = guppy_methcall_out_ch.flatten().toList()


// Combine DeepSignal runs' all results together
process DpSigComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
	file x from deepsignal_combine_in_ch

	output:
	file "${params.dsname}.DeepSignal.combine.tsv.gz" into deepsignal_combine_out_ch

	when:
	x.size() >= 1

	"""
	touch ${params.dsname}.DeepSignal.combine.tsv.gz
	cat ${x} > ${params.dsname}.DeepSignal.combine.tsv.gz

	## gzip ${params.dsname}.DeepSignal.combine.tsv
	"""
}


// Combine Tombo runs' all results together
process TomboComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
	file x from tombo_combine_in_ch // list of tombo bed files

	output:
	file "${params.dsname}.Tombo.combine.bed.gz" into tombo_combine_out_ch

	when:
	x.size() >= 1

	"""
	touch ${params.dsname}.Tombo.combine.bed.gz
	cat ${x} > ${params.dsname}.Tombo.combine.bed.gz

	## gzip ${params.dsname}.Tombo.combine.bed
	"""
}


// Combine Guppy runs' all results together
process GuppyComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
	file x from guppy_combine_in_ch
	each file(reference_genome) from reference_genome_ch8

	output:
	file "${params.dsname}.Guppy.combine.tsv.gz" into guppy_combine_out_ch
	file "guppy.total.meth.bam.tar.gz" into guppy_combine_raw_out_ch

	when:
	x.size() >= 1

	"""
	refGenome=${params.referenceGenome}

	find . -name 'batch_*.bam' -maxdepth 1 |
	  parallel -j8 -N4095 -m --files samtools merge -u - |
	  parallel --xargs samtools merge -@8 total.meth.bam {}

	samtools sort total.meth.bam
	samtools index total.meth.bam
	echo "samtool index is done"

	tar -czf guppy.total.meth.bam.tar.gz total.meth.bam*

	## Ref: https://github.com/nanoporetech/medaka/issues/177
	for i in {1..22} X Y
	do
		fast5mod call total.meth.bam \${refGenome} \
			meth.chr\$i.tsv --meth cpg --quiet \
			--regions chr\$i
	done

	cat  meth.*.tsv > total.meth.tsv
	sort -V -k1,1 -k2,2 total.meth.tsv > ${params.dsname}.Guppy.combine.tsv

	gzip ${params.dsname}.Guppy.combine.tsv

	## Clean
	rm -f meth.chr*.tsv total.meth.tsv
	rm -f total.meth.bam*
	echo "### Guppy combine results DONE. ###"
	"""
}


// Combine Megalodon runs' all results together
process MgldnComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
	file x from megalodon_combine_in_ch

	output:
	file "${params.dsname}.Megalodon.combine.bed.gz" into megalodon_combine_out_ch

	when:
	x.size() >= 1

	"""
	> ${params.dsname}.Megalodon.combine.bed.gz

	##sed -n '1p' \${fn} > ${params.dsname}.Megalodon.combine.bed

	for fn in $x
	do
		cat \${fn} >> ${params.dsname}.Megalodon.combine.bed.gz
	done

	##gzip ${params.dsname}.Megalodon.combine.bed
	"""
}


// Combine Nanopolish runs' all results together
process NplshComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
	file x from nanopolish_combine_in_ch

	output:
	file "${params.dsname}.Nanopolish.combine.tsv.gz" into nanopolish_combine_out_ch

	when:
	x.size() >= 1

	"""
	> ${params.dsname}.Nanopolish.combine.tsv.gz

	## sed -n '1p' \${fn} > ${params.dsname}.Nanopolish.combine.tsv

	for fn in $x
	do
		sed '1d' \${fn} >> ${params.dsname}.Nanopolish.combine.tsv.gz
	done

	## gzip ${params.dsname}.Nanopolish.combine.tsv
	"""
}


// Combine DeepMod runs' all results together
process DpmodComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"
	//save costs, DeepMod running slow even in gpu
	//label 'with_gpus'  // cluster model need gpu

	input:
	file x from deepmod_combine_in_ch
	file deepmod_c_tar_file from Channel.fromPath(params.deepmod_ctar)
	each file("*") from ch_utils6

	output:
	file "${params.dsname}.*.combine.bed.gz" into deepmod_combine_out_ch
	file "${params.dsname}.DeepMod.modoutputs.combined.tar.gz" into deepmod_modoutput_combine_ch

	when:
	x.size() >= 1

	"""
	wget ${params.DeepModGithub} --no-verbose
	tar -xzf v0.1.3.tar.gz
	DeepModProjectDir="DeepMod-0.1.3"

	## Untar deepmod cluster-model inputs
	tar -xzf ${deepmod_c_tar_file}

	## Copy all batch results, then summarize site level outputs by chromosome
	mkdir -p indir
	for dx in $x
	do
		mkdir -p indir/\$dx
		cp -rf \$dx/* indir/\$dx/
	done

	python \${DeepModProjectDir}/DeepMod_tools/sum_chr_mod.py \
		indir/ C ${params.dsname}.deepmod ${params.DeepModSumChrSet}

	> ${params.dsname}.DeepModC.combine.bed

	## Note: for ecoli data, no chr*, but N*
	for f in \$(ls -1 indir/${params.dsname}.deepmod.*.C.bed)
	do
	  cat \$f >> ${params.dsname}.DeepModC.combine.bed
	done
	gzip ${params.dsname}.DeepModC.combine.bed

	if [[ "${params.dataType}" == "human" ]] ; then
		## Only apply to human genome
		echo "### For human, apply cluster model of DeepMod"
		python utils/hm_cluster_predict.py \
			indir/${params.dsname}.deepmod \
			./C \
			\${DeepModProjectDir}/train_deepmod/${params.clusterDeepModModel} || true

		> ${params.dsname}.DeepModC_clusterCpG.combine.bed
		for f in \$(ls -1 indir/${params.dsname}.deepmod_clusterCpG.*.C.bed)
		do
		  cat \$f >> ${params.dsname}.DeepModC_clusterCpG.combine.bed
		done

		gzip ${params.dsname}.DeepModC_clusterCpG.combine.bed
	fi

	tar -czf ${params.dsname}.DeepMod.modoutputs.combined.tar.gz indir/
	echo "###DeepMod combine DONE###"
	"""
}


deepsignal_combine_out_ch.concat(tombo_combine_out_ch,megalodon_combine_out_ch, \
	nanopolish_combine_out_ch,deepmod_combine_out_ch.flatten(), guppy_combine_out_ch)
	.toSortedList()
	.into { readlevel_in_ch; sitelevel_in_ch }


// Read level evaluations
process ReadLevelPerf {
	publishDir "${params.outputDir}/${params.dsname}-nanome-analysis" , mode: "copy"

	input: // TODO: I can not sort fileList by name, seems sorted by date????
	file fileList from readlevel_in_ch
	each file("*") from ch_src1

	output:
	file "MethPerf-*" into readlevel_out_ch

	when:
	params.eval && (fileList.size() >= 1)

	"""
	# Sort files
	flist=(\$(ls *.combine.{tsv,bed}.gz))

	echo \${flist[@]}

	deepsignalFile=\$(find . -maxdepth 1 -name '*.DeepSignal.combine.tsv.gz')
	tomboFile=\$(find . -maxdepth 1 -name '*.Tombo.combine.bed.gz')
	nanopolishFile=\$(find . -maxdepth 1 -name '*.Nanopolish.combine.tsv.gz')
	deepmodFile=\$(find . -maxdepth 1 -name '*.DeepModC.combine.bed.gz')
	megalodonFile=\$(find . -maxdepth 1 -name '*.Megalodon.combine.bed.gz')

	export PYTHONPATH=src:\${PYTHONPATH}

	## Read level evaluations
	### python ${workflow.projectDir}/src/nanocompare/read_level_eval.py
	python src/nanocompare/read_level_eval.py \
		--calls DeepSignal:\${deepsignalFile} \
				Tombo:\${tomboFile} \
				Nanopolish:\${nanopolishFile} \
				DeepMod.C:\${deepmodFile} \
				Megalodon:\${megalodonFile} \
		--bgtruth "${params.bgtruthWithEncode}" \
		--runid MethPerf-${params.runid} \
		--dsname ${params.dsname} \
		--min-bgtruth-cov ${params.bgtruth_cov} \
		--report-joined -mpi -o .

	echo "### Read level analysis DONE"
	"""
}


// Site level correlation analysis
process SiteLevelCorr {
	publishDir "${params.outputDir}/${params.dsname}-nanome-analysis" , mode: "copy"

	input:
	file perfDir from readlevel_out_ch
	file fileList from sitelevel_in_ch
	each file("*") from ch_src2

	output:
	file "MethCorr-*" into sitelevel_out_ch

	when:
	params.eval && (fileList.size() >= 1)

	"""
	# Sort file by my self
	flist=(\$(ls *.combine.{tsv,bed}.gz))
	echo \${flist[@]}

	deepsignalFile=\$(find . -maxdepth 1 -name '*.DeepSignal.combine.tsv.gz')
	tomboFile=\$(find . -maxdepth 1 -name '*.Tombo.combine.bed.gz')
	nanopolishFile=\$(find . -maxdepth 1 -name '*.Nanopolish.combine.tsv.gz')
	deepmodFile=\$(find . -maxdepth 1 -name '*.DeepModC_clusterCpG.combine.bed.gz')
	megalodonFile=\$(find . -maxdepth 1 -name '*.Megalodon.combine.bed.gz')

	export PYTHONPATH=src:\${PYTHONPATH}

	## Site level evaluations
	### python ${workflow.projectDir}/src/nanocompare/site_level_eval.py
	python src/nanocompare/site_level_eval.py \
		--calls DeepSignal:\${deepsignalFile} \
				Tombo:\${tomboFile} \
				Nanopolish:\${nanopolishFile} \
				DeepMod.Cluster:\${deepmodFile} \
				Megalodon:\${megalodonFile} \
		--bgtruth "${params.bgtruthWithEncode}" \
		--runid MethCorr-${params.runid} \
		--dsname ${params.dsname} \
		--min-bgtruth-cov ${params.bgtruth_cov} \
		--toolcov-cutoff ${params.tool_cov} \
		--beddir ${perfDir} \
		-o .

	echo "### Site level analysis DONE"
	"""
}
