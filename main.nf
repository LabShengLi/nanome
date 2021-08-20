#!/usr/bin/env nextflow

log.info """\
NANOME - NF PIPELINE (v1.0)
by Li Lab at The Jackson Laboratory
https://nanome.jax.org
=================================
dsname			:${params.dsname}
input			:${params.input}
reference_genome	:${params.referenceGenome}
runBasecall		:${params.runBasecall}
runMethcall		:${params.runMethcall}
=================================
"""
.stripIndent()

projectDir = workflow.projectDir
ch_utils = Channel.fromPath("${projectDir}/utils",  type: 'dir', followLinks: false)
ch_src   = Channel.fromPath("${projectDir}/src",  type: 'dir', followLinks: false)

// Channel for utils/ and src/ folders
ch_utils.into{ch_utils1; ch_utils2; ch_utils3; ch_utils4; ch_utils5; ch_utils6; ch_utils7; ch_utils8}
ch_src.into{ch_src1; ch_src2; ch_src3; ch_src4}

// Collect all folders of fast5 files, and send into Channels for pipelines
if (params.input.endsWith(".filelist.txt")) { // filelist
	Channel.fromPath( params.input, checkIfExists: true )
		.splitCsv(header: false)
		.map { file(it[0]) }
		.set{ fast5_tar_ch }
} else { // For single file
	Channel.fromPath( params.input, checkIfExists: true ).set {fast5_tar_ch}
}


// Check all tools work well
process EnvCheck {
	tag 'EnvCheck'
	errorStrategy 'terminate'

	input:
	file reference_genome_tar from Channel.fromPath(params.reference_genome_tar)
	each file("*") from ch_utils7

	output:
	file "reference_genome" into reference_genome_ch
	file "gcf52ref" into gcf52ref_code_ch

	"""
	set -x
	which nanopolish
	nanopolish --version

	which megalodon
	megalodon -v

	which deepsignal
	deepsignal

	which guppy_basecaller
	guppy_basecaller -v

	which tombo
	tombo -v

	which DeepMod.py
	DeepMod.py

	which fast5mod
	fast5mod --version

	## Untar to dir reference_genome
	tar -xzf ${reference_genome_tar}

	## gcf52ref file preparation
	### git clone https://github.com/kpalin/gcf52ref.git
	tar xzf utils/gcf52ref.tar.gz -C .
	patch gcf52ref/scripts/extract_methylation_from_rocks.py < utils/gcf52ref.patch

	echo "### Check env DONE"
	"""
}


//duplicate reference_genome dir to all other processes' input
reference_genome_ch.into{
	reference_genome_ch1; reference_genome_ch2;
	reference_genome_ch3; reference_genome_ch4;
	reference_genome_ch5; reference_genome_ch6;
	reference_genome_ch7; reference_genome_ch8;
	reference_genome_ch9 }


// Untar of subfolders named 'M1', ..., 'M10', etc.
process Untar {
	tag "${fast5_tar.baseName}"

	// Disk size is determined by input size, if failed, increase the size
	disk {((fast5_tar.size() * 2.0 as long) >> 30).GB   +  150.GB +   20.GB * task.attempt }

	input:
	file fast5_tar from fast5_tar_ch // using staging, large file suggest firstly using data transfer
	each file("*") from ch_utils5

	output:
	file "${fast5_tar.baseName}.untar" into untar_out_ch
	val "${fast5_tar.size()}" into tar_filesize_ch

	"""
	infn="${fast5_tar}"
	mkdir -p untarTempDir
	if [ "\${infn##*.}" == "tar" ]; then
		### deal with tar
		tar -xf \${infn} -C untarTempDir
		## move fast5 files in tree folders into a single folder
		mkdir -p ${fast5_tar.baseName}.untar
		find untarTempDir -name "*.fast5" -type f -exec mv {} ${fast5_tar.baseName}.untar/ \\;
	elif [ "\${infn##*.}" == "gz" ]; then
		### deal with tar.gz
		tar -xzf \${infn} -C untarTempDir
		## move fast5 files in tree folders into a single folder
		mkdir -p ${fast5_tar.baseName}.untar
		find untarTempDir -name "*.fast5" -type f -exec mv {} ${fast5_tar.baseName}.untar/ \\;
	elif [[ -d ${fast5_tar} ]]; then
		## Copy files, do not change original files such as old analyses data
		cp -rf ${fast5_tar}/* untarTempDir/
		mkdir -p ${fast5_tar.baseName}.untar
		find untarTempDir -name "*.fast5" -type f -exec mv {} ${fast5_tar.baseName}.untar/ \\;
	else
		echo "### Untar error for input=${fast5_tar}"
	fi

	## Clean old analyses in input fast5 files
	if [[ "${params.cleanAnalyses}" = true ]] ; then
		echo "### Start cleaning old analysis"
		python -c 'import h5py; print(h5py.version.info)'
		python utils/clean_old_basecall_in_fast5.py -i ${fast5_tar.baseName}.untar --is-indir --processor ${params.processors * 8}
	fi

	## Clean unused files
	rm -rf untarTempDir

	echo "Total fast5 input files:"
	find ${fast5_tar.baseName}.untar \
		-name "*.fast5" -type f | wc -l
	echo "### Untar DONE"
	"""
}


// Untar output will be used by basecall, megalodon and guppy
untar_out_ch.into{ untar_out_ch1; untar_out_ch2; untar_out_ch3 }
tar_filesize_ch.into{ tar_filesize_ch1; tar_filesize_ch2; tar_filesize_ch3 }


// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
	tag "${fast5_dir.baseName}"

	// Disk size is determined by input size, if failed, increase the size
	disk { (((fast5_tar_size as long)*2.2 as long)>>30).GB   + 100.GB +  20.GB * task.attempt }

	input:
	file fast5_dir from untar_out_ch1
	each file("*") from ch_utils1
	val fast5_tar_size from tar_filesize_ch1
	each file(reference_genome) from reference_genome_ch9

	output:
	file "${fast5_dir.baseName}.basecalled" into basecall_out_ch  // try to fix the christina proposed problems
	file "${fast5_dir.baseName}.basecalled/${fast5_dir.baseName}-sequencing_summary.txt" into qc_ch
	val fast5_tar_size into basecall_filesize_ch
	file "${fast5_dir.baseName}.basecalled.bam" into ont_cov_bam_ch

	when:
	params.runBasecall

	"""
	mkdir -p ${fast5_dir.baseName}.basecalled

	if [[ \${computeName} == "cpu" ]]; then
		## CPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path "${fast5_dir.baseName}.basecalled" \
			--config ${params.GUPPY_BASECALL_MODEL} \
			--num_callers ${params.processors} \
			--fast5_out --recursive \
			--verbose_logs
	elif [[ \${computeName} == "gpu" ]]; then
		## GPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path "${fast5_dir.baseName}.basecalled" \
			--config ${params.GUPPY_BASECALL_MODEL} \
			--num_callers ${params.processors} \
			--fast5_out --recursive \
			--verbose_logs \
			--gpu_runners_per_device ${params.processors} \
			--device auto
	else
		echo "### error value for computeName=\${computeName}"
		exit 255
	fi

	## After basecall, rename and publishe summary file names
	mv ${fast5_dir.baseName}.basecalled/sequencing_summary.txt \
		${fast5_dir.baseName}.basecalled/${fast5_dir.baseName}-sequencing_summary.txt

	## After basecall, we process guppy results for ONT coverage analyses
	# align FASTQ files to reference genome, write sorted alignments to a BAM file
	minimap2 -a -z 600,200 -x map-ont ${params.referenceGenome} ${fast5_dir.baseName}.basecalled/*.fastq \
	    -t ${params.processors*2} > ${fast5_dir.baseName}.basecalled.sam
    echo "Alignment done"

    # Convert the sam file to bam (a binary sam format) using samtools’ view command
    samtools view -u ${fast5_dir.baseName}.basecalled.sam \
        | samtools sort -@ ${params.processors*2} -o ${fast5_dir.baseName}.basecalled.bam --output-fmt BAM
    echo "samtools done"

    # Clean
    rm -f ${fast5_dir.baseName}.basecalled.sam

	echo "### Basecalled by Guppy DONE"
	"""
}


// Collect and output QC results for basecall, and report ONT coverage
process QCExport {
	publishDir "${params.outputDir}/${params.dsname}-qc-report" , mode: "copy"

	input:
	file flist from qc_ch.collect()
	file bamlist from ont_cov_bam_ch.collect()

	output:
	file "${params.dsname}-qc-report.tar.gz" into qc_out_ch
	file "${params.dsname}.coverage.*strand.bed.gz" into ont_cov_ch

	"""
	mkdir -p ${params.dsname}-qc-report
	cp -L -f *-sequencing_summary.txt ${params.dsname}-qc-report/
	tar -czf ${params.dsname}-qc-report.tar.gz ${params.dsname}-qc-report/

    # Merge the bam file
	find . -maxdepth 1 -name "*.bam" | \
	    parallel -j8 -N4095 -m --files samtools merge -u - | \
	    parallel --xargs samtools merge -@ ${params.processors*2} \
	    ${params.dsname}_merged.bam {}";" rm {}

    samtools index ${params.dsname}_merged.bam  -@ ${params.processors*2}
    echo "Samtools merging done!"

    # calculates the sequence coverage at each position/ Reporting genome coverage for all positions in BEDGRAPH format.
    bedtools genomecov -ibam ${params.dsname}_merged.bam -bg -strand + |
        awk '\$4 = \$4 FS "+"' |
        gzip > ${params.dsname}.coverage.positivestrand.bed.gz

    bedtools genomecov -ibam ${params.dsname}_merged.bam -bg -strand - |
        awk '\$4 = \$4 FS "-"' |
        gzip > ${params.dsname}.coverage.negativestrand.bed.gz

    echo "ONT coverage done!"
    echo "### Basecall all DONE"
	"""
}


// Duplicates basecall outputs
basecall_out_ch
	.into { resquiggle_in_ch; nanopolish_in_ch; deepmod_in_ch }


// methylation calling for Guppy
process Guppy {
	tag "${fast5_dir.baseName}"

	// Disk size is determined by input size, if failed, increase the size
	disk { (((fast5_tar_size as long)*2.2 as long) >> 30).GB    + 100.GB +   20.GB * task.attempt }

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/guppy", mode: "copy", pattern: "outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/guppy", mode: "copy", pattern: "batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz"

	input:
	file fast5_dir from untar_out_ch2
	each file(reference_genome) from reference_genome_ch1
	each file("*") from ch_utils4
	each file("*") from gcf52ref_code_ch
	val fast5_tar_size from tar_filesize_ch2

	output:
	file "outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz" into guppy_methcall_gz_out_ch
	file "batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam*" into guppy_methcall_out_ch
	file "batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz" into guppy_gcf52ref_out_ch

	when:
	params.runMethcall && params.runGuppy

	"""
	mkdir -p ${fast5_dir.baseName}.methcalled

	if [[ \${computeName} == "cpu" ]]; then
		## CPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path ${fast5_dir.baseName}.methcalled \
			--config ${params.GUPPY_METHCALL_MODEL} \
			--num_callers ${params.processors} \
			--fast5_out \
			--verbose_logs
	elif [[ \${computeName} == "gpu" ]]; then
		## GPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path ${fast5_dir.baseName}.methcalled \
			--config ${params.GUPPY_METHCALL_MODEL} \
			--num_callers ${params.processors} \
			--fast5_out \
			--verbose_logs \
			--gpu_runners_per_device ${params.processors} \
			--device auto
	else
		echo "### error value for computeName=\${computeName}"
		exit 255
	fi
	echo "### Guppy methylation calling DONE"

	## Extract guppy methylation-callings
	## gcf52ref ways
	minimap2 -t ${params.processors*2} -a -x map-ont ${params.referenceGenome} \
		${fast5_dir.baseName}.methcalled/*.fastq | \
		samtools sort -@ ${params.processors * 2} \
		-T tmp -o gcf52ref.batch.${fast5_dir.baseName}.bam

	samtools index -@ ${params.processors * 2} \
		gcf52ref.batch.${fast5_dir.baseName}.bam
	echo "### gcf52ref minimap2 alignment is done!"

	## Modified version, support dir input, not all fast5 files (too long arguments)
	python utils/extract_methylation_fast5_support_dir.py \
		-p ${params.processors * 2} ${fast5_dir.baseName}.methcalled/workspace
	echo "### gcf52ref extract to db DONE"

	python gcf52ref/scripts/extract_methylation_from_rocks.py \
		-d base_mods.rocksdb \
		-a gcf52ref.batch.${fast5_dir.baseName}.bam \
		-r ${params.referenceGenome} \
		-o tmp.batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv
	echo "### gcf52ref extract to tsv DONE"

	tail -n +2 tmp.batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv | gzip > \
		batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz
	echo "### gcf52ref DONE"

	## fast5mod ways
	FAST5PATH=${fast5_dir.baseName}.methcalled/workspace
	OUTBAM=batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam

	fast5mod guppy2sam \${FAST5PATH} --reference ${params.referenceGenome} \
		--workers 74 --recursive --quiet \
		| samtools sort -@ ${params.processors * 2} | \
		samtools view -b -@ ${params.processors * 2} > \${OUTBAM}

	samtools index -@ ${params.processors * 2}  \${OUTBAM}

	tar -czf outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz \
		batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam*
	echo "### fast5mod DONE"
	echo "### Guppy fast5mod and gcf52ref DONE"
	"""
}


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${fast5_dir.baseName}"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/megalodon", mode: "copy"

	// Disk size is determined by input size, if failed, increase the size
	disk { (((fast5_tar_size as long)*2.2 as long) >> 30).GB    + 100.GB +   20.GB * task.attempt }

	input:
	file fast5_dir from untar_out_ch3
	val fast5_tar_size from tar_filesize_ch3
	each file(reference_genome) from reference_genome_ch2
	each file (megalodonModelTar) from Channel.fromPath(params.megalodon_model_tar)

	output:
	file "batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt.gz" into megalodon_out_ch

	when:
	params.runMethcall && params.runMegalodon

	"""
	## Get megalodon model dir
	tar -xzf ${megalodonModelTar}

	if [[ \${computeName} == "cpu" ]]; then
		## CPU version command
		## Ref: https://github.com/nanoporetech/megalodon
		megalodon \
			${fast5_dir} \
			--overwrite \
			--outputs per_read_mods mods per_read_refs \
			--guppy-server-path guppy_basecall_server \
			--guppy-config ${params.MEGALODON_MODEL_FOR_GUPPY_CONFIG} \
			--guppy-params "-d ./megalodon_model/ --num_callers ${params.processors} --ipc_threads 80" \
			--guppy-timeout ${params.GUPPY_TIMEOUT} \
			--samtools-executable ${params.SAMTOOLS_PATH} \
			--sort-mappings \
			--mappings-format bam \
			--reference ${params.referenceGenome} \
			--mod-motif m CG 0 \
			--mod-output-formats bedmethyl wiggle \
			--write-mods-text \
			--write-mod-log-probs \
			--processes ${params.processors * 2}
	elif [[ \${computeName} == "gpu" ]]; then
		## GPU version command
		## Ref: https://github.com/nanoporetech/megalodon
		megalodon \
			${fast5_dir} \
			--overwrite \
			--outputs per_read_mods mods per_read_refs \
			--guppy-server-path guppy_basecall_server \
			--guppy-config ${params.MEGALODON_MODEL_FOR_GUPPY_CONFIG} \
			--guppy-params "-d ./megalodon_model/ --num_callers ${params.processors} --ipc_threads 80" \
			--guppy-timeout ${params.GUPPY_TIMEOUT} \
			--samtools-executable ${params.SAMTOOLS_PATH} \
			--sort-mappings \
			--mappings-format bam \
			--reference ${params.referenceGenome} \
			--mod-motif m CG 0 \
			--mod-output-formats bedmethyl wiggle \
			--write-mods-text \
			--write-mod-log-probs \
			--processes ${params.processors * 2} \
			--devices 0
	else
		echo "### error value for computeName=\${computeName}"
		exit 255
	fi

	### mv megalodon_results/per_read_modified_base_calls.txt batch_${fast5_dir.baseName}.per_read_modified_base_calls.txt
	tail -n +2 megalodon_results/per_read_modified_base_calls.txt | gzip > \
		batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt.gz
	echo "### Megalodon DONE"
	"""
}


// Resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
// TODO: reduce processors to run resquiggle
process Resquiggle {
	tag "${basecallIndir.baseName}"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/resquiggle" , mode: "copy", pattern: "${basecallIndir.baseName}.resquiggle.run.log"

	// Disk size is determined by input size, if failed, increase the size
	disk { (((file_size as long)*2.2 as long) >> 30).GB    + 100.GB +   20.GB * task.attempt }

	input:
	file basecallIndir from resquiggle_in_ch
	each file(reference_genome) from reference_genome_ch3
	val file_size from basecall_filesize_ch

	output:
	file "${basecallIndir.baseName}.resquiggle" into resquiggle_out_ch
	file "${basecallIndir.baseName}.resquiggle.run.log" into resquiggle_logs

	when:
	params.runMethcall && params.runResquiggle && !params.filterGPUTaskRuns

	"""
	### copy basecall workspace files, due to tombo resquiggle modify base folder
	mkdir -p ${basecallIndir.baseName}.resquiggle

	### original basecalled results will be parrallelly used by other processes
	cp -rf ${basecallIndir}/* ${basecallIndir.baseName}.resquiggle/

	### Need to check Tombo bug of BrokenPipeError, it is very fast even set to 1.
	### ref: https://github.com/nanoporetech/tombo/issues/139
	### ref: https://nanoporetech.github.io/tombo/resquiggle.html?highlight=processes
	tombo resquiggle --dna \
		--processes ${params.processors * 2} \
		--corrected-group ${params.resquiggleCorrectedGroup} \
		--basecall-group ${params.BasecallGroupName} \
		--overwrite \
		${basecallIndir.baseName}.resquiggle/workspace \
		${params.referenceGenome} &> \
		${basecallIndir.baseName}.resquiggle.run.log
	echo "### Resquiggle DONE"
	"""
}


// duplicate resquiggle results for DeepSignal and Tombo
resquiggle_out_ch.into { deepsignal_in_ch; tombo_in_ch }


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${indir.baseName}"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepsignal" , mode: "copy"

	input:
	file indir from deepsignal_in_ch
	each file(deepsignal_model_tar) from Channel.fromPath(params.deepsignal_model_tar)
	each file(reference_genome) from reference_genome_ch4

	output:
	file "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv.gz" into deepsignal_out_ch

	when:
	params.runMethcall && params.runDeepSignal && !params.filterGPUTaskRuns

	"""
	tar -xzf ${deepsignal_model_tar}

	if [[ \${computeName} == "cpu" ]]; then
		## CPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path ./model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt \
			--result_file "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv" \
			--reference_path ${params.referenceGenome} \
			--corrected_group ${params.resquiggleCorrectedGroup} \
			--nproc ${params.processors  * params.deepLearningProcessorTimes} \
			--is_gpu no
	elif [[ \${computeName} == "gpu" ]]; then
		## GPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path ./model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt \
			--result_file "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv" \
			--reference_path ${params.referenceGenome} \
			--corrected_group ${params.resquiggleCorrectedGroup} \
			--nproc ${params.processors  * params.deepLearningProcessorTimes} \
			--is_gpu yes
	else
		echo "### error value for computeName=\${computeName}"
		exit 255
	fi

	gzip batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv
	echo "### DeepSignal methylation DONE"
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
	params.runMethcall && params.runTombo && !params.filterGPUTaskRuns

	"""
	## Check if there is a BrokenPipeError: [Errno 32] Broken pipe
	## Ref: https://github.com/nanoporetech/tombo/issues/183
	## Note 1 is still fast for tombo
	tombo detect_modifications alternative_model \
		--fast5-basedirs ${resquiggleDir}/workspace \
		--dna --standard-log-likelihood-ratio \
		--statistics-file-basename batch_${resquiggleDir.baseName} \
		--per-read-statistics-basename batch_${resquiggleDir.baseName} \
		--alternate-bases CpG \
		--processes 1 \
		--corrected-group ${params.resquiggleCorrectedGroup} \
		--multiprocess-region-size 1000 &> \
		${resquiggleDir.baseName}.tombo.run.log

	if grep -q "BrokenPipeError: \\[Errno 32\\] Broken pipe" ${resquiggleDir.baseName}.tombo.run.log; then
		## Grep the broken pipeline bug for Tombo
		echo "### Tombo bug occur, may need retry to solve it"
		exit 32
	else  ## Tombo was ok
		echo "### Tombo log passed"
	fi

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
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepmod" , mode: "copy", pattern: "batch_${basecallDir.baseName}_num.tar.gz"

	input:
	each file(reference_genome) from reference_genome_ch6
	file basecallDir from deepmod_in_ch

	output:
	file "mod_output/batch_${basecallDir.baseName}_num" into deepmod_out_ch
	file "batch_${basecallDir.baseName}_num.tar.gz" into deepmod_gz_out_ch

	when:
	params.runMethcall && params.runDeepMod && !params.filterGPUTaskRuns

	"""
	wget ${params.DeepModGithub} --no-verbose
	tar -xzf v0.1.3.tar.gz
	DeepModProjectDir="DeepMod-0.1.3"

	if [[ "${params.dataType}" = "human" ]] ; then
		mod_cluster=1 ## Human will use cluster model
	else
		mod_cluster=0 ## Not human will skip cluser model
	fi

	DeepMod.py detect \
			--wrkBase ${basecallDir}/workspace \
			--Ref ${params.referenceGenome} \
			--Base C \
			--modfile \${DeepModProjectDir}/train_deepmod/${params.deepModModel} \
			--FileID batch_${basecallDir.baseName}_num \
			--threads ${params.processors * params.deepLearningProcessorTimes} ${params.DeepModMoveOptions}

	tar -czf batch_${basecallDir.baseName}_num.tar.gz mod_output/batch_${basecallDir.baseName}_num/
	echo "### DeepMod methylation DONE"
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
	params.runMethcall && params.runNanopolish && !params.filterGPUTaskRuns

	"""
	### Put all fq and bam files into working dir, DO NOT affect the basecall dir
	fastqFile=reads.fq
	fastqNoDupFile="\${fastqFile}.noDups.fq"
	bamFileName="${params.dsname}.batch_${basecallDir.baseName}.sorted.bam"

	## Do alignment firstly
	touch \${fastqFile}
	for f in \$(ls -1 ${basecallDir}/*.fastq)
	do
		cat \$f >> \$fastqFile
	done

	python utils/nanopore_nanopolish_preindexing_checkDups.py \${fastqFile} \${fastqNoDupFile}

	# Index the raw read with fastq
	nanopolish index -d ${basecallDir}/workspace \${fastqNoDupFile}

	minimap2 -t ${params.processors} -a -x map-ont ${params.referenceGenome} \${fastqNoDupFile} | \
		samtools sort -@ ${params.processors} -T tmp -o \${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ ${params.processors}  \${bamFileName}
	echo "### samtools finished"
	echo "### Alignment step DONE"

	nanopolish call-methylation -t ${params.processors} -r \${fastqNoDupFile} \
		-b \${bamFileName} -g ${params.referenceGenome} > tmp.tsv

	tail -n +2 tmp.tsv | gzip > batch_${basecallDir.baseName}.nanopolish.methylation_calls.tsv.gz
	echo "### Nanopolish methylation calling DONE"
	"""
}


// prepare combining results
nanopolish_combine_in_ch = nanopolish_out_ch.toList()
deepmod_combine_in_ch=deepmod_out_ch.toList()
megalodon_combine_in_ch = megalodon_out_ch.toList()
tombo_combine_in_ch = tombo_out_ch.toList()
deepsignal_combine_in_ch = deepsignal_out_ch.toList()
guppy_combine_in_ch = guppy_methcall_out_ch.collect()


// Combine DeepSignal runs' all results together
process DpSigComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
	file x from deepsignal_combine_in_ch

	output:
	file "${params.dsname}.deepsignal.per_read.combine.tsv.gz" into deepsignal_combine_out_ch

	when:
	x.size() >= 1 && params.runCombine

	"""
	touch ${params.dsname}.deepsignal.per_read.combine.tsv.gz
	cat ${x} > ${params.dsname}.deepsignal.per_read.combine.tsv.gz
	echo "### DeepSignal combine DONE"
	"""
}


// Combine Tombo runs' all results together
process TomboComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
	file x from tombo_combine_in_ch // list of tombo bed files

	output:
	file "${params.dsname}.tombo.per_read.combine.bed.gz" into tombo_combine_out_ch

	when:
	x.size() >= 1 && params.runCombine

	"""
	touch ${params.dsname}.tombo.per_read.combine.bed.gz
	cat ${x} > ${params.dsname}.tombo.per_read.combine.bed.gz
	echo "### Tombo combine DONE"
	"""
}


// Combine Guppy runs' all results together
process GuppyComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy", pattern: "${params.dsname}.guppy.*.combine.tsv.gz"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/guppy" , mode: "copy", pattern: "${params.dsname}.guppy_fast5mod.combined.bam.tar.gz"

	input:
	file x from guppy_combine_in_ch
	file y from guppy_gcf52ref_out_ch.collect()
	file reference_genome from reference_genome_ch8

	output:
	file "${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz" into guppy_fast5mod_combine_out_ch
	file "${params.dsname}.guppy.gcf52ref_per_read.combine.tsv.gz" into guppy_gcf52ref_combine_out_ch
	file "${params.dsname}.guppy_fast5mod.combined.bam.tar.gz" into guppy_combine_raw_out_ch

	when:
	x.size() >= 1 && params.runCombine

	"""
	## gcf52ref ways
	cat batch_*.guppy.gcf52ref_per_read.tsv.gz > ${params.dsname}.guppy.gcf52ref_per_read.combine.tsv.gz
	echo "### gcf52ref combine DONE"

	## fast5mod ways combine
	## find name like batch_\${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam*
	find . -name 'batch_*.guppy.fast5mod_guppy2sam.bam' -maxdepth 1 |
	  parallel -j8 -N4095 -m --files samtools merge -u - |
	  parallel --xargs samtools merge -@${params.processors * 2} total.meth.bam {}

	### sort is not used
	### samtools sort -@ ${params.processors * 2} total.meth.bam
	samtools index -@ ${params.processors * 2} total.meth.bam
	echo "samtool index is done"

	tar -czf ${params.dsname}.guppy_fast5mod.combined.bam.tar.gz total.meth.bam*

	if [[ "${params.dataType}" == "human" ]] ; then
		echo "### For human, extract chr1-22, X and Y"
		## Ref: https://github.com/nanoporetech/medaka/issues/177
		for i in {1..22} X Y
		do
			fast5mod call total.meth.bam ${params.referenceGenome} \
				meth.chr\$i.tsv --meth cpg --quiet \
				--regions chr\$i &
		done
	elif [[ "${params.dataType}" == "ecoli" ]] ; then
		echo "### For ecoli, chr=${params.chrSet}"
		fast5mod call total.meth.bam ${params.referenceGenome} \
			meth.chr_${params.chrSet}.tsv \
			--meth cpg --quiet \
			--regions ${params.chrSet}
	fi
	wait

	cat  meth.chr*.tsv > ${params.dsname}.guppy.fast5mod_per_site.combine.tsv
	gzip ${params.dsname}.guppy.fast5mod_per_site.combine.tsv

	## Clean
	rm -f meth.chr*.tsv
	rm -f total.meth.bam*
	echo "### Guppy combine DONE"
	"""
}


// Combine Megalodon runs' all results together
process MgldnComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
	file x from megalodon_combine_in_ch

	output:
	file "${params.dsname}.megalodon.per_read.combine.bed.gz" into megalodon_combine_out_ch

	when:
	x.size() >= 1  && params.runCombine

	"""
	> ${params.dsname}.megalodon.per_read.combine.bed.gz
	cat ${x} > ${params.dsname}.megalodon.per_read.combine.bed.gz
	echo "### Megalodon combine DONE"
	"""
}


// Combine Nanopolish runs' all results together
process NplshComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
	file x from nanopolish_combine_in_ch

	output:
	file "${params.dsname}.nanopolish.per_read.combine.tsv.gz" into nanopolish_combine_out_ch

	when:
	x.size() >= 1 && params.runCombine

	"""
	> ${params.dsname}.nanopolish.per_read.combine.tsv.gz
	cat ${x} > ${params.dsname}.nanopolish.per_read.combine.tsv.gz
	echo "### Nanopolish combine DONE"
	"""
}


// Combine DeepMod runs' all results together
process DpmodComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy", pattern: "${params.dsname}.deepmod.*.combine.bed.gz"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepmod" , mode: "copy", pattern: "${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepmod" , mode: "copy", pattern: "${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz"
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepmod" , mode: "copy", pattern: "${params.dsname}.deepmod.all_batch.C.bed.tar.gz"

	input:
	file x from deepmod_combine_in_ch
	file deepmod_c_tar_file from Channel.fromPath(params.deepmod_ctar)
	each file("*") from ch_utils6

	output:
	file "${params.dsname}.deepmod.*.combine.bed.gz" into deepmod_combine_out_ch
	file "${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz" into deepmod_combine_sum_chrs_mod_ch
	file "${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz" into deepmod_combine_c_cluster_all_chrs_ch
	file "${params.dsname}.deepmod.all_batch.C.bed.tar.gz" into deepmod_combine_all_batch_c_ch

	when:
	x.size() >= 1 && params.runCombine

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
		indir/ C ${params.dsname}.deepmod ${params.chrSet}

	> ${params.dsname}.deepmod.C_per_site.combine.bed

	## Note: for ecoli data, no chr*, but N*
	for f in \$(ls -1 indir/${params.dsname}.deepmod.*.C.bed)
	do
	  cat \$f >> ${params.dsname}.deepmod.C_per_site.combine.bed
	done
	gzip ${params.dsname}.deepmod.C_per_site.combine.bed

	if [[ "${params.dataType}" == "human" ]] ; then
		## Only apply to human genome
		echo "### For human, apply cluster model of DeepMod"
		python utils/hm_cluster_predict.py \
			indir/${params.dsname}.deepmod \
			./C \
			\${DeepModProjectDir}/train_deepmod/${params.clusterDeepModModel} || true

		> ${params.dsname}.deepmod.C_clusterCpG_per_site.combine.bed
		for f in \$(ls -1 indir/${params.dsname}.deepmod_clusterCpG.*.C.bed)
		do
		  cat \$f >> ${params.dsname}.deepmod.C_clusterCpG_per_site.combine.bed
		done

		gzip ${params.dsname}.deepmod.C_clusterCpG_per_site.combine.bed
		tar -czf ${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz \
			indir/${params.dsname}.deepmod_clusterCpG.chr*.C.bed
	else
		touch ${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz
	fi

	tar -czf ${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz \
		indir/${params.dsname}.deepmod.*.C.bed
	tar -czf ${params.dsname}.deepmod.all_batch.C.bed.tar.gz \
		indir/batch_*_num/mod_pos.*.C.bed

	echo "### DeepMod combine DONE"
	"""
}


deepsignal_combine_out_ch.concat(tombo_combine_out_ch,megalodon_combine_out_ch, \
	nanopolish_combine_out_ch,deepmod_combine_out_ch.flatten(), guppy_fast5mod_combine_out_ch, guppy_gcf52ref_combine_out_ch)
	.toSortedList()
	.into { readlevel_unify_in; ; sitelevel_unify_in }


// Read level unified output, and get METEORE output
process METEORE {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy", pattern: "${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz"
	publishDir "${params.outputDir}/nanome-analysis-${params.dsname}" , mode: "copy", pattern: "Read_Level-${params.dsname}/${params.dsname}_*-METEORE-perRead-score.tsv.gz"

	input:
	file fileList from readlevel_unify_in
	each file("*") from ch_src3
	each file("*") from ch_utils8

	output:
	file "${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz" into meteore_combine_out_ch
	file "Read_Level-${params.dsname}/${params.dsname}_*-METEORE-perRead-score.tsv.gz" into unify_read_level_out_ch

	when:
	fileList.size() >= 1

	"""
	# Show all files
	flist=(\$(ls *.combine.{tsv,bed}.gz))
	echo \${flist[@]}

	nanopolishFile=\$(find . -maxdepth 1 -name '*.nanopolish.per_read.combine.*.gz')
	megalodonFile=\$(find . -maxdepth 1 -name '*.megalodon.per_read.combine.*.gz')
	deepsignalFile=\$(find . -maxdepth 1 -name '*.deepsignal.per_read.combine.*.gz')
	guppyFile=\$(find . -maxdepth 1 -name '*.guppy.*per_read.combine.*.gz')
	tomboFile=\$(find . -maxdepth 1 -name '*.tombo.per_read.combine.*.gz')

	tss_more_options=""
	if [[ "${params.dataType}" == "ecoli" ]] ; then
		tss_more_options="--chrs ${params.chrSet}"
	fi

	## Read level unify
	PYTHONPATH=src python src/nanocompare/tss_eval.py \
		--calls \
			Nanopolish:\${nanopolishFile} \
			Megalodon:\${megalodonFile} \
			DeepSignal:\${deepsignalFile} \
			Guppy.gcf52ref:\${guppyFile} \
			Tombo:\${tomboFile} \
		--runid Read_Level-${params.dsname} \
		--dsname ${params.dsname} --output-unified-format \
		--processors 8	-o . \${tss_more_options}

	## METEORE outputs by combining other tools
	nanopolishFileName=\$(find Read_Level-${params.dsname} -name "${params.dsname}_Nanopolish-METEORE-perRead-score.tsv.gz")
	deepsignalFileName=\$(find Read_Level-${params.dsname} -name "${params.dsname}_DeepSignal-METEORE-perRead-score.tsv.gz")
	megalodonFileName=\$(find Read_Level-${params.dsname} -name "${params.dsname}_Megalodon-METEORE-perRead-score.tsv.gz")

	outFileName=${params.dsname}_Megalodon_DeepSignal_combine.model_content.tsv
	> \$outFileName
	printf '%s\t%s\n' deepsignal \${deepsignalFileName} >> \$outFileName
	printf '%s\t%s\n' megalodon \${megalodonFileName} >> \$outFileName

	outFileName=${params.dsname}_Nanopolish_Megalodon_combine.model_content.tsv
	> \$outFileName
	printf '%s\t%s\n' megalodon \${megalodonFileName} >> \$outFileName
	printf '%s\t%s\n' nanopolish \${nanopolishFileName} >> \$outFileName

	wget https://github.com/comprna/METEORE/archive/refs/tags/v1.0.0.tar.gz --no-verbose
	tar -xzf v1.0.0.tar.gz
	METEORE_Dir="METEORE-1.0.0"

	## Degrade sk-learn for METEORE program
	pip install -U scikit-learn==0.21.3

	combineScript=utils/combination_model_prediction.py

	modelContentFileName=\$(find . -name "${params.dsname}_Megalodon_DeepSignal_combine.model_content.tsv")
	# Use the optimized model
	python \${combineScript} \
		-i \${modelContentFileName} -m optimized -b \${METEORE_Dir} \
		-o ${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz
	# Use the default model
	python \${combineScript} \
		-i \${modelContentFileName} -m default -b \${METEORE_Dir} \
		-o ${params.dsname}.meteore.megalodon_deepsignal_default_rf_model_per_read.combine.tsv.gz

	echo "### METEORE post combine DONE"
	"""
}


meteore_combine_out_ch.concat(sitelevel_unify_in.flatten())
	.toSortedList()
	.into { sitelevel_unify_in1; readlevel_in_ch; sitelevel_in_ch}


process SiteLevelUnify {
	publishDir "${params.outputDir}/nanome-analysis-${params.dsname}" , mode: "copy", pattern: "Site_Level-${params.dsname}/*.tss.*.cov1.bed.gz"

	input:
	file fileList from sitelevel_unify_in1
	each file("*") from ch_src4

	output:
	file "Site_Level-${params.dsname}/*.tss.*.cov1.bed.gz" into site_unify_out_ch

	when:
	fileList.size() >= 1

	"""
	# Show all files
	flist=(\$(ls *.combine.{tsv,bed}.gz))
	echo \${flist[@]}

	nanopolishFile=\$(find . -maxdepth 1 -name '*.nanopolish.per_read.combine.*.gz')
	megalodonFile=\$(find . -maxdepth 1 -name '*.megalodon.per_read.combine.*.gz')
	deepsignalFile=\$(find . -maxdepth 1 -name '*.deepsignal.per_read.combine.*.gz')
	guppyFile=\$(find . -maxdepth 1 -name '*.guppy.*per_site.combine.*.gz')
	tomboFile=\$(find . -maxdepth 1 -name '*.tombo.per_read.combine.*.gz')
	deepmodFile=\$(find . -maxdepth 1 -name '*.deepmod.C_clusterCpG_per_site.combine.*.gz')
	meteoreFile=\$(find . -maxdepth 1 -name '*.meteore.megalodon_deepsignal_optimized_model_per_read.combine.*.gz')

	## Site level unify
	PYTHONPATH=src python src/nanocompare/tss_eval.py \
		--calls \
			Nanopolish:\${nanopolishFile} \
			Megalodon:\${megalodonFile} \
			DeepSignal:\${deepsignalFile} \
			Guppy:\${guppyFile} \
			Tombo:\${tomboFile} \
			METEORE:\${meteoreFile} \
			DeepMod.Cluster:\${deepmodFile} \
		--runid Site_Level-${params.dsname} \
		--dsname ${params.dsname} \
		--processors 8 -o .
	"""
}


// Read level evaluations
process ReadLevelPerf {
	publishDir "${params.outputDir}/nanome-analysis-${params.dsname}" , mode: "copy"

	input:
	file fileList from readlevel_in_ch
	each file("*") from ch_src1

	output:
	file "MethPerf-*" into readlevel_out_ch

	when:
	params.eval && (fileList.size() >= 1)

	"""
	# Show all files
	flist=(\$(ls *.combine.{tsv,bed}.gz))
	echo \${flist[@]}

	nanopolishFile=\$(find . -maxdepth 1 -name '*.nanopolish.per_read.combine.*.gz')
	megalodonFile=\$(find . -maxdepth 1 -name '*.megalodon.per_read.combine.*.gz')
	deepsignalFile=\$(find . -maxdepth 1 -name '*.deepsignal.per_read.combine.*.gz')
	guppyFile=\$(find . -maxdepth 1 -name '*.guppy.fast5mod_per_site.combine.*.gz')
	tomboFile=\$(find . -maxdepth 1 -name '*.tombo.per_read.combine.*.gz')
	meteoreFile=\$(find . -maxdepth 1 -name '*.meteore.megalodon_deepsignal_optimized_model_per_read.combine.*.gz')
	deepmodFile=\$(find . -maxdepth 1 -name '*.deepmod.C_per_site.combine.*.gz')

	## Read level evaluations
	PYTHONPATH=src python src/nanocompare/read_level_eval.py \
		--calls \
				Nanopolish:\${nanopolishFile} \
				Megalodon:\${megalodonFile} \
				DeepSignal:\${deepsignalFile} \
				Guppy:\${guppyFile} \
				Tombo:\${tomboFile} \
				METEORE:\${meteoreFile} \
				DeepMod.C:\${deepmodFile} \
		--bgtruth "${params.bgtruthWithEncode}" \
		--runid MethPerf-${params.runid} \
		--dsname ${params.dsname} \
		--min-bgtruth-cov ${params.bgtruth_cov} \
		--processors ${params.processors * 2} \
		--report-joined -o . ## --distribution

	echo "### Read level analysis DONE"
	"""
}


// Site level correlation analysis
process SiteLevelCorr {
	publishDir "${params.outputDir}/nanome-analysis-${params.dsname}" , mode: "copy"

	input:
	file perfDir from readlevel_out_ch
	file fileList from sitelevel_in_ch
	each file("*") from ch_src2

	output:
	file "MethCorr-*" into sitelevel_out_ch

	when:
	params.eval && (fileList.size() >= 1)

	"""
	# Show all files
	flist=(\$(ls *.combine.{tsv,bed}.gz))
	echo \${flist[@]}

	nanopolishFile=\$(find . -maxdepth 1 -name '*.nanopolish.per_read.combine.*.gz')
	megalodonFile=\$(find . -maxdepth 1 -name '*.megalodon.per_read.combine.*.gz')
	deepsignalFile=\$(find . -maxdepth 1 -name '*.deepsignal.per_read.combine.*.gz')
	guppyFile=\$(find . -maxdepth 1 -name '*.guppy.fast5mod_per_site.combine.*.gz')
	tomboFile=\$(find . -maxdepth 1 -name '*.tombo.per_read.combine.*.gz')
	meteoreFile=\$(find . -maxdepth 1 -name '*.meteore.megalodon_deepsignal_optimized_model_per_read.combine.*.gz')
	deepmodFile=\$(find . -maxdepth 1 -name '*.deepmod.C_clusterCpG_per_site.combine.*.gz')

	## Site level evaluations
	PYTHONPATH=src  python src/nanocompare/site_level_eval.py \
		--calls \
				Nanopolish:\${nanopolishFile} \
				Megalodon:\${megalodonFile} \
				DeepSignal:\${deepsignalFile} \
				Tombo:\${tomboFile} \
				Guppy:\${guppyFile} \
				METEORE:\${meteoreFile} \
				DeepMod.Cluster:\${deepmodFile} \
		--bgtruth "${params.bgtruthWithEncode}" \
		--runid MethCorr-${params.runid} \
		--dsname ${params.dsname} \
		--min-bgtruth-cov ${params.bgtruth_cov} \
		--toolcov-cutoff ${params.tool_cov} \
		--beddir ${perfDir} \
		--processors ${params.processors * 2} \
		-o . --gen-venn ## --summary-coverage

	echo "### Site level analysis DONE"
	"""
}
