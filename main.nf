#!/usr/bin/env nextflow

def helpMessage() {
	log.info"""
	Usage:
	The typical command for running the pipeline is as follows:
	nextflow run https://github.com/liuyangzzu/nanome.git -profile singularity,hpc --dsname TestData --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt

	Mandatory arguments:
	  --dsname		Dataset name
	  --input		Input path for raw fast5 folders/tar/tar.gz files

	Options:
	  --processors		Processors used for each task
	  --outputDir		Output dir, default is 'outputs'
	  --dataType		Data type, default is 'human', can be also 'ecoli'
	  --chrSet		Chromosomes used in analysis, default is chr1-22, X and Y, seperated by comma. For E. coli data, it is set to 'NC_000913.3'

	  --cleanCache		True if clean work dir after complete
	  --commandType		Command used for tools, default is 'gpu', can be also 'cpu'

	  --queueName		SLURM job submission queue name for cluster running, default is 'gpu'
	  --qosName		SLURM job submission qos name for cluster running, default is 'inference'
	  --gresGPUOptions	SLURM job submission GPU allocation options for cluster running, default is '--gres=gpu:v100:1'
	  --jobMaxTime		SLURM job submission time allocation options for cluster running, default is '05:00:00'
	  --jobMaxMem		SLURM job submission memory allocation options for cluster running, default is '64G'

	  --conda_name			Conda name used for pipeline
	  --docker_name			Docker name used for pipeline
	  --singularity_name		Singularity name used for pipeline
	  --singularity_cache_dir	Singularity cache dir

	Other options:
	  --guppyDir		Guppy installation dir

	-profile options:
	  Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

	  docker 	A generic configuration profile to be used with Docker, pulls software from Docker Hub: quay.io/liuyangzzu/nanome:v1.5
	  singulairy	A generic configuration profile to be used with Singularity, pulls software from: docker://quay.io/liuyangzzu/nanome:v1.5
	  conda		Please only use conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity. Create conda enviroment by 'conda env create -f environment.yml'
	  hpc		A generic configuration profile to be used on HPC cluster with SLURM job submission support.
	  google	A generic configuration profile to be used on Google Cloud platform with 'google-lifesciences' support.

	Contact to https://github.com/liuyangzzu/nanome or https://nanome.jax.org for bug report.
	""".stripIndent()
}

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Check mandatory params
assert params.dsname != null : "Missing --dsname option, for command help use --help"
assert params.input != null : "Missing --input option, for command help use --help"

if (params.dataType == 'human') {
	referenceGenome="reference_genome/hg38/hg38.fasta"
	chromSizesFile="reference_genome/hg38/hg38.chrom.sizes"
} else if (params.dataType == 'ecoli') {
	referenceGenome="reference_genome/ecoli/Ecoli_k12_mg1655.fasta"
	chromSizesFile="reference_genome/ecoli/Ecoli_k12_mg1655.fasta.genome.sizes"
} else {
	println "Param dataType=${params.dataType} is not support"
	exit 1
}


log.info """\
NANOME - NF PIPELINE (v$workflow.manifest.version)
by Li Lab at The Jackson Laboratory
https://nanome.jax.org
=================================
dsname			:${params.dsname}
input			:${params.input}
output			:${params.outputDir}
work			:${workflow.workDir}
reference_genome	:${referenceGenome}
runBasecall		:${params.runBasecall}
runMethcall		:${params.runMethcall}
=================================
"""
.stripIndent()

projectDir = workflow.projectDir
ch_utils = Channel.fromPath("${projectDir}/utils",  type: 'dir', followLinks: false)
ch_src   = Channel.fromPath("${projectDir}/src",  type: 'dir', followLinks: false)


workflow.onComplete {
	if (workflow.success && params.cleanCache) {
		def workDir = new File("${workflow.workDir}")
		println "rm -rf ${workflow.workDir}".execute().text
	}
}

// Channel for utils/ and src/ folders
ch_utils.into{ch_utils1; ch_utils2; ch_utils3; ch_utils4; ch_utils5; ch_utils6; ch_utils7; ch_utils8; ch_utils9}
ch_src.into{ch_src1; ch_src2; ch_src3; ch_src4; ch_src5}

// Collect all folders of fast5 files, and send into Channels for pipelines
if (params.input.endsWith(".filelist.txt")) {
	// list of files in filelist.txt
	Channel.fromPath( params.input, checkIfExists: true )
		.splitCsv(header: false)
		.map { file(it[0]) }
		.set{ fast5_tar_ch }
} else if(params.input.endsWith("/*")) {
	// match all files in the folder
	Channel.fromPath(params.input, type: 'any').set{fast5_tar_ch}
} else {
	// For single file
	Channel.fromPath( params.input, checkIfExists: true ).set{fast5_tar_ch}
}

if (params.bgTruth) {
	Channel
    	.from(params.bgTruth.split(';') )
		.map { file(it) }
		.collect().into{in_bg_ch1; in_bg_ch2}
} else {
	in_bg_ch1 = Channel.empty()
	in_bg_ch2 = Channel.empty()
}


// Check all tools work well
process EnvCheck {
	tag 'EnvCheck'
	errorStrategy 'terminate'

	input:
	path reference_genome_tar 	from 	Channel.fromPath(params.reference_genome_tar, type: 'any')

	output:
	path "reference_genome" into reference_genome_ch

	"""
	set -x
	date; hostname; pwd

	which python
	python --version

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

	## Get dir for reference_genome
	if [[ "${reference_genome_tar}" == *.tar.gz ]] ; then
		tar -xzf ${reference_genome_tar}
	else
		if [[ "${reference_genome_tar}" != "reference_genome" ]] ; then
			mv ${reference_genome_tar} reference_genome
		fi
	fi
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

	input:
	path fast5_tar from 	fast5_tar_ch
	each path("*") from 	ch_utils5

	output:
	path "${fast5_tar.baseName}.untar" optional true into untar_out_ch

	"""
	infn="${fast5_tar}"
	mkdir -p untarTempDir
	if [ "\${infn##*.}" == "tar" ]; then
		### deal with tar
		tar -xf \${infn} -C untarTempDir
	elif [ "\${infn##*.}" == "gz" ]; then
		### deal with tar.gz
		tar -xzf \${infn} -C untarTempDir
	elif [[ -d ${fast5_tar} ]]; then
		## Copy files, do not change original files such as old analyses data
		cp -rf ${fast5_tar}/* untarTempDir/  || true # failed means nothing in this folder
	else
		echo "### Untar error for input=${fast5_tar}"
	fi

	## move fast5 files in tree folders into a single folder
	mkdir -p ${fast5_tar.baseName}.untar
	find untarTempDir -name "*.fast5" -type f -exec mv {} ${fast5_tar.baseName}.untar/ \\;

	## Clean old analyses in input fast5 files
	if [[ "${params.cleanAnalyses}" = true ]] ; then
		echo "### Start cleaning old analysis"
		python -c 'import h5py; print(h5py.version.info)'
		python utils/clean_old_basecall_in_fast5.py \
			-i ${fast5_tar.baseName}.untar --is-indir \
			--processor \$(( numProcessor*8 ))
	fi

	## Clean unused files
	rm -rf untarTempDir

	totalFiles=\$( find ${fast5_tar.baseName}.untar -name "*.fast5" -type f | wc -l )
	echo "### Total fast5 input files:\${totalFiles}"
	if (( totalFiles==0 )); then
		echo "### no fast5 files at ${fast5_tar.baseName}.untar, skip this job"
		rm -rf ${fast5_tar.baseName}.untar
	fi
	echo "### Untar DONE"
	"""
}


// Untar output will be used by basecall, megalodon and guppy
untar_out_ch.into{ untar_out_ch1; untar_out_ch2; untar_out_ch3 }


// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
	tag "${fast5_dir.baseName}"

	// errorStrategy = {task.attempt == params.maxRetries ? 'terminate' :  'retry'}

	input:
	path fast5_dir 				from 	untar_out_ch1
	each path("*") 				from 	ch_utils1
	each path(reference_genome) from 	reference_genome_ch9

	output:
	path "${fast5_dir.baseName}.basecalled" into basecall_out_ch  // try to fix the christina proposed problems
	path "${fast5_dir.baseName}.basecalled/${fast5_dir.baseName}-sequencing_summary.txt" into qc_ch
	path "${fast5_dir.baseName}.basecalled.bam" into ont_cov_bam_ch

	when:
	params.runBasecall

	"""
	date; hostname; pwd
	which guppy_basecaller
	mkdir -p ${fast5_dir.baseName}.basecalled

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path "${fast5_dir.baseName}.basecalled" \
			--config ${params.GUPPY_BASECALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out \
			--verbose_logs
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path "${fast5_dir.baseName}.basecalled" \
			--config ${params.GUPPY_BASECALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out \
			--verbose_logs \
			-x auto
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi

	## After basecall, rename and publishe summary filenames
	mv ${fast5_dir.baseName}.basecalled/sequencing_summary.txt \
		${fast5_dir.baseName}.basecalled/${fast5_dir.baseName}-sequencing_summary.txt

	## After basecall, we process guppy results for ONT coverage analyses
	# align FASTQ files to reference genome, write sorted alignments to a BAM file
	minimap2 -a -z 600,200 -x map-ont ${referenceGenome} ${fast5_dir.baseName}.basecalled/*.fastq \
	    -t \$(( numProcessor*2 )) > ${fast5_dir.baseName}.basecalled.sam
    echo "Alignment done"

    # Convert the sam files to bam (a binary sam format) using samtools’ view command
    samtools view -u ${fast5_dir.baseName}.basecalled.sam \
        | samtools sort -@ \$(( numProcessor*2 )) -o ${fast5_dir.baseName}.basecalled.bam --output-fmt BAM
    echo "samtools done"

    # Clean
    rm -f ${fast5_dir.baseName}.basecalled.sam

	echo "### Basecalled by Guppy DONE"
	"""
}


// Collect and output QC results for basecall, and report ONT coverage
process QCExport {
	publishDir "${params.outputDir}/${params.dsname}-basecallings",
		mode: "copy", enabled: params.outputQC

	input:
	path flist 		from 	qc_ch.collect()
	path bamlist 	from 	ont_cov_bam_ch.collect()

	output:
	path "${params.dsname}_basecall_report.html" into qc_out_ch
	path "${params.dsname}_QCReport" into qc_report_out_ch

	"""
	## Combine all sequencing summary files
	touch ${params.dsname}_combine_sequencing_summary.txt.gz
	fnlist=\$(find . -name '*-sequencing_summary.txt')
	firstFile=true
	for fn in \$fnlist; do
		if \$firstFile ; then
			awk '(NR>=1)' \$fn | \
				gzip >> ${params.dsname}_combine_sequencing_summary.txt.gz
			firstFile=false
		else
			awk '(NR>1)' \$fn | \
				gzip >> ${params.dsname}_combine_sequencing_summary.txt.gz
		fi
	done

	## QC report generation
	if ! command -v NanoComp &> /dev/null
	then
		echo "NanoComp could not be found, install it..."
		conda install -c bioconda nanocomp
	else
		NanoComp -v
	fi
	NanoComp --summary ${params.dsname}_combine_sequencing_summary.txt.gz  \
		--names ${params.dsname} --outdir ${params.dsname}_QCReport -t \$(( numProcessor )) \
		--verbose  --raw  -f pdf -p ${params.dsname}_

    ## Merge the bam file
	find . -maxdepth 1 -name "*.bam" | \
	    parallel -j8 -N4095 -m --files samtools merge -u - | \
	    parallel --xargs samtools merge -@ \$(( numProcessor*2 )) \
	    ${params.dsname}_merged.bam {}";" rm {}

    samtools index ${params.dsname}_merged.bam  -@ \$(( numProcessor*2 ))
    echo "Samtools merging done!"

    ## calculates the sequence coverage at each position/ Reporting genome coverage for all positions in BEDGRAPH format.
    bedtools genomecov -ibam ${params.dsname}_merged.bam -bg -strand + |
        awk '\$4 = \$4 FS "+"' |
        gzip > ${params.dsname}.coverage.positivestrand.bed.gz

    bedtools genomecov -ibam ${params.dsname}_merged.bam -bg -strand - |
        awk '\$4 = \$4 FS "-"' |
        gzip > ${params.dsname}.coverage.negativestrand.bed.gz

    cat ${params.dsname}.coverage.positivestrand.bed.gz > ${params.dsname}_ONT_coverage_combine.bed.gz
	cat ${params.dsname}.coverage.negativestrand.bed.gz >> ${params.dsname}_ONT_coverage_combine.bed.gz

	mv ${params.dsname}_ONT_coverage_combine.bed.gz ${params.dsname}_QCReport/
	mv ${params.dsname}_combine_sequencing_summary.txt.gz ${params.dsname}_QCReport/
	mv ${params.dsname}_QCReport/${params.dsname}_NanoComp-report.html ${params.dsname}_basecall_report.html

	## Clean
	rm -f ${params.dsname}.coverage.positivestrand.bed.gz ${params.dsname}.coverage.negativestrand.bed.gz
	rm -f ${params.dsname}_merged.bam*

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

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/guppy",
		mode: "copy",
		pattern: "outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz",
		enabled: params.outputIntermediate
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/guppy",
		mode: "copy",
		pattern: "batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz",
		enabled: params.outputIntermediate

	input:
	path fast5_dir 					from 	untar_out_ch2
	each path(reference_genome) 	from 	reference_genome_ch1
	each path("*") 					from 	ch_utils4

	output:
	path "outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz" into guppy_methcall_gz_out_ch
	path "batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam*" into guppy_methcall_out_ch
	path "batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz" into guppy_gcf52ref_out_ch

	when:
	params.runMethcall && params.runGuppy

	"""
	date; hostname; pwd

	mkdir -p ${fast5_dir.baseName}.methcalled

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path ${fast5_dir.baseName}.methcalled \
			--config ${params.GUPPY_METHCALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out \
			--verbose_logs
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path ${fast5_dir.baseName}.methcalled \
			--config ${params.GUPPY_METHCALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out \
			--verbose_logs \
			--device auto
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi
	echo "### Guppy methylation calling DONE"

	## Extract guppy methylation-callings
	## gcf52ref ways
	minimap2 -t \$(( numProcessor*2 )) -a -x map-ont ${referenceGenome} \
		${fast5_dir.baseName}.methcalled/*.fastq | \
		samtools sort -@ \$(( numProcessor*2 )) \
		-T tmp -o gcf52ref.batch.${fast5_dir.baseName}.bam

	samtools index -@ \$(( numProcessor*2 )) \
		gcf52ref.batch.${fast5_dir.baseName}.bam
	echo "### gcf52ref minimap2 alignment is done!"

	## Modified version, support dir input, not all fast5 files (too long arguments)
	python utils/extract_methylation_fast5_support_dir.py \
		-p \$(( numProcessor*2 )) ${fast5_dir.baseName}.methcalled/workspace
	echo "### gcf52ref extract to db DONE"

	## gcf52ref files preparation
	### git clone https://github.com/kpalin/gcf52ref.git
	tar -xzf utils/gcf52ref.tar.gz -C .
	patch gcf52ref/scripts/extract_methylation_from_rocks.py < utils/gcf52ref.patch

	python gcf52ref/scripts/extract_methylation_from_rocks.py \
		-d base_mods.rocksdb \
		-a gcf52ref.batch.${fast5_dir.baseName}.bam \
		-r ${referenceGenome} \
		-o tmp.batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv
	echo "### gcf52ref extract to tsv DONE"

	tail -n +2 tmp.batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv | gzip > \
		batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz
	echo "### gcf52ref DONE"

	## fast5mod ways
	FAST5PATH=${fast5_dir.baseName}.methcalled/workspace
	OUTBAM=batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam

	fast5mod guppy2sam \${FAST5PATH} --reference ${referenceGenome} \
		--workers 74 --recursive --quiet \
		| samtools sort -@ \$(( numProcessor*2 )) | \
		samtools view -b -@ \$(( numProcessor*2 )) > \${OUTBAM}

	samtools index -@ \$(( numProcessor*2 ))  \${OUTBAM}

	tar -czf outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz \
		batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam*
	echo "### fast5mod DONE"
	echo "### Guppy fast5mod and gcf52ref DONE"
	"""
}


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${fast5_dir.baseName}"

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/megalodon",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path fast5_dir 					from untar_out_ch3
	each path(reference_genome) 	from reference_genome_ch2
	each path(megalodonModelTar) 	from Channel.fromPath(params.megalodon_model_tar)

	output:
	path "batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt.gz" into megalodon_out_ch

	when:
	params.runMethcall && params.runMegalodon

	"""
	date; hostname; pwd

	## Get megalodon model dir
	tar -xzf ${megalodonModelTar}

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		## Ref: https://github.com/nanoporetech/megalodon
		megalodon \
			${fast5_dir} \
			--overwrite \
			--outputs per_read_mods mods per_read_refs \
			--guppy-server-path guppy_basecall_server \
			--guppy-config ${params.MEGALODON_MODEL_FOR_GUPPY_CONFIG} \
			--guppy-params "-d ./megalodon_model/ --num_callers \$(( numProcessor )) --ipc_threads 80" \
			--guppy-timeout ${params.GUPPY_TIMEOUT} \
			--samtools-executable ${params.SAMTOOLS_PATH} \
			--sort-mappings \
			--mappings-format bam \
			--reference ${referenceGenome} \
			--mod-motif m CG 0 \
			--mod-output-formats bedmethyl wiggle \
			--write-mods-text \
			--write-mod-log-probs \
			--processes \$(( numProcessor*2 ))
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		## Ref: https://github.com/nanoporetech/megalodon
		megalodon \
			${fast5_dir} \
			--overwrite \
			--outputs per_read_mods mods per_read_refs \
			--guppy-server-path guppy_basecall_server \
			--guppy-config ${params.MEGALODON_MODEL_FOR_GUPPY_CONFIG} \
			--guppy-params "-d ./megalodon_model/ --num_callers \$(( numProcessor )) --ipc_threads 80" \
			--guppy-timeout ${params.GUPPY_TIMEOUT} \
			--samtools-executable ${params.SAMTOOLS_PATH} \
			--sort-mappings \
			--mappings-format bam \
			--reference ${referenceGenome} \
			--mod-motif m CG 0 \
			--mod-output-formats bedmethyl wiggle \
			--write-mods-text \
			--write-mod-log-probs \
			--processes \$(( numProcessor*2 )) \
			--devices 0
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi

	### mv megalodon_results/per_read_modified_base_calls.txt batch_${fast5_dir.baseName}.per_read_modified_base_calls.txt
	tail -n +2 megalodon_results/per_read_modified_base_calls.txt | gzip > \
		batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt.gz
	echo "### Megalodon DONE"
	"""
}


// Resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
	tag "${basecallIndir.baseName}"

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/resquiggle",
		mode: "copy",
		pattern: "${basecallIndir.baseName}.resquiggle.run.log",
		enabled: params.outputIntermediate

	input:
	path 	basecallIndir 			from resquiggle_in_ch
	each 	path(reference_genome) 	from reference_genome_ch3

	output:
	path "${basecallIndir.baseName}.resquiggle" into resquiggle_out_ch
	path "${basecallIndir.baseName}.resquiggle.run.log" into resquiggle_logs

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
		--processes \$(( numProcessor*2 )) \
		--corrected-group ${params.resquiggleCorrectedGroup} \
		--basecall-group ${params.BasecallGroupName} \
		--overwrite \
		${basecallIndir.baseName}.resquiggle/workspace \
		${referenceGenome} &> \
		${basecallIndir.baseName}.resquiggle.run.log
	echo "### Resquiggle DONE"
	"""
}


// duplicate resquiggle results for DeepSignal and Tombo
resquiggle_out_ch.into { deepsignal_in_ch; tombo_in_ch }


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${indir.baseName}"

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepsignal",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path indir 						from deepsignal_in_ch
	each path(deepsignal_model_tar) from Channel.fromPath(params.deepsignal_model_tar)
	each path(reference_genome) 	from reference_genome_ch4

	output:
	path "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv.gz" into deepsignal_out_ch

	when:
	params.runMethcall && params.runDeepSignal && !params.filterGPUTaskRuns

	"""
	tar -xzf ${deepsignal_model_tar}

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path "./${params.DEEPSIGNAL_MODEL}" \
			--result_file "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv" \
			--reference_path ${referenceGenome} \
			--corrected_group ${params.resquiggleCorrectedGroup} \
			--nproc \$(( numProcessor * ${params.deepLearningProcessorTimes}  )) \
			--is_gpu no
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path "./${params.DEEPSIGNAL_MODEL}" \
			--result_file "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv" \
			--reference_path ${referenceGenome} \
			--corrected_group ${params.resquiggleCorrectedGroup} \
			--nproc \$(( numProcessor * ${params.deepLearningProcessorTimes}  )) \
			--is_gpu yes
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi

	gzip batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv
	echo "### DeepSignal methylation DONE"
	"""
}


// Tombo runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Tombo {
	tag "${resquiggleDir.baseName}"

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/tombo",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	each path(reference_genome) from reference_genome_ch5
	each path("*") 				from ch_utils2
	path resquiggleDir	 		from tombo_in_ch

	output:
	path "batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats.bed.gz" into tombo_out_ch
	path "${resquiggleDir.baseName}.tombo.run.log" into tombo_log_ch

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
		--processes \$(( numProcessor )) \
		--corrected-group ${params.resquiggleCorrectedGroup} \
		--multiprocess-region-size 1000 &> \
		${resquiggleDir.baseName}.tombo.run.log

	retry=1
	while grep -q "BrokenPipeError:" ${resquiggleDir.baseName}.tombo.run.log
	do
		echo "### Found error 32, repeat tombo running again!!!"
		tombo detect_modifications alternative_model \
			--fast5-basedirs ${resquiggleDir}/workspace \
			--dna --standard-log-likelihood-ratio \
			--statistics-file-basename batch_${resquiggleDir.baseName} \
			--per-read-statistics-basename batch_${resquiggleDir.baseName} \
			--alternate-bases CpG \
			--processes \$(( numProcessor )) \
			--corrected-group ${params.resquiggleCorrectedGroup} \
			--multiprocess-region-size 1000 &> \
			${resquiggleDir.baseName}.tombo.run.log
		retry=\$(( retry+1 ))
		if (( retry >= 8 )); then
			break
		fi
	done

	if grep -q "BrokenPipeError: \\[Errno 32\\] Broken pipe" ${resquiggleDir.baseName}.tombo.run.log; then
		## Grep the broken pipeline bug for Tombo
		echo "### Tombo bug occur, max retry reached at \${retry} times, return error"
		exit 32
	else
		## Tombo was ok
		echo "### Tombo log passed"
	fi

	python utils/tombo_extract_per_read_stats.py \
		${chromSizesFile} \
		"batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats" \
		"batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats.bed"

	gzip batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats.bed
	echo "### Tombo methylation calling DONE"
	"""
}


// DeepMod runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepMod {
	tag "${basecallDir.baseName}"

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepmod",
		mode: "copy", pattern: "batch_${basecallDir.baseName}_num.tar.gz",
		enabled: params.outputIntermediate

	input:
	each path(reference_genome) from reference_genome_ch6
	path basecallDir 			from deepmod_in_ch

	output:
	path "mod_output/batch_${basecallDir.baseName}_num" into deepmod_out_ch
	path "batch_${basecallDir.baseName}_num.tar.gz" into deepmod_gz_out_ch

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
			--Ref ${referenceGenome} \
			--Base C \
			--modfile \${DeepModProjectDir}/train_deepmod/${params.DEEPMOD_RNN_MODEL} \
			--FileID batch_${basecallDir.baseName}_num \
			--threads \$(( numProcessor*${params.deepLearningProcessorTimes} ))  ${params.DeepModMoveOptions}

	tar -czf batch_${basecallDir.baseName}_num.tar.gz mod_output/batch_${basecallDir.baseName}_num/
	echo "### DeepMod methylation DONE"
	"""
}


// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Nanopolish {
	tag "${basecallDir.baseName}"

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/nanopolish",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path basecallDir 			from nanopolish_in_ch
	each path(reference_genome) from reference_genome_ch7
	each path("*") 				from ch_utils3

	output:
	path "batch_${basecallDir.baseName}.nanopolish.methylation_calls.tsv.gz" into nanopolish_out_ch

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

	minimap2 -t \$(( numProcessor*2 )) -a -x map-ont ${referenceGenome} \${fastqNoDupFile} | \
		samtools sort -@ \$(( numProcessor*2 )) -T tmp -o \${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ \$(( numProcessor*2 ))  \${bamFileName}
	echo "### samtools finished"
	echo "### Alignment step DONE"

	nanopolish call-methylation -t \$(( numProcessor*2 )) -r \${fastqNoDupFile} \
		-b \${bamFileName} -g ${referenceGenome} > tmp.tsv

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
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		enabled: params.outputRaw

	input:
	path x from deepsignal_combine_in_ch

	output:
	path "${params.dsname}.deepsignal.per_read.combine.tsv.gz" into deepsignal_combine_out_ch

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
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		enabled: params.outputRaw

	input:
	path x from tombo_combine_in_ch // list of tombo bed files

	output:
	path "${params.dsname}.tombo.per_read.combine.bed.gz" into tombo_combine_out_ch

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
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}.guppy.*.combine.tsv.gz",
		enabled: params.outputRaw

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/guppy",
		mode: "copy", pattern: "${params.dsname}.guppy_fast5mod.combined.bam.tar.gz",
		enabled: params.outputIntermediate

	input:
	path x 					from guppy_combine_in_ch
	path y 					from guppy_gcf52ref_out_ch.collect()
	path reference_genome 	from reference_genome_ch8

	output:
	path "${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz" into guppy_fast5mod_combine_out_ch
	path "${params.dsname}.guppy.gcf52ref_per_read.combine.tsv.gz" into guppy_gcf52ref_combine_out_ch
	path "${params.dsname}.guppy_fast5mod.combined.bam.tar.gz" into guppy_combine_raw_out_ch

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
	  parallel --xargs samtools merge -@\$(( numProcessor*2 )) total.meth.bam {}

	### sort is not used
	### samtools sort -@ \$(( numProcessor*2 )) total.meth.bam
	samtools index -@ \$(( numProcessor*2 )) total.meth.bam
	echo "samtool index is done"

	tar -czf ${params.dsname}.guppy_fast5mod.combined.bam.tar.gz total.meth.bam*

	touch ${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz
	if [[ "${params.dataType}" == "human" ]] ; then
		echo "### For human, extract chr1-22, X and Y"
		## Ref: https://github.com/nanoporetech/medaka/issues/177
		for i in {1..22} X Y
		do
			fast5mod call total.meth.bam ${referenceGenome} \
				meth.chr_\$i.tsv \
				--meth cpg --quiet \
				--regions chr\$i &
		done
	elif [[ "${params.dataType}" == "ecoli" ]] ; then
		echo "### For ecoli, chr=${params.chrSet}"
		fast5mod call total.meth.bam ${referenceGenome} \
			meth.chr_${params.chrSet}.tsv \
			--meth cpg --quiet \
			--regions ${params.chrSet}
		cat  meth.chr_${params.chrSet}.tsv | gzip >> ${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz
	fi

	if [[ "${params.dataType}" == "human" ]] ; then
		wait
		echo "### For human, extract chr1-22, X and Y"
		for i in {1..22} X Y
		do
			cat  meth.chr_\$i.tsv | \
				gzip >> ${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz
		done
	fi

	## Clean
	rm -f meth.chr*.tsv
	rm -f total.meth.bam*
	echo "### Guppy combine DONE"
	"""
}


// Combine Megalodon runs' all results together
process MgldnComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		enabled: params.outputRaw

	input:
	path x from megalodon_combine_in_ch

	output:
	path "${params.dsname}.megalodon.per_read.combine.bed.gz" into megalodon_combine_out_ch

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
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		enabled: params.outputRaw

	input:
	path x from nanopolish_combine_in_ch

	output:
	path "${params.dsname}.nanopolish.per_read.combine.tsv.gz" into nanopolish_combine_out_ch

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
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}.deepmod.*.combine.bed.gz",
		enabled: params.outputRaw
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepmod",
		mode: "copy", pattern: "${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz",
		enabled: params.outputIntermediate
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepmod",
		mode: "copy", pattern: "${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz",
		enabled: params.outputIntermediate
	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepmod",
		mode: "copy", pattern: "${params.dsname}.deepmod.all_batch.C.bed.tar.gz",
		enabled: params.outputIntermediate

	input:
	path x 					from deepmod_combine_in_ch
	path deepmod_c_tar_file from Channel.fromPath(params.deepmod_ctar)
	each path("*") 			from ch_utils6

	output:
	path "${params.dsname}.deepmod.*.combine.bed.gz" into deepmod_combine_out_ch
	path "${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz" into deepmod_combine_sum_chrs_mod_ch
	path "${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz" into deepmod_combine_c_cluster_all_chrs_ch
	path "${params.dsname}.deepmod.all_batch.C.bed.tar.gz" into deepmod_combine_all_batch_c_ch

	when:
	x.size() >= 1 && params.runCombine

	"""
	wget ${params.DeepModGithub} --no-verbose
	tar -xzf v0.1.3.tar.gz
	DeepModProjectDir="DeepMod-0.1.3"

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

	## Note: for ecoli data, no pattern for chr*, but N*
	for f in \$(ls -1 indir/${params.dsname}.deepmod.*.C.bed)
	do
	  cat \$f >> ${params.dsname}.deepmod.C_per_site.combine.bed
	done
	gzip ${params.dsname}.deepmod.C_per_site.combine.bed

	if [[ "${params.dataType}" == "human" ]] ; then
		## Only apply to human genome
		echo "### For human, apply cluster model of DeepMod"

		## Get dir for deepmod cluster-model inputs
		if [[ "${deepmod_c_tar_file}" == *.tar.gz ]] ; then
			tar -xzf ${deepmod_c_tar_file}
		else
			if [[ "${deepmod_c_tar_file}" != "C" ]] ; then
				mv ${deepmod_c_tar_file} C
			fi
		fi

		python utils/hm_cluster_predict.py \
			indir/${params.dsname}.deepmod \
			./C \
			\${DeepModProjectDir}/train_deepmod/${params.DEEPMOD_CLUSTER_MODEL} ## || true

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


deepsignal_combine_out_ch
	.concat(tombo_combine_out_ch,
			megalodon_combine_out_ch,
			nanopolish_combine_out_ch,
			guppy_gcf52ref_combine_out_ch)
	.toSortedList()
	.into { raw_results_five_tools1; raw_results_five_tools2 }


// Read level unified output, and get METEORE output
process METEORE {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz",
		enabled: params.outputRaw
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	input:
	path fileList 	from 	raw_results_five_tools1
	each path("*") 	from 	ch_src3
	each path("*") 	from 	ch_utils8

	output:
	path "${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz" into meteore_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz" into read_unify_out_ch

	when:
	fileList.size() >= 1 && !params.filterGPUTaskRuns

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
		--processors \$(( numProcessor*2 ))	\
		-o . \${tss_more_options}

	## METEORE outputs by combining other tools
	nanopolishFileName=\$(find Read_Level-${params.dsname} -name "${params.dsname}_Nanopolish-perRead-score.tsv.gz")
	deepsignalFileName=\$(find Read_Level-${params.dsname} -name "${params.dsname}_DeepSignal-perRead-score.tsv.gz")
	megalodonFileName=\$(find Read_Level-${params.dsname} -name "${params.dsname}_Megalodon-perRead-score.tsv.gz")

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

	##cp ${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz \
	##	Read_Level-${params.dsname}/TestData_METEORE-perRead-score.tsv.gz

	zcat ${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz | \
		awk -F '\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$6,\$5}' |
		gzip > Read_Level-${params.dsname}/${params.dsname}_METEORE-perRead-score.tsv.gz

	echo "### METEORE post combine DONE"
	"""
}


meteore_combine_out_ch
	.concat(raw_results_five_tools2.flatten(),
			deepmod_combine_out_ch.flatten(),
			guppy_fast5mod_combine_out_ch)
	.toSortedList()
	.into { all_raw_results1; all_raw_results2; all_raw_results3; all_raw_results4}


process SiteLevelUnify {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	publishDir "${params.outputDir}",
		mode: "copy", pattern: "README.txt"

	input:
	path fileList 	from 	all_raw_results1
	each path("*") 	from 	ch_src4
	each path("*") 	from 	ch_utils9

	output:
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz" into site_unify_out_ch
	path "README.txt" into readme_out_ch

	when:
	fileList.size() >= 1 && !params.filterGPUTaskRuns

	"""
	# Show all files
	flist=(\$(ls *.combine.{tsv,bed}.gz))
	echo \${flist[@]}

	nanopolishFile=\$(find . -maxdepth 1 -name '*.nanopolish.per_read.combine.*.gz')
	megalodonFile=\$(find . -maxdepth 1 -name '*.megalodon.per_read.combine.*.gz')
	deepsignalFile=\$(find . -maxdepth 1 -name '*.deepsignal.per_read.combine.*.gz')
	guppyFile=\$(find . -maxdepth 1 -name '*.guppy.*per_site.combine.*.gz')
	tomboFile=\$(find . -maxdepth 1 -name '*.tombo.per_read.combine.*.gz')
	deepmodFile=\$(find . -maxdepth 1 -name '*.deepmod.C_per_site.combine.*.gz')
	meteoreFile=\$(find . -maxdepth 1 -name '*.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.*.gz')
	deepmodEncode="DeepMod.C"

	if [[ "${params.useDeepModCluster}" == true ]] ; then
		deepmodFile=\$(find . -maxdepth 1 -name '*.deepmod.C_clusterCpG_per_site.combine.*.gz')
		deepmodEncode="DeepMod.Cluster"
	fi

	tss_more_options=""
	if [[ "${params.dataType}" == "ecoli" ]] ; then
		tss_more_options="--chrs ${params.chrSet}"
		deepmodFile=\$(find . -maxdepth 1 -name '*.deepmod.C_per_site.combine.*.gz')
		deepmodEncode="DeepMod.C"
	fi

	## Site level unify
	PYTHONPATH=src python src/nanocompare/tss_eval.py \
		--calls \
			Nanopolish:\${nanopolishFile} \
			Megalodon:\${megalodonFile} \
			DeepSignal:\${deepsignalFile} \
			Guppy:\${guppyFile} \
			Tombo:\${tomboFile} \
			METEORE:\${meteoreFile} \
			\${deepmodEncode}:\${deepmodFile} \
		--runid Site_Level-${params.dsname} \
		--dsname ${params.dsname} \
		--processors \$(( numProcessor*2 )) \
		-o .  \${tss_more_options}

	## Sort site level results
	fnlist=\$(find Site_Level-${params.dsname} -name '*-perSite-cov1.bed.gz')
	for fn in \$fnlist ; do
		outfn=\${fn/cov1.bed.gz/cov1.sort.bed.gz}
		bedtools sort -i \${fn} | gzip > \${outfn}
		rm \${fn}
	done

	PYTHONIOENCODING=UTF-8 python utils/gen_readme.py \
		utils/readme.txt.template ${params.dsname} ${params.outputDir} \
		${workflow.projectDir} ${workflow.workDir} "${workflow.commandLine}" ${workflow.runName} "${workflow.start}" \
		> README.txt
	"""
}


qc_report_out_ch
	.concat(all_raw_results2.flatten(), site_unify_out_ch.flatten(), read_unify_out_ch.flatten())
	.toSortedList()
	.set {report_in_ch}


process ReportHtml {
	publishDir "${params.outputDir}", mode: "copy"

	input:
	path fileList 	from 	report_in_ch
	each path("*") 	from 	ch_src5

	output:
	path "report" 	into	report_out_ch

	when:
	fileList.size() >= 1

	"""
	## Generate running information tsv
	> running_information.tsv
	printf '%s\t%s\n' Title Information >> running_information.tsv
	printf '%s\t%s\n' dsname ${params.dsname} >> running_information.tsv
	printf '%s\t%s\n' projectDir ${workflow.projectDir} >> running_information.tsv
	printf '%s\t%s\n' workDir ${workflow.workDir} >> running_information.tsv
	printf '%s\t%s\n' commandLine "${workflow.commandLine}" >> running_information.tsv
	printf '%s\t%s\n' runName ${workflow.runName} >> running_information.tsv
	printf '%s\t%s\n' start ${workflow.start} >> running_information.tsv
	printf '%s\t%s\n' input "${params.input}" >> running_information.tsv
	printf '%s\t%s\n' outputDir ${params.outputDir} >> running_information.tsv

	## Get basecalling results from NanoComp
	basecallOutputFile=\$(find ${params.dsname}_QCReport/ -name "*NanoStats.txt" -type f)

	## Generate report dir and html utilities
	mkdir -p report
	cp src/nanocompare/report/style.css report/
	cp -rf src/nanocompare/report/icons report/
	cp -rf src/nanocompare/report/js report/

	PYTHONPATH=src  python src/nanocompare/report/gen_html_report.py \
		${params.dsname} \
		running_information.tsv \
		\${basecallOutputFile} \
		. \
		report

	echo "### report html DONE"
	"""
}


// Read level evaluations
process ReadLevelPerf {
	publishDir "${params.outputDir}/nanome-analysis-${params.dsname}", mode: "copy"

	input:
	path fileList 	from 	all_raw_results3
	each path("*") 	from 	ch_src1
	path bgFiles 	from 	in_bg_ch1

	output:
	path "MethPerf-*" into readlevel_out_ch

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
	meteoreFile=\$(find . -maxdepth 1 -name '*.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.*.gz')
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
		--bgtruth "${params.bgEncode}:${bgFiles.join(';')}" \
		--runid MethPerf-${params.runid} \
		--dsname ${params.dsname} \
		--min-bgtruth-cov ${params.bgtruth_cov} \
		--processors \$(( numProcessor*2 )) \
		--report-joined -o . ## --distribution

	echo "### Read level analysis DONE"
	"""
}


// Site level correlation analysis
process SiteLevelCorr {
	publishDir "${params.outputDir}/nanome-analysis-${params.dsname}", mode: "copy"

	input:
	path perfDir 	from 	readlevel_out_ch
	path fileList 	from 	all_raw_results4
	each path("*") 	from 	ch_src2
	path bgFiles 	from 	in_bg_ch2

	output:
	path "MethCorr-*" into sitelevel_out_ch

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
	meteoreFile=\$(find . -maxdepth 1 -name '*.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.*.gz')
	deepmodFile=\$(find . -maxdepth 1 -name '*.deepmod.C_clusterCpG_per_site.combine.*.gz')

	## Site level evaluations
	PYTHONPATH=src python src/nanocompare/site_level_eval.py \
		--calls \
				Nanopolish:\${nanopolishFile} \
				Megalodon:\${megalodonFile} \
				DeepSignal:\${deepsignalFile} \
				Tombo:\${tomboFile} \
				Guppy:\${guppyFile} \
				METEORE:\${meteoreFile} \
				DeepMod.Cluster:\${deepmodFile} \
		--bgtruth "${params.bgEncode}:${bgFiles.join(';')}" \
		--runid MethCorr-${params.runid} \
		--dsname ${params.dsname} \
		--min-bgtruth-cov ${params.bgtruth_cov} \
		--toolcov-cutoff ${params.tool_cov} \
		--beddir ${perfDir} \
		--processors \$(( numProcessor*2 )) \
		-o . --gen-venn ## --summary-coverage

	echo "### Site level analysis DONE"
	"""
}
