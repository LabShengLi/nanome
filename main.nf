#!/usr/bin/env nextflow
// @Author   : Yang Liu
// @FileName : main.nf
// @Software : NANOME project
// @Organization : JAX Li Lab
// @Website  : https://github.com/TheJacksonLaboratory/nanome

def helpMessage() { // print help message
	log.info"""
	NANOME - NF PIPELINE (v$workflow.manifest.version)
	by Li Lab at The Jackson Laboratory
	https://github.com/TheJacksonLaboratory/nanome
	=================================
	Usage:
	The typical command for running the pipeline is as follows:

	nextflow run TheJacksonLaboratory/nanome -profile ci,conda
	nextflow run TheJacksonLaboratory/nanome -profile ci,singularity

	Mandatory arguments:
	  --dsname		Dataset name
	  --input		Input path for raw fast5 folders/tar/tar.gz files

	General options:
	  --processors		Processors used for each task
	  --outputDir		Output dir, default is 'outputs'
	  --dataType		Data type, default is 'human', can be also 'ecoli'
	  --chrSet		Chromosomes used in analysis, default is true, means chr1-22, X and Y, seperated by comma. For E. coli data, it needs be set to 'NC_000913.3'

	  --cleanCache		If clean work dir after complete, default is true

	Running environment options:
	  --conda_name			Conda name used for pipeline, default is 'nanome'
	  --conda_base_dir		Conda base directory, default is '/opt/conda'
	  --docker_name			Docker name used for pipeline, default is 'liuyangzzu/nanome:latest'
	  --singularity_name		Singularity name used for pipeline, default is 'docker://liuyangzzu/nanome:latest'
	  --singularity_cache		Singularity cache dir, default is 'local_singularity_cache'

	Platform specific options:
	  --queueName		SLURM job submission queue name for cluster running, default is 'gpu'
	  --qosName		SLURM job submission qos name for cluster running, default is 'inference'
	  --gresGPUOptions	SLURM job submission GPU allocation options for cluster running, default is '--gres=gpu:v100:1'
	  --jobMaxTime		SLURM job submission time allocation options for cluster running, default is '05:00:00'
	  --jobMaxMem		SLURM job submission memory allocation options for cluster running, default is '32G'

	  --googleProjectName	Google Cloud project name for google-lifesciences task running

	Other options:
	  --guppyDir		Guppy installation local directory, used only for conda environment

	-profile options:
	  Use this parameter to choose a predefined configuration profile. Profiles can give configuration presets for different compute environments.

	  docker 	A generic configuration profile to be used with Docker, pulls software from Docker Hub: liuyangzzu/nanome:latest
	  singulairy	A generic configuration profile to be used with Singularity, pulls software from: docker://liuyangzzu/nanome:latest
	  conda		Please only use conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity. Create conda enviroment by 'conda env create -f environment.yml'
	  hpc		A generic configuration profile to be used on HPC cluster with SLURM job submission support.
	  google	A generic configuration profile to be used on Google Cloud platform with 'google-lifesciences' support.

	Contact to https://github.com/TheJacksonLaboratory/nanome/issues for bug report.
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

projectDir = workflow.projectDir
ch_utils = Channel.fromPath("${projectDir}/utils",  type: 'dir', followLinks: false)
ch_src   = Channel.fromPath("${projectDir}/src",  type: 'dir', followLinks: false)

// Reference genome, deepmod cluster settings
deepmod_tar_file = "${projectDir}/README.md"

referenceGenome = 'reference_genome/ref.fasta'
chromSizesFile = 'reference_genome/chrom.sizes'

// TODO: auto detect based on ref.fasta
if (params.dataType == 'human') {
	isDeepModCluster = params.useDeepModCluster
	if (isDeepModCluster) {
		deepmod_tar_file = params.deepmod_ctar
	}
} else if (params.dataType == 'ecoli') {
	isDeepModCluster = false
} else {
	println "Param dataType=${params.dataType} is not support"
	exit 1
}

// If need, preload C.tar.gz file in advance
deepmod_c_tar_ch = Channel.fromPath(deepmod_tar_file)


// if is true or 'true' (string), using '  '
chrSet = params.chrSet.toBoolean() ? '  ' : params.chrSet

log.info """\
NANOME - NF PIPELINE (v$workflow.manifest.version)
by Li Lab at The Jackson Laboratory
https://github.com/TheJacksonLaboratory/nanome
=================================
dsname			:${params.dsname}
input			:${params.input}
output			:${params.outputDir}
work			:${workflow.workDir}
dataType		:${params.dataType}
runBasecall		:${params.runBasecall}
runMethcall		:${params.runMethcall}
=================================
"""
.stripIndent()


workflow.onComplete {
	if (workflow.success && params.cleanCache) {
		def workDir = new File("${workflow.workDir}")
		println "rm -rf ${workflow.workDir}".execute().text
	}
}

// Channel for utils/ and src/ folders
ch_utils
	.into{	ch_utils1; ch_utils2; ch_utils3; ch_utils4;
			ch_utils5; ch_utils6; ch_utils7; ch_utils8;
			ch_utils9; ch_utils10}
ch_src
	.into{	ch_src1; ch_src2; ch_src3; ch_src4;
			ch_src5; ch_src_c1;ch_src_c2;ch_src_c3;
			ch_src_c4;ch_src_c5;ch_src_c6}

// Collect all folders of fast5 files, and send into Channels for pipelines
if (params.input.endsWith(".filelist.txt")) {
	// list of files in filelist.txt
	Channel.fromPath( params.input, checkIfExists: true )
		.splitCsv(header: false)
		.map { file(it[0]) }
		.set{ fast5_tar_ch }
} else if(params.input.endsWith("/*")) {
	// match all files in the folder, note: input must use '', prevent expand in advance
	// such as --input '/fastscratch/liuya/nanome/NA12878/NA12878_CHR22/input_chr22/*'
	Channel.fromPath(params.input, type: 'any').set{fast5_tar_ch}
} else {
	// For single file
	Channel.fromPath( params.input, checkIfExists: true ).set{fast5_tar_ch}
}

// Check all tools work well
process EnvCheck {
	tag 'EnvCheck'
	errorStrategy 'terminate'

	input:
	path reference_genome 	from 	Channel.fromPath(params.reference_genome, type: 'any')
	path("*") 				from 	ch_utils10
	path megalodonModelTar 	from 	Channel.fromPath(params.megalodon_model_tar)

	output:
	path "reference_genome" into reference_genome_ch
	path "megalodon_model"	optional true into megalodon_model_ch

	"""
	date; hostname; pwd
	echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"

	## Validate nanome container/environment is correct
	bash utils/validate_nanome_container.sh

	## Untar and prepare megalodon model
	if [[ ${params.runMegalodon} == "true" ]]; then
		tar -xzf ${megalodonModelTar}
		ls -lh megalodon_model
	fi

	## Get dir for reference_genome
	mkdir -p reference_genome
	find_dir="\$PWD/reference_genome"
	if [[ ${reference_genome} == *.tar.gz ]] ; then
		tar -xzf ${reference_genome} -C reference_genome
	elif [[ ${reference_genome} == *.tar ]] ; then
		tar -xf ${reference_genome} -C reference_genome
	else
		## for folder, use ln, note this is a symbolic link to a folder
		find_dir=\$( readlink -f ${reference_genome} )
	fi

	find \${find_dir} -name '*.fasta*' | \
		 parallel -j0 -v  'fn={/} ; ln -s -f  {}   reference_genome/\${fn/*.fasta/ref.fasta}'
	find \${find_dir} -name '*.sizes' | \
			parallel -j1 -v ln -s -f {} reference_genome/chrom.sizes

	ls -lh ${referenceGenome}
	ls -lh ${chromSizesFile}

	echo "### Check reference genome and chrSet"
	echo "referenceGenome=${referenceGenome}"
	echo "chromSizesFile=${chromSizesFile}"
	echo "chrSet=${chrSet}"
	echo "params.dataType=${params.dataType}"

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
	date; hostname; pwd
	echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"

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
	## find untarTempDir -name "*.fast5" -type f -exec mv {} ${fast5_tar.baseName}.untar/ \\;
	find untarTempDir -name "*.fast5" -type f | \
		parallel -j\$(( numProcessor ))  mv {}  ${fast5_tar.baseName}.untar/

	## Clean unused files
	rm -rf untarTempDir

	## Clean old analyses in input fast5 files
	if [[ "${params.cleanAnalyses}" == true ]] ; then
		echo "### Start cleaning old analysis"
		## python -c 'import h5py; print(h5py.version.info)'
		python utils/clean_old_basecall_in_fast5.py \
			-i ${fast5_tar.baseName}.untar --is-indir --verbose\
			--processor \$(( numProcessor * ${params.highProcTimes} ))
	fi

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

	input:
	path fast5_dir 				from 	untar_out_ch1
	each path("*") 				from 	ch_utils1

	output:
	path "${fast5_dir.baseName}.basecalled" into basecall_out_ch

	when:
	params.runBasecall

	"""
	date; hostname; pwd

	echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"
	if [[ "\${CUDA_VISIBLE_DEVICES:-}" == "" ]] ; then
		echo "Detect no GPU, using CPU commandType"
		commandType='cpu'
	else
		echo "Detect GPU, using GPU commandType"
		commandType='gpu'
	fi

	which guppy_basecaller
	guppy_basecaller -v
	mkdir -p ${fast5_dir.baseName}.basecalled

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path "${fast5_dir.baseName}.basecalled" \
			--config ${params.GUPPY_BASECALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out --compress_fastq\
			--verbose_logs
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path "${fast5_dir.baseName}.basecalled" \
			--config ${params.GUPPY_BASECALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out --compress_fastq\
			--verbose_logs \
			-x auto
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi

	## Combine fastq
	touch "${fast5_dir.baseName}.basecalled"/batch_basecall_combine_fq_${fast5_dir.baseName}.fq.gz

	## Below is compatable with both Guppy v4.2.2 (old) and newest directory structures
	find "${fast5_dir.baseName}.basecalled/" "${fast5_dir.baseName}.basecalled/pass/"\
	 	"${fast5_dir.baseName}.basecalled/fail/" -maxdepth 1 -name '*.fastq.gz' -type f\
	 	-print0 2>/dev/null | \
	 	while read -d \$'\0' file ; do
	 		cat \$file >> "${fast5_dir.baseName}.basecalled"/batch_basecall_combine_fq_${fast5_dir.baseName}.fq.gz
	 	done
	echo "### Combine fastq.gz DONE"

	## Remove fastq.gz
	find "${fast5_dir.baseName}.basecalled/"   "${fast5_dir.baseName}.basecalled/pass/"\
	 	"${fast5_dir.baseName}.basecalled/fail/" -maxdepth 1 -name '*.fastq.gz' -type f 2>/dev/null |\
	 	parallel -j\$(( numProcessor )) 'rm -f {}'

	## After basecall, rename and publish summary filenames, summary may also be used by resquiggle
	mv ${fast5_dir.baseName}.basecalled/sequencing_summary.txt \
		${fast5_dir.baseName}.basecalled/${fast5_dir.baseName}-sequencing_summary.txt

    # Clean
    if [[ ${params.cleanStep} == "true" ]]; then
    	rm -f ${fast5_dir.baseName}.basecalled.sam
    fi
	echo "### Basecalled by Guppy DONE"
	"""
}


// Duplicates basecall outputs
basecall_out_ch
	.into { qc_basecall_in_ch; resquiggle_in_ch; nanopolish_in_ch; deepmod_in_ch }


// Collect and output QC results for basecall, and report ONT coverage
process QCExport {
	tag "${params.dsname}"

	publishDir "${params.outputDir}/${params.dsname}-basecallings",
		mode: "copy", enabled: params.outputQC, overwrite: true

	input:
	path qc_basecall_list 		from 	qc_basecall_in_ch.collect()
	each path(reference_genome) from 	reference_genome_ch9

	output:
	path "${params.dsname}_basecall_report.html" into qc_out_ch
	path "${params.dsname}_QCReport" into qc_report_out_ch

	"""
	## Combine all sequencing summary files
	touch ${params.dsname}_combine_sequencing_summary.txt.gz
	firstFile=true
	find *.basecalled/ -name '*-sequencing_summary.txt' -type f -print0 |\
		while read -d \$'\0' file ; do
			if \$firstFile ; then
				awk 'NR>=1' \$file | \
					gzip -f >> ${params.dsname}_combine_sequencing_summary.txt.gz
				firstFile=false
			else
				awk 'NR>1' \$file | \
					gzip -f >> ${params.dsname}_combine_sequencing_summary.txt.gz
			fi
		done

	## Perform QC report by NanoComp
	NanoComp --summary ${params.dsname}_combine_sequencing_summary.txt.gz  \
		--names ${params.dsname} --outdir ${params.dsname}_QCReport -t \$(( numProcessor )) \
		--raw  -f pdf -p ${params.dsname}_

	## Combine all batch fq.gz
	> merge_all_fq.fq.gz
	cat *.basecalled/batch_basecall_combine_fq_*.fq.gz > merge_all_fq.fq.gz
	echo "### Fastq merge from all batches done!"

	## After basecall, we align results for ONT coverage analyses
	# align FASTQ files to reference genome, write sorted alignments to a BAM file
	minimap2 -t \$(( numProcessor * ${params.mediumProcTimes} )) -a  -x map-ont \
		${referenceGenome} \
		merge_all_fq.fq.gz | \
		samtools sort -@ \$(( numProcessor * ${params.mediumProcTimes} )) -T tmp -o \
			merge_all_bam.bam &&\
		samtools index -@ \$(( numProcessor * ${params.mediumProcTimes} ))  merge_all_bam.bam
    echo "### Samtools alignment done"

    ## calculates the sequence coverage at each position
    ## reporting genome coverage for all positions in BEDGRAPH format.
    bedtools genomecov -ibam merge_all_bam.bam -bg -strand + |
        awk '\$4 = \$4 FS "+"' |
        gzip -f > ${params.dsname}.coverage.positivestrand.bed.gz

    bedtools genomecov -ibam merge_all_bam.bam -bg -strand - |
        awk '\$4 = \$4 FS "-"' |
        gzip -f > ${params.dsname}.coverage.negativestrand.bed.gz

    cat ${params.dsname}.coverage.positivestrand.bed.gz > ${params.dsname}_ONT_coverage_combine.bed.gz
	cat ${params.dsname}.coverage.negativestrand.bed.gz >> ${params.dsname}_ONT_coverage_combine.bed.gz

	mv ${params.dsname}_ONT_coverage_combine.bed.gz ${params.dsname}_QCReport/
	mv ${params.dsname}_combine_sequencing_summary.txt.gz ${params.dsname}_QCReport/
	mv ${params.dsname}_QCReport/${params.dsname}_NanoComp-report.html ${params.dsname}_basecall_report.html

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f ${params.dsname}.coverage.positivestrand.bed.gz ${params.dsname}.coverage.negativestrand.bed.gz
		rm -f ${params.dsname}_merged.bam*
		rm -f merge_all_bam.bam*  merge_all_fq.fq.gz
	fi
    echo "### ONT coverage done!"
    echo "### QCReport all DONE"
	"""
}


// Resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
	tag "${basecallIndir.baseName}"

	input:
	path 	basecallIndir 			from resquiggle_in_ch
	each 	path(reference_genome) 	from reference_genome_ch3

	output:
	path "${basecallIndir.baseName}.resquiggle" into resquiggle_out_ch

	when:
	params.runMethcall && params.runResquiggle && !params.filterGPUTaskRuns

	"""
	### copy basecall workspace files, due to tombo resquiggle modify base folder
	rm -rf ${basecallIndir.baseName}.resquiggle
	mkdir -p ${basecallIndir.baseName}.resquiggle/workspace

	### original basecalled results will be parrallelly used by other processes
	cp -f ${basecallIndir}/batch_basecall_combine_fq_*.fq.gz  ${basecallIndir.baseName}.resquiggle/

	## cp -rf ${basecallIndir}/workspace  ${basecallIndir.baseName}.resquiggle/
	find ${basecallIndir}/workspace -name '*.fast5' -type f| \
		parallel -j\$(( numProcessor * ${params.mediumProcTimes} ))  \
		'cp {}   ${basecallIndir.baseName}.resquiggle/workspace/'
	echo "### Duplicate from basecall DONE"

	### Prerocessing, using combined fq.gz
	### ref: https://github.com/bioinfomaticsCSU/deepsignal#quick-start
	gunzip ${basecallIndir.baseName}.resquiggle/batch_basecall_combine_fq_*.fq.gz
	tombo preprocess annotate_raw_with_fastqs\
	 	--fast5-basedir ${basecallIndir.baseName}.resquiggle/workspace\
	 	--fastq-filenames ${basecallIndir.baseName}.resquiggle/batch_basecall_combine_fq_*.fq\
	 	--basecall-group ${params.BasecallGroupName}\
	 	--basecall-subgroup ${params.BasecallSubGroupName}\
	 	--overwrite --processes \$(( numProcessor * ${params.mediumProcTimes} ))  &>> Resquiggle.run.log
	echo "### tombo preprocess DONE"

	### Need to check Tombo resquiggle bugs, lots of users report long runtime and hang at nearly completion for large data
	### ref: https://github.com/nanoporetech/tombo/issues/139, https://github.com/nanoporetech/tombo/issues/111
	### ref: https://github.com/nanoporetech/tombo/issues/365, https://github.com/nanoporetech/tombo/issues/167
	### ref: https://nanoporetech.github.io/tombo/resquiggle.html?highlight=processes
	### Out of memory solution for large data: --tomboResquiggleOptions '--signal-length-range 0 500000  --sequence-length-range 0 50000'
	tombo resquiggle\
		--processes \$(( numProcessor * ${params.lowProcTimes} )) \
		--corrected-group ${params.ResquiggleCorrectedGroup} \
		--basecall-group ${params.BasecallGroupName} \
		--basecall-subgroup ${params.BasecallSubGroupName}\
		--ignore-read-locks ${params.tomboResquiggleOptions}\
		--overwrite \
		${basecallIndir.baseName}.resquiggle/workspace \
		${referenceGenome} &>> Resquiggle.run.log

	echo "### tombo resquiggle DONE"
	"""
}


// duplicate resquiggle results for DeepSignal and Tombo
resquiggle_out_ch.into { deepsignal_in_ch; tombo_in_ch }


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
	bamFileName="${params.dsname}.batch_${basecallDir.baseName}.sorted.bam"

	## Do alignment firstly, find the combined fastq file
	fastqFile=\$(find ${basecallDir}/ -name 'batch_basecall_combine_fq_*.fq.gz' -type f)

	# Index the raw read with fastq, we do not index in basecalled dir, in case of cache can be work
	ln -s \${fastqFile}  \${fastqFile##*/}

	## Index, ref: https://github.com/jts/nanopolish#data-preprocessing
	nanopolish index -d ${basecallDir}/workspace -s ${basecallDir}/${basecallDir.baseName}-sequencing_summary.txt  \${fastqFile##*/}

	## Aligning reads to the reference genome, ref: https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html#aligning-reads-to-the-reference-genome
	minimap2 -t \$(( numProcessor * ${params.mediumProcTimes} )) -a -x map-ont ${referenceGenome} \${fastqFile##*/} | \
		samtools sort -@ \$(( numProcessor * ${params.mediumProcTimes} )) -T tmp -o \${bamFileName} &&\
		samtools index -@ \$(( numProcessor * ${params.mediumProcTimes} ))  \${bamFileName}
	echo "### Alignment step: minimap2 and samtools DONE"

	## Calling methylation, ref: https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html#calling-methylation
	## there are segment fault issues, if set -t to a large number or use low memory,
	## ref: https://github.com/jts/nanopolish/issues/872
	## ref: https://github.com/jts/nanopolish/issues/683, https://github.com/jts/nanopolish/issues/580
	nanopolish call-methylation -t \$(( numProcessor/2 )) -r \${fastqFile##*/} \
		-b \${bamFileName} -g ${referenceGenome} | awk 'NR>1' | \
		gzip -f > batch_${basecallDir.baseName}.nanopolish.methylation_calls.tsv.gz

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f *.sorted.bam *.sorted.bam.bai
		rm -f *.fq.gz.index*
	fi
	echo "### Nanopolish methylation calling DONE"
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
	each path(megalodon_model_dir)	from megalodon_model_ch

	output:
	path "batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt.gz" into megalodon_out_ch

	when:
	params.runMethcall && params.runMegalodon

	"""
	date; hostname; pwd

	echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"
	if [[ "\${CUDA_VISIBLE_DEVICES:-}" == "" ]] ; then
		echo "Detect no GPU, using CPU commandType"
		commandType='cpu'
	else
		echo "Detect GPU, using GPU commandType"
		commandType='gpu'
	fi

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		## Ref: https://github.com/nanoporetech/megalodon
		## CPU issues: https://github.com/nanoporetech/megalodon/issues/172
		megalodon \
			${fast5_dir} \
			--overwrite \
			--outputs per_read_mods mods per_read_refs \
			--guppy-server-path guppy_basecall_server \
			--guppy-config ${params.MEGALODON_MODEL_FOR_GUPPY_CONFIG} \
			--guppy-params "-d ${megalodon_model_dir}/ --num_callers \$(( numProcessor )) --ipc_threads \$(( numProcessor * ${params.lowProcTimes} ))" \
			--guppy-timeout ${params.GUPPY_TIMEOUT} \
			--samtools-executable ${params.SAMTOOLS_PATH} \
			--sort-mappings \
			--mappings-format bam \
			--reference ${referenceGenome} \
			--mod-motif m CG 0 \
			--mod-output-formats bedmethyl wiggle \
			--write-mods-text \
			--write-mod-log-probs \
			--processes \$(( numProcessor ))  &>> Megalodon.run.log
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		## Ref: https://github.com/nanoporetech/megalodon
		megalodon \
			${fast5_dir} \
			--overwrite \
			--outputs per_read_mods mods per_read_refs \
			--guppy-server-path guppy_basecall_server \
			--guppy-config ${params.MEGALODON_MODEL_FOR_GUPPY_CONFIG} \
			--guppy-params "-d ${megalodon_model_dir}/ --num_callers \$(( numProcessor )) --ipc_threads \$(( numProcessor * ${params.highProcTimes} ))" \
			--guppy-timeout ${params.GUPPY_TIMEOUT} \
			--samtools-executable ${params.SAMTOOLS_PATH} \
			--sort-mappings \
			--mappings-format bam \
			--reference ${referenceGenome} \
			--mod-motif m CG 0 \
			--mod-output-formats bedmethyl wiggle \
			--write-mods-text \
			--write-mod-log-probs \
			--processes \$(( numProcessor * ${params.mediumProcTimes} )) \
			--devices 0  &>> Megalodon.run.log
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi

	awk 'NR>1' megalodon_results/per_read_modified_base_calls.txt | gzip -f > \
		batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt.gz

	### Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		### keep guppy server log, due to it may fail when remove that folder, rm -rf megalodon_results
		find megalodon_results/  -maxdepth 1 -type f |\
		 	parallel -j\$(( numProcessor )) 'rm {}'
	fi
	echo "### Megalodon DONE"
	"""
}


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${indir.baseName}"

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/deepsignal",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path indir 						from deepsignal_in_ch
	// each path(deepsignal_model_tar) from Channel.fromPath(params.deepsignal_model_tar)
	each path(reference_genome) 	from reference_genome_ch4

	output:
	path "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv.gz" into deepsignal_out_ch

	when:
	params.runMethcall && params.runDeepSignal && !params.filterGPUTaskRuns

	"""
	if ls /data/${params.DEEPSIGNAL_MODEL}* 1> /dev/null 2>&1; then
		DeepSignalModelBaseDir=/data
	else
		wget ${params.deepsignal_model_tar}  --no-verbose
		tar -xzf ${params.DEEPSIGNAL_MODEL_TAR_GZ} &&
			rm -f ${params.DEEPSIGNAL_MODEL_TAR_GZ}
		DeepSignalModelBaseDir="."
	fi

	commandType='gpu'
	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path "\${DeepSignalModelBaseDir}/${params.DEEPSIGNAL_MODEL}" \
			--result_file "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv" \
			--reference_path ${referenceGenome} \
			--corrected_group ${params.ResquiggleCorrectedGroup} \
			--nproc \$(( numProcessor * ${params.highProcTimes}  )) \
			--is_gpu no
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path "\${DeepSignalModelBaseDir}/${params.DEEPSIGNAL_MODEL}" \
			--result_file "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv" \
			--reference_path ${referenceGenome} \
			--corrected_group ${params.ResquiggleCorrectedGroup} \
			--nproc \$(( numProcessor * ${params.highProcTimes}  )) \
			--is_gpu yes
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi

	gzip -f batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv
	echo "### DeepSignal methylation DONE"
	"""
}


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
	path "outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz" optional true into guppy_methcall_gz_out_ch
	path "batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam*" into guppy_methcall_out_ch
	path "batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz" into guppy_gcf52ref_out_ch

	when:
	params.runMethcall && params.runGuppy

	"""
	date; hostname; pwd

	echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"
	if [[ "\${CUDA_VISIBLE_DEVICES:-}" == "" ]] ; then
		echo "Detect no GPU, using CPU commandType"
		commandType='cpu'
	else
		echo "Detect GPU, using GPU commandType"
		commandType='gpu'
	fi

	mkdir -p ${fast5_dir.baseName}.methcalled

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path ${fast5_dir.baseName}.methcalled \
			--config ${params.GUPPY_METHCALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out --compress_fastq\
			--verbose_logs
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path ${fast5_dir.baseName}.methcalled \
			--config ${params.GUPPY_METHCALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out --compress_fastq\
			--verbose_logs \
			--device auto
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi
	echo "### Guppy methylation calling DONE"

	## Extract guppy methylation-callings
	## Combine fastq
	touch batch_combine_fq.fq.gz

	find "${fast5_dir.baseName}.methcalled/" "${fast5_dir.baseName}.methcalled/pass/"\
	 	"${fast5_dir.baseName}.methcalled/fail/" -maxdepth 1 -name '*.fastq.gz' -type f\
	 	-print0 2>/dev/null | \
	 	while read -d \$'\0' file ; do
	 		cat \$file >> batch_combine_fq.fq.gz
	 	done

	## Remove fastq.gz
	find "${fast5_dir.baseName}.methcalled/"   "${fast5_dir.baseName}.methcalled/pass/"\
	 	"${fast5_dir.baseName}.methcalled/fail/" -maxdepth 1 -name '*.fastq.gz' \
	 	-type f 2>/dev/null |\
	 	parallel -j\$(( numProcessor )) 'rm -f {}'

	## gcf52ref ways
	minimap2 -t \$(( numProcessor * ${params.mediumProcTimes} )) -a -x map-ont ${referenceGenome} \
		batch_combine_fq.fq.gz | \
		samtools sort -@ \$(( numProcessor * ${params.mediumProcTimes} )) \
			-T tmp -o gcf52ref.batch.${fast5_dir.baseName}.bam &&\
		samtools index -@ \$(( numProcessor * ${params.mediumProcTimes} )) \
			gcf52ref.batch.${fast5_dir.baseName}.bam
	echo "### gcf52ref minimap2 alignment is done"

	## Modified version, support dir input, not all fast5 files (too long arguments)
	python utils/extract_methylation_fast5_support_dir.py \
		-p \$(( numProcessor * ${params.mediumProcTimes} )) ${fast5_dir.baseName}.methcalled/workspace
	echo "### gcf52ref extract to db done"

	## gcf52ref files preparation
	### git clone https://github.com/kpalin/gcf52ref.git
	tar -xzf utils/gcf52ref.tar.gz -C .
	patch gcf52ref/scripts/extract_methylation_from_rocks.py < utils/gcf52ref.patch

	python gcf52ref/scripts/extract_methylation_from_rocks.py \
		-d base_mods.rocksdb \
		-a gcf52ref.batch.${fast5_dir.baseName}.bam \
		-r ${referenceGenome} \
		-o tmp.batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv
	echo "### gcf52ref extract to tsv done"

	awk 'NR>1' tmp.batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv | gzip -f > \
		batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz
	echo "### gcf52ref extraction DONE"

	## fast5mod ways
	FAST5PATH=${fast5_dir.baseName}.methcalled/workspace
	OUTBAM=batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam

	fast5mod guppy2sam \${FAST5PATH} --reference ${referenceGenome} \
		--workers \$(( numProcessor * ${params.highProcTimes} )) --recursive --quiet \
		| samtools sort -@ \$(( numProcessor * ${params.mediumProcTimes} )) | \
		samtools view -b -@ \$(( numProcessor * ${params.mediumProcTimes} )) > \${OUTBAM} &&\
		samtools index -@ \$(( numProcessor * ${params.mediumProcTimes} ))  \${OUTBAM}

	if [[ "${params.outputIntermediate}" == true ]] ; then
		tar -czf outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz \
			batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam*
	fi
	echo "### fast5mod extraction DONE"

	## Clean
	## methcalled folder is no need, keep only gcf52ref's tsv and fast5mod's bam for combine step
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -rf ${fast5_dir.baseName}.methcalled
		rm -f gcf52ref.*.bam gcf52ref.*.bam.bai tmp*.tsv batch_combine_fq.fq.gz
		rm -rf gcf52ref/
		rm -rf base_mods.rocksdb/
		echo "### Clean DONE"
	fi
	echo "### Guppy methcall, extracted by fast5mod and gcf52ref DONE"
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

	when:
	params.runMethcall && params.runTombo && !params.filterGPUTaskRuns

	"""
	## Check if there is a BrokenPipeError: [Errno 32] Broken pipe
	## Ref: https://github.com/nanoporetech/tombo/issues/183
	## Note 1 is still fast for tombo
	tombo detect_modifications alternative_model \
		--fast5-basedirs ${resquiggleDir}/workspace \
		--dna\
		--statistics-file-basename batch_${resquiggleDir.baseName} \
		--per-read-statistics-basename batch_${resquiggleDir.baseName} \
		--alternate-bases CpG \
		--processes \$(( numProcessor )) \
		--corrected-group ${params.ResquiggleCorrectedGroup} \
		--multiprocess-region-size ${params.tomboMultiprocessRegionSize} &> \
		${resquiggleDir.baseName}.Tombo.run.log

	retry=1
	## while grep -q "BrokenPipeError:" ${resquiggleDir.baseName}.Tombo.run.log
	while ! tail -n 1 ${resquiggleDir.baseName}.Tombo.run.log |  grep -q "100%"
	do
		echo "### Found error in tombo detect_modifications, repeat tombo running again!!!"
		tombo detect_modifications alternative_model \
			--fast5-basedirs ${resquiggleDir}/workspace \
			--dna\
			--statistics-file-basename batch_${resquiggleDir.baseName} \
			--per-read-statistics-basename batch_${resquiggleDir.baseName} \
			--alternate-bases CpG \
			--processes \$(( numProcessor )) \
			--corrected-group ${params.ResquiggleCorrectedGroup} \
			--multiprocess-region-size ${params.tomboMultiprocessRegionSize} &> \
			${resquiggleDir.baseName}.Tombo.run.log
		retry=\$(( retry+1 ))
		if (( retry >= 5 )); then
			break
		fi
	done

	## if grep -q "BrokenPipeError: \\[Errno 32\\] Broken pipe" ${resquiggleDir.baseName}.Tombo.run.log; then
	if ! tail -n 1 ${resquiggleDir.baseName}.Tombo.run.log |  grep -q "100%" ; then
		## Grep the broken pipeline bug for Tombo
		echo "### Tombo seems not finish 100% after retry reached at \${retry} times, please check by yourself, it may be software or genome reference problem."
	else
		## Tombo was ok
		echo "### Tombo log passed, OK"
	fi

	## Tombo lib need h5py lower than 3.0
	## Error may occur with higher version of h5py: AttributeError: 'Dataset' object has no attribute 'value'
	## ref: https://github.com/nanoporetech/tombo/issues/325
	python utils/tombo_extract_per_read_stats.py \
		${chromSizesFile} \
		"batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats" \
		"batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats.bed"

	gzip -f batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats.bed

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f *.tombo.per_read_stats   *.tombo.stats
	fi
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
	path "batch_${basecallDir.baseName}_num.tar.gz" optional true  into deepmod_gz_out_ch

	when:
	params.runMethcall && params.runDeepMod && !params.filterGPUTaskRuns

	"""
	set +u
	source ${params.conda_base_dir}/etc/profile.d/conda.sh
	conda activate ${params.conda_name}
	set -u
	echo "### Env set ok"
	## Find the model dir
	DeepModTrainModelDir=\$(find \$CONDA_PREFIX -name 'train_deepmod' -type d)
	if [[ \${DeepModTrainModelDir:-} == "" ]]; then
		wget ${params.DeepModGithub} --no-verbose
		tar -xzf v0.1.3.tar.gz
		DeepModTrainModelDir="DeepMod-0.1.3/train_deepmod"
	fi

	if [[ "${params.dataType}" = "human" ]] ; then
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
			--Ref ${referenceGenome} \
			--Base C \
			--modfile \${DeepModTrainModelDir}/${params.DEEPMOD_RNN_MODEL} \
			--FileID batch_${basecallDir.baseName}_num \
			--threads \$(( numProcessor * ${params.mediumProcTimes} )) \
			${params.moveOption ? '--move' : ' '} &>> DeepMod.run.log

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


// prepare combining results
nanopolish_combine_in_ch = nanopolish_out_ch.toList()
deepmod_combine_in_ch=deepmod_out_ch.toList()
megalodon_combine_in_ch = megalodon_out_ch.toList()
tombo_combine_in_ch = tombo_out_ch.toList()
deepsignal_combine_in_ch = deepsignal_out_ch.toList()
guppy_combine_in_ch = guppy_methcall_out_ch.collect()


// Combine Nanopolish runs' all results together
process NplshComb {
	tag "${params.dsname}"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.*.combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x 		from nanopolish_combine_in_ch
	path ("*")	from ch_src_c1

	output:
	path "${params.dsname}.nanopolish.per_read.combine.tsv.gz" into nanopolish_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz" into read_unify_nanopolish
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz" into site_unify_nanopolish

	when:
	x.size() >= 1 && params.runCombine

	"""
	> ${params.dsname}.nanopolish.per_read.combine.tsv.gz
	cat ${x} > ${params.dsname}.nanopolish.per_read.combine.tsv.gz

	## Unify format output
	bash src/unify_format_for_calls.sh \
		${params.dsname}  Nanopolish ${params.dsname}.nanopolish.per_read.combine.tsv.gz \
		.  \$((numProcessor))  12  ${chrSet}

	echo "### Nanopolish combine DONE"
	"""
}


// Combine Megalodon runs' all results together
process MgldnComb {
	tag "${params.dsname}"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.*.combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x from megalodon_combine_in_ch
	path ("*")	from ch_src_c2

	output:
	path "${params.dsname}.megalodon.per_read.combine.bed.gz" into megalodon_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz" into read_unify_megalodon
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz" into site_unify_megalodon

	when:
	x.size() >= 1  && params.runCombine

	"""
	> ${params.dsname}.megalodon.per_read.combine.bed.gz
	cat ${x} > ${params.dsname}.megalodon.per_read.combine.bed.gz

	## Unify format output
	bash src/unify_format_for_calls.sh \
		${params.dsname}  Megalodon ${params.dsname}.megalodon.per_read.combine.bed.gz \
		.  \$((numProcessor))  12  ${chrSet}

	echo "### Megalodon combine DONE"
	"""
}


// Combine DeepSignal runs' all results together
process DpSigComb {
	tag "${params.dsname}"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.*.combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x from deepsignal_combine_in_ch
	path ("*")	from ch_src_c3

	output:
	path "${params.dsname}.deepsignal.per_read.combine.tsv.gz" into deepsignal_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz" into read_unify_deepsignal
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz" into site_unify_deepsignal

	when:
	x.size() >= 1 && params.runCombine

	"""
	touch ${params.dsname}.deepsignal.per_read.combine.tsv.gz
	cat ${x} > ${params.dsname}.deepsignal.per_read.combine.tsv.gz

	## Unify format output
	bash src/unify_format_for_calls.sh \
		${params.dsname}  DeepSignal ${params.dsname}.deepsignal.per_read.combine.tsv.gz \
		.  \$((numProcessor))  12  ${chrSet}
	echo "### DeepSignal combine DONE"
	"""
}


// Combine Guppy runs' all results together
process GuppyComb {
	tag "${params.dsname}"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}.guppy.*.combine.tsv.gz",
		enabled: params.outputRaw

	publishDir "${params.outputDir}/${params.dsname}_raw_outputs/guppy",
		mode: "copy", pattern: "${params.dsname}.guppy_fast5mod.combined.bam.tar.gz",
		enabled: params.outputIntermediate

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x 					from guppy_combine_in_ch
	path y 					from guppy_gcf52ref_out_ch.collect()
	path reference_genome 	from reference_genome_ch8
	path ("*")				from ch_src_c4

	output:
	path "${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz" into guppy_fast5mod_combine_out_ch
	path "${params.dsname}.guppy.gcf52ref_per_read.combine.tsv.gz" into guppy_gcf52ref_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz" into read_unify_guppy
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz" into site_unify_guppy
	path "${params.dsname}.guppy_fast5mod.combined.bam.tar.gz" optional true into guppy_combine_raw_out_ch

	when:
	x.size() >= 1 && params.runCombine

	"""
	## gcf52ref ways
	cat batch_*.guppy.gcf52ref_per_read.tsv.gz > ${params.dsname}.guppy.gcf52ref_per_read.combine.tsv.gz
	echo "### gcf52ref combine DONE"

	## fast5mod ways combine
	## find name like batch_*.guppy.fast5mod_guppy2sam.bam*
	find . -maxdepth 1  -name 'batch_*.guppy.fast5mod_guppy2sam.bam' |
		parallel -j\$(( numProcessor )) --xargs -v \
		samtools merge -@\$(( numProcessor * ${params.mediumProcTimes} )) total.meth.bam {}

	### sort is not needed due to merge the sorted bam, ref: http://www.htslib.org/doc/samtools-merge.html
	### samtools sort -@ \$(( numProcessor * ${params.mediumProcTimes} )) total.meth.bam
	samtools index -@ \$(( numProcessor * ${params.mediumProcTimes} )) total.meth.bam
	echo "Samtool merge and index for fast5mod DONE"

	if [[ ${params.outputIntermediate} == true ]] ; then
		tar -czf ${params.dsname}.guppy_fast5mod.combined.bam.tar.gz total.meth.bam*
	fi

	awk '/^>/' ${referenceGenome} | awk '{print \$1}' \
		> rf_chr_all_list.txt
	if [[ "${params.dataType}" == "human" ]] ; then
		echo "### For human, extract chr1-22, X and Y"
		> chr_all_list.txt
		for i in {1..22} X Y
		do
			if cat rf_chr_all_list.txt | grep -w ">chr\${i}" -q ; then
				echo chr\${i} >> chr_all_list.txt
			fi
		done
		rm  -f rf_chr_all_list.txt
		echo "### Chomosome list"
		cat chr_all_list.txt

		## Ref: https://github.com/nanoporetech/medaka/issues/177
		parallel -j\$(( numProcessor * ${params.highProcTimes} )) -v \
			"fast5mod call total.meth.bam ${referenceGenome} \
        		meth.chr_{}.tsv  --meth cpg --quiet --regions {} ; \
        		gzip -f meth.chr_{}.tsv" :::: chr_all_list.txt

		touch ${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz
        for i in {1..22} X Y
		do
			if [ -f "meth.chr_chr\$i.tsv.gz" ]; then
				cat  meth.chr_chr\$i.tsv.gz >> \
					${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz
			fi
		done
	elif [[ "${params.dataType}" == "ecoli" ]] ; then
		echo "### For ecoli, chr=${chrSet}"
		fast5mod call total.meth.bam ${referenceGenome} \
			meth.chr_${chrSet}.tsv \
			--meth cpg --quiet \
			--regions ${chrSet}
		gzip -f  meth.chr_${chrSet}.tsv && \
			mv meth.chr_${chrSet}.tsv.gz ${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz
	fi

	## Unify format output for read level
	bash src/unify_format_for_calls.sh \
		${params.dsname}  Guppy.gcf52ref ${params.dsname}.guppy.gcf52ref_per_read.combine.tsv.gz \
		.  \$((numProcessor))  1  ${chrSet}

	## Unify format output for site level
	bash src/unify_format_for_calls.sh \
		${params.dsname}  Guppy ${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz \
		.  \$((numProcessor))  2  ${chrSet}

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f meth.chr*.tsv.gz
		rm -f total.meth.bam*
	fi
	echo "### Guppy combine DONE"
	"""
}


// Combine Tombo runs' all results together
process TomboComb {
	tag "${params.dsname}"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.*.combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x from tombo_combine_in_ch // list of tombo bed files
	path ("*")	from ch_src_c5

	output:
	path "${params.dsname}.tombo.per_read.combine.bed.gz" into tombo_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz" into read_unify_tombo
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz" into site_unify_tombo

	when:
	x.size() >= 1 && params.runCombine

	"""
	touch ${params.dsname}.tombo.per_read.combine.bed.gz
	cat ${x} > ${params.dsname}.tombo.per_read.combine.bed.gz

	## Unify format output
	bash src/unify_format_for_calls.sh \
		${params.dsname}  Tombo ${params.dsname}.tombo.per_read.combine.bed.gz \
		.  \$((numProcessor))  12 ${chrSet}
	echo "### Tombo combine DONE"
	"""
}


// Combine DeepMod runs' all results together
process DpmodComb {
	tag "${params.dsname}"

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

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x 					from deepmod_combine_in_ch
	path deepmod_c_tar_file from deepmod_c_tar_ch
	each path("*") 			from ch_utils6
	path ("*")				from ch_src_c6

	output:
	path "${params.dsname}.deepmod.*.combine.bed.gz" into deepmod_combine_out_ch
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz" into site_unify_deepmod
	path "${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz" optional true into deepmod_combine_sum_chrs_mod_ch
	path "${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz" optional true into deepmod_combine_c_cluster_all_chrs_ch
	path "${params.dsname}.deepmod.all_batch.C.bed.tar.gz" optional true  into deepmod_combine_all_batch_c_ch

	when:
	x.size() >= 1 && params.runCombine

	"""
	set +u
	source ${params.conda_base_dir}/etc/profile.d/conda.sh
	conda activate ${params.conda_name}
	set -u
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
		indir/ C ${params.dsname}.deepmod ${chrSet}  &>> DpmodComb.run.log

	> ${params.dsname}.deepmod.C_per_site.combine.bed

	## Note: for ecoli data, no pattern for chr*, but N*
	for f in \$(ls -1 indir/${params.dsname}.deepmod.*.C.bed)
	do
	  cat \$f >> ${params.dsname}.deepmod.C_per_site.combine.bed
	done
	gzip -f ${params.dsname}.deepmod.C_per_site.combine.bed

	if [[ "${params.dataType}" == "human" && "${isDeepModCluster}" == "true" ]] ; then
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

		## consider modification cluster effect.
		## ref: https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md#3-how-to-consider-modification-cluster-effect
		python utils/hm_cluster_predict.py \
			indir/${params.dsname}.deepmod \
			./C \
			\${DeepModTrainModelDir}/${params.DEEPMOD_CLUSTER_MODEL} &>> DpmodComb.run.log

		> ${params.dsname}.deepmod.C_clusterCpG_per_site.combine.bed
		for f in \$(ls -1 indir/${params.dsname}.deepmod_clusterCpG.*.C.bed)
		do
		  cat \$f >> ${params.dsname}.deepmod.C_clusterCpG_per_site.combine.bed
		done

		gzip -f ${params.dsname}.deepmod.C_clusterCpG_per_site.combine.bed

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

	if [[ "${isDeepModCluster}" == "true" ]] ; then
		callfn=${params.dsname}.deepmod.C_clusterCpG_per_site.combine.bed.gz
	else
		callfn=${params.dsname}.deepmod.C_per_site.combine.bed.gz
	fi

	## Unify format output
	bash src/unify_format_for_calls.sh \
		${params.dsname}  DeepMod \${callfn} \
		.  \$((numProcessor))  2  ${chrSet}  &>> DpmodComb.run.log

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -rf indir/
		echo "### Clean DONE"
	fi
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
	tag "${params.dsname}"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz",
		enabled: params.outputRaw

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path naonopolish	from 	read_unify_nanopolish
	path megalodon 		from 	read_unify_megalodon
	path deepsignal 	from 	read_unify_deepsignal
	path("*") 		from 	ch_src3
	path("*") 		from 	ch_utils8

	output:
	path "${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz" optional true into meteore_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz" optional true into read_unify_meteore
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz" optional true into site_unify_meteore

	when:
	params.runMETEORE && ! params.filterGPUTaskRuns

	"""
	## METEORE outputs by combining other tools
	outFileName=${params.dsname}_Megalodon_DeepSignal_combine.model_content.tsv
	> \$outFileName
	printf '%s\t%s\n' deepsignal ${deepsignal} >> \$outFileName
	printf '%s\t%s\n' megalodon ${megalodon} >> \$outFileName

	METEOREDIR=\$(find /data -maxdepth 1 -name "${params.METEORE_Dir}" -type d)
	if [[ \$METEOREDIR == "" ]]; then
		wget ${params.METEOREGithub}  --no-verbose
		tar -xzf v1.0.0.tar.gz
		METEOREDIR=\$(find . -name "${params.METEORE_Dir}" -type d)
	fi

	## Degrade sk-learn for METEORE program if needed, it's model load need lower version
	## pip install -U scikit-learn==0.21.3
	combineScript=utils/combination_model_prediction.py

	modelContentFileName=\$(find . -name "${params.dsname}_Megalodon_DeepSignal_combine.model_content.tsv" -type f)
	# Use the optimized model (n_estimator = 3 and max_dep = 10)
	# Please note this optimized model is reported in METEORE paper, ref: https://github.com/comprna/METEORE#command
	# paper: https://doi.org/10.1101/2020.10.14.340315, text:  random forest (RF) (parameters: max_depth=3 and n_estimator=10)
	# therefore, we output optimized model results
	python \${combineScript} \
		-i \${modelContentFileName} \
		-m optimized -b \${METEOREDIR} \
		-o ${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz

	# Use default parameters from the sklearn library (n_estimator = 100 and max_dep = None)
	python \${combineScript} \
		-i \${modelContentFileName} -m default -b \${METEOREDIR} \
		-o ${params.dsname}.meteore.megalodon_deepsignal_default_rf_model_per_read.combine.tsv.gz

	# Read level and site level output
	if [ -f ${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz ] ; then
		mkdir -p Read_Level-${params.dsname}
		zcat ${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz | \
			awk -F '\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$6,\$5}' |
			gzip -f > Read_Level-${params.dsname}/${params.dsname}_METEORE-perRead-score.tsv.gz

		## Unify format output for site level
		bash src/unify_format_for_calls.sh \
			${params.dsname}  METEORE\
			${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz \
			.  \$((numProcessor))  2  ${chrSet}
	fi
	echo "### METEORE consensus DONE"
	"""
}


meteore_combine_out_ch
	.concat(raw_results_five_tools2.flatten(),
			deepmod_combine_out_ch.flatten(),
			guppy_fast5mod_combine_out_ch)
	.toList()
	.into { all_raw_results1; all_raw_results2; all_raw_results3; all_raw_results4}


qc_report_out_ch
	.concat(
		site_unify_nanopolish, site_unify_megalodon, site_unify_deepsignal,
		site_unify_guppy, site_unify_tombo, site_unify_deepmod, site_unify_meteore
	)
	.toList()
	.set {report_in_ch}


// Not cache due to the script contains run information, each time of resume run will need updated
process Report {
	tag "${params.dsname}"

	publishDir "${params.outputDir}", mode: "copy"

	input:
	path fileList 	from 	report_in_ch
	path("*") 		from 	ch_src5
	path("*") 		from 	ch_utils9

	output:
	path "${params.dsname}_nanome_report.html" 	into	report_out_ch
	path "README_${params.dsname}.txt" 			into 	readme_out_ch

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
	mkdir -p ${params.dsname}_NANOME_report
	cp src/nanocompare/report/style.css ${params.dsname}_NANOME_report/
	cp -rf src/nanocompare/report/icons ${params.dsname}_NANOME_report/
	cp -rf src/nanocompare/report/js ${params.dsname}_NANOME_report/

	PYTHONPATH=src  python src/nanocompare/report/gen_html_report.py \
		${params.dsname} \
		running_information.tsv \
		\${basecallOutputFile} \
		. \
		${params.dsname}_NANOME_report

	## No tty usage, ref: https://github.com/remy/inliner/issues/151
	script -qec "inliner ${params.dsname}_NANOME_report/nanome_report.html" /dev/null  > ${params.dsname}_nanome_report.html

	PYTHONIOENCODING=UTF-8 python utils/gen_readme.py\
		utils/readme.txt.template ${params.dsname} ${params.outputDir}\
		${workflow.projectDir} ${workflow.workDir} "${workflow.commandLine}"\
		${workflow.runName} "${workflow.start}"\
		> README_${params.dsname}.txt

	echo "### report html DONE"
	"""
}
