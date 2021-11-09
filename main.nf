#!/usr/bin/env nextflow
/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/TheJacksonLaboratory/nanome
 @Author   : Yang Liu
 @FileName : main.nf
 @Software : NANOME project
 @Organization : JAX Li Lab
----------------------------------------------------------------------------------------
*/
// We now support both latest and lower versions, due to Lifebit CloudOS is only support 21.04
// Note: NXF_VER=20.04.1 nextflow run main.nf -profile test,singularity
if( ! nextflow.version.matches(">= 20.07.1") ){
	nextflow.preview.dsl=2
} else {
	// Support lower version of nextflow
	nextflow.enable.dsl=2
}

def helpMessage() {
	log.info"""
	NANOME - Nextflow PIPELINE (v$workflow.manifest.version)
	by Li Lab at The Jackson Laboratory
	https://github.com/TheJacksonLaboratory/nanome
	=================================
	Usage:
	The typical command for running the pipeline is as follows:

	nextflow run TheJacksonLaboratory/nanome -profile test,docker
	nextflow run TheJacksonLaboratory/nanome -profile test,singularity

	Mandatory arguments:
	  --dsname		Dataset name
	  --input		Input path for raw fast5 folders/tar/tar.gz files
	  --genome		Genome reference name ('hg38', 'hg38_chr22', or 'ecoli'), or a directory, or a tar.gz file

	General options:
	  --processors		Processors used for each task
	  --outdir		Output dir, default is 'outputs'
	  --type		Data type, default is 'human', can be also 'ecoli'
	  --chrSet		Chromosomes used in analysis, default is true, means chr1-22, X and Y, seperated by comma. For E. coli data, it needs be set to 'NC_000913.3'

	  --cleanCache		If clean work dir after complete, default is true

	Running environment options:
	  --docker_name		Docker name used for pipeline, default is 'liuyangzzu/nanome:latest'
	  --singularity_name	Singularity name used for pipeline, default is 'docker://liuyangzzu/nanome:latest'
	  --singularity_cache	Singularity cache dir, default is 'local_singularity_cache'
	  --conda_name		Conda name used for pipeline, default is 'nanome'
	  --conda_base_dir	Conda base directory, default is '/opt/conda'

	Platform specific options:
	  --queue		SLURM job submission queue name for cluster running, default is 'gpu'
	  --qos			SLURM job submission qos name for cluster running, default is 'inference'
	  --gresOptions		SLURM job submission GPU allocation options for cluster running, default is 'gpu:v100:1'
	  --time		SLURM job submission time allocation options for cluster running, default is '2h'
	  --memory		SLURM job submission memory allocation options for cluster running, default is '32GB'

	  --googleProjectName	Google Cloud project name for google-lifesciences task running

	Tools's specific configurations:
	  --run[Tool-name]	Default we run top four performers in nanome paper, specify '--run[Tool-name]' can include other tool, supported tools: Megalodon, Nanopolish, DeepSignal, Guppy, Tombo, METEORE, and DeepMod

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
if (params.help){
    helpMessage()
    exit 0
}

// Check mandatory params
if (params.dsname == false) { exit 1, "Missing --dsname option for dataset name, check command help use --help" }
if (params.input == false) { exit 1, "Missing --input option for input data, check command help use --help" }

// Parse genome params
genome_map = params.genome_map
// online input, or google storage input
megalodon_model_tar = params.megalodon_model_tar

if (genome_map[params.genome] != null) { genome_path = genome_map[params.genome] } else { 	genome_path = params.genome }

// Get src and utils dir
projectDir = workflow.projectDir
ch_utils = Channel.fromPath("${projectDir}/utils",  type: 'dir', followLinks: false)
ch_src   = Channel.fromPath("${projectDir}/src",  type: 'dir', followLinks: false)

// Reference genome, deepmod cluster settings
def deepmod_tar_file = "${projectDir}/README.md"
def referenceGenome = 'reference_genome/ref.fasta'
def chromSizesFile = 'reference_genome/chrom.sizes'

if (params.type == 'human') {
	isDeepModCluster = params.useDeepModCluster
	if (isDeepModCluster && params.runDeepMod) {
		deepmod_tar_file = params.deepmod_ctar
	}
} else if (params.type == 'ecoli') { isDeepModCluster = false }
else { 	exit 1, "Param type=${params.type} is not support" }

// if is true or 'true' (string), using '  '
chrSet = params.chrSet.toBoolean() ? '  ' : params.chrSet

workflow.onComplete {
	if (workflow.success && params.cleanCache) {
		def workDir = new File("${workflow.workDir}")
		println "rm -rf ${workflow.workDir}".execute().text
	}
}


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
	Channel.fromPath(params.input, type: 'any').set{ fast5_tar_ch }
} else {
	// For single file/wildcard matched files
	Channel.fromPath( params.input, checkIfExists: true ).set{ fast5_tar_ch }
}

// Header log info
def summary = [:]
summary['dsname'] 			= params.dsname
summary['input'] 			= params.input
summary['genome'] 			= params.genome

summary['\nRunning settings']         = "--------"
summary['processors'] 		= params.processors
if (params.runNanopolish) summary['runNanopolish'] = 'Yes'
if (params.runMegalodon) summary['runMegalodon'] = 'Yes'
if (params.runDeepSignal) summary['runDeepSignal'] = 'Yes'
if (params.runGuppy) summary['runGuppy'] = 'Yes'
if (params.runTombo) summary['runTombo'] = 'Yes'
if (params.runMETEORE) summary['runMETEORE'] = 'Yes'
if (params.runDeepMod) summary['runDeepMod'] = 'Yes'

// summary['megalodon_model_tar'] = megalodon_model_tar

summary['\nPipeline settings']         = "--------"
summary['Working dir'] 		= workflow.workDir
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Profile']          = workflow.profile
summary['Config Files'] 	= workflow.configFiles.join(',')
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
if (workflow.profile.contains('hpc') || workflow.profile.contains('winter') || workflow.profile.contains('sumner') ) {
	summary['\nHPC settings']         = "--------"
    summary['queue']         = params.queue
    summary['qos']          = params.qos
    summary['memory']            = params.memory
    summary['time']            = params.time
    if (params.gresOptions != false) {summary['gresOptions'] = params.gresOptions }
}
if (workflow.profile.contains('google') || params.config.contains('lifebit')) {
	summary['\nGCP settings']         = "--------"
	if (!params.config.contains('lifebit')) {
		summary['googleProjectName']    = params.googleProjectName
	} else { // lifebit specific settings
		summary['config']       		= params.config
		summary['zoneCloud']       		= params.zoneCloud
		summary['networkLifebit']       = params.networkLifebit
		summary['subnetworkLifebit']	= params.subnetworkLifebit
	}
    summary['googleLocation']          = params.googleLocation
    summary['googleRegion']            = params.googleRegion
    summary['bootDiskSizeCloud']       = params.bootDiskSizeCloud

	summary['machineType']         	= params.machineType
	summary['highmemMachineType']  	= params.highmemMachineType
	summary['gpuType']         	= params.gpuType
	summary['gpuNumber']        = params.gpuNumber

	summary['lowDiskSize']      = params.lowDiskSize
	summary['midDiskSize']      = params.midDiskSize
	summary['highDiskSize']     = params.highDiskSize

	summary['errorStrategy']    = params.errorStrategy
	summary['maxRetries']       = params.maxRetries
	summary['echo']         	= params.echo
}

log.info """\
NANOME - NF PIPELINE (v$workflow.manifest.version)
by Li Lab at The Jackson Laboratory
https://github.com/TheJacksonLaboratory/nanome
================================="""
.stripIndent()

log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "================================="

// Check all tools work well
process EnvCheck {
	tag 'EnvCheck'
	errorStrategy 'terminate'

	input:
	path reference_genome
	path megalodonModelTar

	output:
	path "reference_genome",	emit: reference_genome
	path "megalodon_model",		emit: megalodon_model, optional: true

	script:
	"""
	date; hostname; pwd
	echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"

	## Validate nanome container/environment is correct
	validate_nanome_container.sh

	## Untar and prepare megalodon model
	if [[ ${params.runMegalodon} == "true" ]]; then
		ls -lh ${megalodonModelTar}

		tar -xzf ${megalodonModelTar}
		## Check Megalodon model
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
	echo "params.type=${params.type}"

	echo "### Check env DONE"
	"""
}


// Untar of subfolders named 'M1', ..., 'M10', etc.
process Untar {
	tag "${fast5_tar.baseName}"

	input:
	path fast5_tar

	output:
	path "${fast5_tar.baseName}.untar", emit:untar,  optional: true

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
		## cp -rf ${fast5_tar}/* untarTempDir/  || true # failed means nothing in this folder
		find ${fast5_tar}/ -name '*.fast5' | \
			parallel -j\$(( numProcessor ))  cp {} untarTempDir/
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
		clean_old_basecall_in_fast5.py \
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


// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
	tag "${fast5_dir.baseName}"

	input:
	path fast5_dir

	output:
	path "${fast5_dir.baseName}.basecalled", 	emit: basecall

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
			--verbose_logs  &>> Basecall.run.log
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path "${fast5_dir.baseName}.basecalled" \
			--config ${params.GUPPY_BASECALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out --compress_fastq\
			--verbose_logs \
			-x auto  &>> Basecall.run.log
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


// Collect and output QC results for basecall, and report ONT coverage
process QCExport {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-basecallings",
		mode: "copy", enabled: params.outputQC, overwrite: true

	input:
	path qc_basecall_list
	path reference_genome

	output:
	path "${params.dsname}_basecall_report.html",	emit: qc_html
	path "${params.dsname}_QCReport",				emit: qc_report
	path "${params.dsname}_merge_all_bam.bam*",		optional: true,	 emit: merge_all_bam

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
		--raw  -f pdf -p ${params.dsname}_   &>> QCReport.run.log

	if [[ ${params.outputBam} == true  || ${params.outputONTCoverage} == true ]]; then
		## Combine all batch fq.gz
		> merge_all_fq.fq.gz
		cat *.basecalled/batch_basecall_combine_fq_*.fq.gz > merge_all_fq.fq.gz
		echo "### Fastq merge from all batches done!"

		## After basecall, we align results to merged, sorted bam, can be for ONT coverage analyses/output bam
		# align FASTQ files to reference genome, write sorted alignments to a BAM file
		minimap2 -t \$(( numProcessor * ${params.mediumProcTimes} )) -a  -x map-ont \
			${referenceGenome} \
			merge_all_fq.fq.gz | \
			samtools sort -@ \$(( numProcessor * ${params.mediumProcTimes} )) -T tmp -o \
				${params.dsname}_merge_all_bam.bam &&\
			samtools index -@ \$(( numProcessor * ${params.mediumProcTimes} ))  ${params.dsname}_merge_all_bam.bam
		echo "### Samtools alignment done"
	fi

	if [[ ${params.outputONTCoverage} == true ]]; then
		## calculates the sequence coverage at each position
		## reporting genome coverage for all positions in BEDGRAPH format.
		bedtools genomecov -ibam ${params.dsname}_merge_all_bam.bam -bg -strand + |
			awk '\$4 = \$4 FS "+"' |
			gzip -f > ${params.dsname}.coverage.positivestrand.bed.gz

		bedtools genomecov -ibam ${params.dsname}_merge_all_bam.bam -bg -strand - |
			awk '\$4 = \$4 FS "-"' |
			gzip -f > ${params.dsname}.coverage.negativestrand.bed.gz

		cat ${params.dsname}.coverage.positivestrand.bed.gz > ${params.dsname}_ONT_coverage_combine.bed.gz
		cat ${params.dsname}.coverage.negativestrand.bed.gz >> ${params.dsname}_ONT_coverage_combine.bed.gz

		mv ${params.dsname}_ONT_coverage_combine.bed.gz ${params.dsname}_QCReport/
	fi

	mv ${params.dsname}_combine_sequencing_summary.txt.gz ${params.dsname}_QCReport/
	mv ${params.dsname}_QCReport/${params.dsname}_NanoComp-report.html ${params.dsname}_basecall_report.html

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f ${params.dsname}.coverage.positivestrand.bed.gz ${params.dsname}.coverage.negativestrand.bed.gz
		rm -f merge_all_fq.fq.gz
		if [[ ${params.outputBam} == false ]]; then
			rm -f ${params.dsname}_merge_all_bam.bam*
		fi
	fi
    echo "### ONT coverage done!"
    echo "### QCReport all DONE"
	"""
}


// Resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
	tag "${basecallIndir.baseName}"

	input:
	path 	basecallIndir
	each 	path(reference_genome)

	output:
	path "${basecallIndir.baseName}.resquiggle", 	emit: resquiggle

	when:
	params.runMethcall && (params.runDeepSignal || params.runTombo)

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
		--processes \$( echo "print(int( \$numProcessor * ${params.reduceProcTimes} ))" | python3 ) \
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


// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Nanopolish {
	tag "${basecallDir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_raw_outputs/nanopolish",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path basecallDir
	each path(reference_genome)

	output:
	path "batch_${basecallDir.baseName}.nanopolish.methylation_calls.tsv.gz", 	emit: nanopolish_out

	when:
	params.runMethcall && params.runNanopolish

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
	nanopolish call-methylation -t \$( echo "print(int( \$numProcessor * ${params.reduceProcTimes} ))" | python3 ) -r \${fastqFile##*/} \
		-b \${bamFileName} -g ${referenceGenome} -q cpg | awk 'NR>1' | \
		gzip -f > batch_${basecallDir.baseName}.nanopolish.methylation_calls.tsv.gz

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f *.sorted.bam *.sorted.bam.bai
		rm -f *.fq.gz.index*
	fi
	echo "### Nanopolish methylation calling DONE"
	"""
}


// Combine Nanopolish runs' all results together
process NplshComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.*.combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x

	output:
	path "${params.dsname}.nanopolish.per_read.combine.tsv.gz",	emit: nanopolish_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz",	emit: site_unify

	when:
	x.size() >= 1 && params.runCombine

	"""
	> ${params.dsname}.nanopolish.per_read.combine.tsv.gz
	cat ${x} > ${params.dsname}.nanopolish.per_read.combine.tsv.gz

	## Unify format output
	unify_format_for_calls.sh \
		${params.dsname}  Nanopolish ${params.dsname}.nanopolish.per_read.combine.tsv.gz \
		.  \$((numProcessor))  12  ${chrSet}

	echo "### Nanopolish combine DONE"
	"""
}


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${fast5_dir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_raw_outputs/megalodon",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path fast5_dir
	each path(reference_genome)
	each path(megalodon_model_dir)

	output:
	path "batch_${fast5_dir.baseName}.megalodon.per_read_modified_base_calls.txt.gz", emit: megalodon_out

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


// Combine Megalodon runs' all results together
process MgldnComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.*.combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x

	output:
	path "${params.dsname}.megalodon.per_read.combine.bed.gz",	emit: megalodon_combine
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz",	emit: site_unify

	when:
	x.size() >= 1  && params.runCombine

	"""
	> ${params.dsname}.megalodon.per_read.combine.bed.gz
	cat ${x} > ${params.dsname}.megalodon.per_read.combine.bed.gz

	## Unify format output
	unify_format_for_calls.sh \
		${params.dsname}  Megalodon ${params.dsname}.megalodon.per_read.combine.bed.gz \
		.  \$((numProcessor))  12  ${chrSet}

	echo "### Megalodon combine DONE"
	"""
}


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${indir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_raw_outputs/deepsignal",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path indir
	each path(reference_genome)

	output:
	path "batch_${indir.baseName}.CpG.deepsignal.call_mods.tsv.gz",	emit: deepsignal_out

	when:
	params.runMethcall && params.runDeepSignal

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


// Combine DeepSignal runs' all results together
process DpSigComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.*.combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x

	output:
	path "${params.dsname}.deepsignal.per_read.combine.tsv.gz",	emit: deepsignal_combine_out
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz",	emit: site_unify

	when:
	x.size() >= 1 && params.runCombine

	"""
	touch ${params.dsname}.deepsignal.per_read.combine.tsv.gz
	cat ${x} > ${params.dsname}.deepsignal.per_read.combine.tsv.gz

	## Unify format output
	unify_format_for_calls.sh \
		${params.dsname}  DeepSignal ${params.dsname}.deepsignal.per_read.combine.tsv.gz \
		.  \$((numProcessor))  12  ${chrSet}
	echo "### DeepSignal combine DONE"
	"""
}


// methylation calling for Guppy
process Guppy {
	tag "${fast5_dir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_raw_outputs/guppy",
		mode: "copy",
		pattern: "outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz",
		enabled: params.outputIntermediate
	publishDir "${params.outdir}/${params.dsname}_raw_outputs/guppy",
		mode: "copy",
		pattern: "batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz",
		enabled: params.outputIntermediate

	input:
	path fast5_dir
	each path(reference_genome)
	each path(utils)

	output:
	path "outbatch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam.tar.gz",	optional: true,	emit: guppy_fast5mod_bam_gz
	path "batch_${fast5_dir.baseName}.guppy.fast5mod_guppy2sam.bam*",	emit: guppy_fast5mod_bam
	path "batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv.gz",	emit: guppy_gcf52ref_tsv

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
			--verbose_logs  &>> Guppy.run.log
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		guppy_basecaller --input_path ${fast5_dir} \
			--save_path ${fast5_dir.baseName}.methcalled \
			--config ${params.GUPPY_METHCALL_MODEL} \
			--num_callers \$(( numProcessor )) \
			--fast5_out --compress_fastq\
			--verbose_logs \
			--device auto  &>> Guppy.run.log
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
	extract_methylation_fast5_support_dir.py \
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


// Combine Guppy runs' all results together
process GuppyComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}.guppy.*.combine.tsv.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}_raw_outputs/guppy",
		mode: "copy", pattern: "${params.dsname}.guppy_fast5mod.combined.bam.tar.gz",
		enabled: params.outputIntermediate

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x
	path y
	path reference_genome

	output:
	path "${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz", emit: guppy_fast5mod_combine_out
	path "${params.dsname}.guppy.gcf52ref_per_read.combine.tsv.gz", emit: guppy_gcf52ref_combine_out
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz", emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz", emit: site_unify
	path "${params.dsname}.guppy_fast5mod.combined.bam.tar.gz", emit: guppy_combine_raw_out_ch,  optional: true

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
	if [[ "${params.type}" == "human" ]] ; then
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
	elif [[ "${params.type}" == "ecoli" ]] ; then
		echo "### For ecoli, chr=${chrSet}"
		fast5mod call total.meth.bam ${referenceGenome} \
			meth.chr_${chrSet}.tsv \
			--meth cpg --quiet \
			--regions ${chrSet}
		gzip -f  meth.chr_${chrSet}.tsv && \
			mv meth.chr_${chrSet}.tsv.gz ${params.dsname}.guppy.fast5mod_per_site.combine.tsv.gz
	fi

	## Unify format output for read level
	unify_format_for_calls.sh \
		${params.dsname}  Guppy.gcf52ref ${params.dsname}.guppy.gcf52ref_per_read.combine.tsv.gz \
		.  \$((numProcessor))  1  ${chrSet}

	## Unify format output for site level
	unify_format_for_calls.sh \
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


// Tombo runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Tombo {
	tag "${resquiggleDir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_raw_outputs/tombo",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path resquiggleDir
	each path(reference_genome)


	output:
	path "batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats.bed.gz",	emit: tombo_out

	when:
	params.runMethcall && params.runTombo

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
	tombo_extract_per_read_stats.py \
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


// Combine Tombo runs' all results together
process TomboComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.*.combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x

	output:
	path "${params.dsname}.tombo.per_read.combine.bed.gz",	emit: tombo_combine_out
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz",	emit: site_unify

	when:
	x.size() >= 1 && params.runCombine

	"""
	touch ${params.dsname}.tombo.per_read.combine.bed.gz
	cat ${x} > ${params.dsname}.tombo.per_read.combine.bed.gz

	## Unify format output
	unify_format_for_calls.sh \
		${params.dsname}  Tombo ${params.dsname}.tombo.per_read.combine.bed.gz \
		.  \$((numProcessor))  12 ${chrSet}
	echo "### Tombo combine DONE"
	"""
}


// DeepMod runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepMod {
	tag "${basecallDir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_raw_outputs/deepmod",
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

	if [[ "${params.type}" = "human" ]] ; then
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


// Combine DeepMod runs' all results together
process DpmodComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}.deepmod.*.combine.bed.gz",
		enabled: params.outputRaw
	publishDir "${params.outdir}/${params.dsname}_raw_outputs/deepmod",
		mode: "copy", pattern: "${params.dsname}.deepmod.sum_chrs_mod.C.bed.tar.gz",
		enabled: params.outputIntermediate
	publishDir "${params.outdir}/${params.dsname}_raw_outputs/deepmod",
		mode: "copy", pattern: "${params.dsname}.deepmod_clusterCpG.all_chrs.C.bed.tar.gz",
		enabled: params.outputIntermediate
	publishDir "${params.outdir}/${params.dsname}_raw_outputs/deepmod",
		mode: "copy", pattern: "${params.dsname}.deepmod.all_batch.C.bed.tar.gz",
		enabled: params.outputIntermediate

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path x
	path deepmod_c_tar_file

	output:
	path "${params.dsname}.deepmod.*.combine.bed.gz",	emit: deepmod_combine_out
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz",	emit: site_unify
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
	sum_chr_mod.py \
		indir/ C ${params.dsname}.deepmod ${chrSet}  &>> DpmodComb.run.log

	> ${params.dsname}.deepmod.C_per_site.combine.bed

	## Note: for ecoli data, no pattern for chr*, but N*
	for f in \$(ls -1 indir/${params.dsname}.deepmod.*.C.bed)
	do
	  cat \$f >> ${params.dsname}.deepmod.C_per_site.combine.bed
	done
	gzip -f ${params.dsname}.deepmod.C_per_site.combine.bed

	if [[ "${params.type}" == "human" && "${isDeepModCluster}" == "true" ]] ; then
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
		hm_cluster_predict.py \
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
	unify_format_for_calls.sh \
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


// Read level unified output, and get METEORE output
process METEORE {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz"

	input:
	path naonopolish
	path megalodon
	path deepsignal

	output:
	path "${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz",	emit: meteore_combine_out, optional: true
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score.tsv.gz",	emit: read_unify, 	optional: true
	path "Site_Level-${params.dsname}/*-perSite-cov1.sort.bed.gz",	emit: site_unify,	optional: true

	when:
	params.runMethcall && params.runMETEORE

	"""
	## METEORE outputs by combining other tools
	outFileName=${params.dsname}_Megalodon_DeepSignal_combine.model_content.tsv
	> \$outFileName
	printf '%s\t%s\n' deepsignal ${deepsignal} >> \$outFileName
	printf '%s\t%s\n' megalodon ${megalodon} >> \$outFileName

	if [ -d /data ]; then
		METEOREDIR=\$(find /data -maxdepth 1 -name "${params.METEORE_Dir}" -type d)
	else
		METEOREDIR=""
	fi
	if [[ \$METEOREDIR == "" ]]; then
		wget ${params.METEOREGithub}  --no-verbose
		tar -xzf v1.0.0.tar.gz
		METEOREDIR=\$(find . -name "${params.METEORE_Dir}" -type d)
	fi

	## Degrade sk-learn for METEORE program if needed, it's model load need lower version
	## pip install -U scikit-learn==0.21.3
	combineScript=combination_model_prediction.py

	modelContentFileName=\$(find . -name "${params.dsname}_Megalodon_DeepSignal_combine.model_content.tsv" -type f)
	# Use the optimized model (n_estimator = 3 and max_dep = 10)
	# Please note this optimized model is reported in METEORE paper, ref: https://github.com/comprna/METEORE#command
	# paper: https://doi.org/10.1101/2020.10.14.340315, text:  random forest (RF) (parameters: max_depth=3 and n_estimator=10)
	# therefore, we output optimized model results
	\${combineScript}\
		-i \${modelContentFileName} \
		-m optimized -b \${METEOREDIR} \
		-o ${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz\
		&>> METEORE.run.log

	# Use default parameters from the sklearn library (n_estimator = 100 and max_dep = None)
	\${combineScript}\
		-i \${modelContentFileName} -m default -b \${METEOREDIR} \
		-o ${params.dsname}.meteore.megalodon_deepsignal_default_rf_model_per_read.combine.tsv.gz\
		&>> METEORE.run.log

	# Read level and site level output
	if [ -f ${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz ] ; then
		mkdir -p Read_Level-${params.dsname}
		zcat ${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz | \
			awk -F '\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$6,\$5}' |
			gzip -f > Read_Level-${params.dsname}/${params.dsname}_METEORE-perRead-score.tsv.gz

		## Unify format output for site level
		unify_format_for_calls.sh \
			${params.dsname}  METEORE\
			${params.dsname}.meteore.megalodon_deepsignal_optimized_rf_model_per_read.combine.tsv.gz \
			.  \$((numProcessor))  2  ${chrSet}\
			&>> METEORE.run.log
	fi
	echo "### METEORE consensus DONE"
	"""
}


// Not cache due to the script contains run information, each time of resume run will need updated
process Report {
	tag "${params.dsname}"

	publishDir "${params.outdir}",
		mode: "copy", pattern: "README_${params.dsname}.txt"

	publishDir "${params.outdir}",
		mode: "copy", pattern: "${params.dsname}_nanome_report.html"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Site_Level-${params.dsname}",
		mode: "copy", pattern: "${params.dsname}_NANOME-*.sort.bed.gz"

	publishDir "${params.outdir}/MultiQC",
		mode: "copy", pattern: "multiqc_report.html"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "GenomeBrowser-${params.dsname}"

	input:
	path fileList
	path qc_report
	path src
	path utils
	path reference_genome

	output:
	path "${params.dsname}_nanome_report.html",	emit:	report_out_ch
	path "README_${params.dsname}.txt",	emit: 	readme_out_ch
	path "${params.dsname}_NANOME-*.sort.bed.gz",	emit: 	nanome_consensus_ch, optional: true
	path "multiqc_report.html",	emit: 	lbt_report_ch
	path "GenomeBrowser-${params.dsname}", emit:  genome_browser_ch, optional: true

	when:
	fileList.size() >= 1

	"""
	## NANOME consensus method
	if [[ ${params.nanomeNanopolish} == true ]]; then
		NanopolishSiteReport=\$(find . -maxdepth 1 -name '*Nanopolish-perSite-*.sort.bed.gz')
	fi
	if [[ ${params.nanomeMegalodon} == true ]]; then
		MegalodonSiteReport=\$(find . -maxdepth 1 -name '*Megalodon-perSite-*.sort.bed.gz')
	fi
	if [[ ${params.nanomeDeepSignal} == true ]]; then
		DeepSignalSiteReport=\$(find . -maxdepth 1 -name '*DeepSignal-perSite-*.sort.bed.gz')
	fi
	if [[ ${params.nanomeGuppy} == true ]]; then
		GuppySiteReport=\$(find . -maxdepth 1 -name '*Guppy-perSite-*.sort.bed.gz')
	fi
	nanome_consensus.py\
	 	--site-reports   \${NanopolishSiteReport:-} \${MegalodonSiteReport:-}\
	 		\${DeepSignalSiteReport:-} \${GuppySiteReport:-}\
	 	--union -o ${params.dsname}_NANOME-perSite-cov1.sort.bed.gz &>> Report.run.log  || true

	##nanome_consensus.py\
	## 	--site-reports   \${NanopolishSiteReport:-} \${MegalodonSiteReport:-}\
	## 		\${DeepSignalSiteReport:-} \${GuppySiteReport:-}\
	## 	--join -o ${params.dsname}_NANOMEJoin-perSite-cov1.sort.bed.gz  &>> Report.run.log  || true

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
	printf '%s\t%s\n' outputs ${params.outdir} >> running_information.tsv

	## Note that the reason of report process can not be cached, is due to
	## Above script codes will be changed each time, so report can not apply old cached script

	## Get basecalling results from NanoComp
	basecallOutputFile=\$(find ${params.dsname}_QCReport/ -name "*NanoStats.txt" -type f)

	## Generate report dir and html utilities
	if [ -d /opt/nanome ]; then
		nanome_dir=/opt/nanome
	else
		nanome_dir="."
	fi
	mkdir -p ${params.dsname}_NANOME_report
	cp \${nanome_dir}/src/nanocompare/report/style.css ${params.dsname}_NANOME_report/
	cp -rf \${nanome_dir}/src/nanocompare/report/icons ${params.dsname}_NANOME_report/
	cp -rf \${nanome_dir}/src/nanocompare/report/js ${params.dsname}_NANOME_report/

	## Generate html report
	gen_html_report.py \
		${params.dsname} \
		running_information.tsv \
		\${basecallOutputFile} \
		. \
		${params.dsname}_NANOME_report \
		./src/nanocompare/report  &>> Report.run.log

	## Combine a single html report
	## No tty usage, ref: https://github.com/remy/inliner/issues/151
	script -qec "inliner ${params.dsname}_NANOME_report/nanome_report.html" /dev/null  > ${params.dsname}_nanome_report.html

	## Test on lifebit only
	cp ${params.dsname}_nanome_report.html   multiqc_report.html

	PYTHONIOENCODING=UTF-8 gen_readme.py\
		utils/readme.txt.template ${params.dsname} ${params.outdir}\
		${workflow.projectDir} ${workflow.workDir} "${workflow.commandLine}"\
		${workflow.runName} "${workflow.start}"\
		> README_${params.dsname}.txt   2>> Report.run.log

	if [[ ${params.outputGenomeBrowser} == true ]]; then
		mkdir -p GenomeBrowser-${params.dsname}
		find . -name '*-perSite-cov1.sort.bed.gz' | \
			parallel -j\$((numProcessor)) -v \
				"basefn={/}  && \
					zcat {} | \
					awk '{printf \\"%s\\t%d\\t%d\\t%2.3f\\n\\" , \\\$1,\\\$2,\\\$3,\\\$7}' > \
					GenomeBrowser-${params.dsname}/\\\${basefn/-perSite-cov1.sort.bed.gz/.bedgraph} && \
					LC_COLLATE=C sort -u -k1,1 -k2,2n \
						GenomeBrowser-${params.dsname}/\\\${basefn/-perSite-cov1.sort.bed.gz/.bedgraph} > \
							GenomeBrowser-${params.dsname}/\\\${basefn/-perSite-cov1.sort.bed.gz/.sorted.bedgraph} && \
					bedGraphToBigWig GenomeBrowser-${params.dsname}/\\\${basefn/-perSite-cov1.sort.bed.gz/.sorted.bedgraph} \
						reference_genome/chrom.sizes   GenomeBrowser-${params.dsname}/\\\${basefn/-perSite-cov1.sort.bed.gz/.bw} && \
						rm -f GenomeBrowser-${params.dsname}/\\\${basefn/-perSite-cov1.sort.bed.gz/.bedgraph}"
	fi
	echo "### report html DONE"
	"""
}


workflow {
	genome_ch = Channel.fromPath(genome_path, type: 'any', checkIfExists: true)
	megalodon_model_ch = Channel.fromPath(megalodon_model_tar, type: 'any', checkIfExists: true)

	EnvCheck(genome_ch, megalodon_model_ch)
	Untar(fast5_tar_ch)
	if (params.runBasecall) {
		Basecall(Untar.out.untar)
		QCExport(Basecall.out.basecall.collect(), EnvCheck.out.reference_genome)
	}

	// Resquiggle running if use Tombo or DeepSignal
	if ((params.runDeepSignal || params.runTombo) && params.runMethcall) {
		Resquiggle(Basecall.out.basecall, EnvCheck.out.reference_genome)
	}

	if (params.runNanopolish && params.runMethcall) {
		Nanopolish(Basecall.out.basecall, EnvCheck.out.reference_genome)
		comb_nanopolish = NplshComb(Nanopolish.out.nanopolish_out.collect())
		s1 = comb_nanopolish.site_unify
		r1 = comb_nanopolish.read_unify
	} else {
		s1 = Channel.empty()
		r1 = Channel.empty()
	}

	if (params.runMegalodon && params.runMethcall) {
		Megalodon(Untar.out.untar, EnvCheck.out.reference_genome, EnvCheck.out.megalodon_model)
		comb_megalodon = MgldnComb(Megalodon.out.megalodon_out.collect())
		s2 = comb_megalodon.site_unify
		r2 = comb_megalodon.read_unify
	} else {
		s2 = Channel.empty()
		r2 = Channel.empty()
	}

	if (params.runDeepSignal && params.runMethcall) {
		DeepSignal(Resquiggle.out.resquiggle, EnvCheck.out.reference_genome)
		comb_deepsignal = DpSigComb(DeepSignal.out.deepsignal_out.collect())
		s3 = comb_deepsignal.site_unify
		r3 = comb_deepsignal.read_unify
	} else {
		s3 = Channel.empty()
		r3 = Channel.empty()
	}

	if (params.runGuppy && params.runMethcall) {
		Guppy(Untar.out.untar, EnvCheck.out.reference_genome, ch_utils)
		comb_guppy = GuppyComb(Guppy.out.guppy_fast5mod_bam.collect(), Guppy.out.guppy_gcf52ref_tsv.collect(), EnvCheck.out.reference_genome)
		s4 = comb_guppy.site_unify
		r4 = comb_guppy.read_unify
	} else {
		s4 = Channel.empty()
		r4 = Channel.empty()
	}

	if (params.runTombo && params.runMethcall) {
		Tombo(Resquiggle.out.resquiggle, EnvCheck.out.reference_genome)
		comb_tombo = TomboComb(Tombo.out.tombo_out.collect())
		s5 = comb_tombo.site_unify
		r5 = comb_tombo.read_unify
	} else {
		s5 = Channel.empty()
		r5 = Channel.empty()
	}

	if (params.runDeepMod && params.runMethcall) {
		DeepMod(Basecall.out.basecall, EnvCheck.out.reference_genome)
		comb_deepmod = DpmodComb(DeepMod.out.deepmod_out.collect(), Channel.fromPath(deepmod_tar_file))
		s6 = comb_deepmod.site_unify
	} else {
		s6 = Channel.empty()
	}

	if (params.runMETEORE && params.runMethcall) {
		// Read level combine a list for top3 used by METEORE
		METEORE(r1, r2, r3)
		s7 = METEORE.out.site_unify
		r7 = METEORE.out.read_unify
	} else {
		s7 = Channel.empty()
		r7 = Channel.empty()
	}

	// Site level combine a list
	Channel.fromPath("${projectDir}/README.md").concat(
		s1, s2, s3, s4, s5, s6, s7
		).toList().set { tools_site_unify }
	Report(tools_site_unify, QCExport.out.qc_report, ch_src, ch_utils, EnvCheck.out.reference_genome)
}
