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
if( nextflow.version.matches(">= 20.07.1") ){
	nextflow.enable.dsl=2
} else {
	// Support lower version of nextflow
	nextflow.preview.dsl=2
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
	  --dsname		Dataset/analysis name
	  --input		Input path for raw fast5 files (folders, tar/tar.gz files)
	  --genome		Genome reference name ('hg38', 'ecoli', or 'hg38_chr22') or a directory, the directory must contain only one .fasta file with .fasta.fai index file. Default is hg38

	General options:
	  --processors		Processors used for each task
	  --outdir		Output dir, default is 'results'
	  --chrSet		Chromosomes used in analysis, default is chr1-22, X and Y, for human. For E. coli data, it is default as 'NC_000913.3'. For other reference genome, please specify each chromosome with space seperated.
	  --cleanAnalyses	If clean old basecalling info in fast5 files
	  --skipBasecall	Skip redo basecalling if users provide basecalled inputs

	  --cleanup		If clean work dir after complete, default is false

	Tools specific options:
	  --run[Tool-name]	By default, we run top four performers in nanome paper, specify '--run[Tool-name]' can include other tool, supported tools: NANOME, Megalodon, Nanopolish, DeepSignal, Guppy, Tombo, METEORE, and DeepMod
	  --rerioDir		Rerio dir for Megalodon model, default will get online
	  --MEGALODON_MODEL	Megalodon model name, default is 'res_dna_r941_min_modbases_5mC_v001.cfg'
	  --guppyDir		Guppy installation local directory, used only for conda environment
	  --GUPPY_BASECALL_MODEL	Guppy basecalling model, default is 'dna_r9.4.1_450bps_hac.cfg'
	  --GUPPY_METHCALL_MODEL	Guppy methylation calling model, default is 'dna_r9.4.1_450bps_modbases_5mc_hac.cfg'
	  --deepsignalDir	DeepSignal model dir, default will get online
	  --tomboResquiggleOptions	Tombo resquiggle options for super long/damaged sequencing, set to '--signal-length-range 0 500000  --sequence-length-range 0 50000'
	  --moveOption	If using move table for DeepMod, default is true
	  --useDeepModCluster	If using DeepMod cluster model for human, default is false
	  --METEOREDir	METEORE model dir, default will get online

	Running environment options:
	  --docker_name		Docker name used for pipeline, default is 'liuyangzzu/nanome:latest'
	  --singularity_name	Singularity name used for pipeline, default is 'docker://liuyangzzu/nanome:latest'
	  --singularity_cache	Singularity cache dir, default is 'local_singularity_cache'
	  --conda_name		Conda name used for pipeline, default is 'nanome'
	  --conda_base_dir	Conda base directory, default is '/opt/conda'

	Platform specific options:
	  --queue		SLURM job submission queue name for cluster running, e.g., 'gpu'
	  --qos			SLURM job submission qos name for cluster running, e.g., 'inference'
	  --gresOptions		SLURM job submission GPU allocation options for cluster running, e.g., 'gpu:v100:1'
	  --time		SLURM job submission time allocation options for cluster running, e.g., '2h', '1d'
	  --memory		SLURM job submission memory allocation options for cluster running, e.g., '32GB'

	  --googleProjectName	Google Cloud Platform (GCP) project name for google-lifesciences
	  --config		Lifebit CloudOS config file, e.g., 'conf/executors/lifebit.config'

	-profile options:
	  Use this parameter to choose a predefined configuration profile. Profiles can give configuration presets for different compute environments.

	  test		A bundle of input params for ecoli test
	  test_human	A bundle of input params for human test
	  docker 	A generic configuration profile to be used with Docker, pulls software from Docker Hub: liuyangzzu/nanome:latest
	  singulairy	A generic configuration profile to be used with Singularity, pulls software from: docker://liuyangzzu/nanome:latest
	  conda		Please only use conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity. Check our GitHub for how to install local conda enviroment
	  hpc		A generic configuration profile to be used on HPC cluster with SLURM
	  google	A generic configuration profile to be used on Google Cloud platform with 'google-lifesciences'

	Contact to https://github.com/TheJacksonLaboratory/nanome/issues for bug report.
	""".stripIndent()
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Check mandatory params
if (! params.dsname)  exit 1, "Missing --dsname option for dataset name, check command help use --help"
if (! params.input)  exit 1, "Missing --input option for input data, check command help use --help"
if ( !file(params.input.toString()).exists() )   exit 1, "input does not exist, check params: --input ${params.input}"

// Parse genome params
genome_map = params.genome_map

if (genome_map[params.genome]) { genome_path = genome_map[params.genome] } else { 	genome_path = params.genome }

// infer dataType, chrSet based on reference genome name, hg - human, ecoli - ecoli, otherwise is other reference genome
if (params.genome.contains('hg') || (params.dataType && params.dataType == 'human')) {
	dataType = "human"
	if (!params.chrSet) {
		// default for human, if false or 'false' (string), using '  '
		chrSet = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'
	} else {
		chrSet = params.chrSet
	}
} else if (params.genome.contains('ecoli') || (params.dataType && params.dataType == 'ecoli')) {
	dataType = "ecoli"
	if (!params.chrSet) {
		// default for ecoli
		chrSet = 'NC_000913.3'
	} else {
		chrSet = params.chrSet
	}
} else {
	// default will not found name, use other
	if (!params.dataType) { dataType = 'other' } else { dataType = params.dataType }
	if (!params.chrSet) {
		// No default value for other reference genome
		exit 1, "Missing --chrSet option for other reference genome, please sepecify chromsomes used in reference genome [${params.genome}]"
	}
	chrSet = params.chrSet
}


// Get src and utils dir
projectDir = workflow.projectDir
ch_utils = Channel.fromPath("${projectDir}/utils",  type: 'dir', followLinks: false)
ch_src   = Channel.fromPath("${projectDir}/src",  type: 'dir', followLinks: false)

// Reference genome, deepmod cluster settings
def referenceGenome = 'reference_genome/ref.fasta'
def chromSizesFile = 'reference_genome/chrom.sizes'

if (dataType == 'human') { isDeepModCluster = params.useDeepModCluster } else { 	isDeepModCluster = false }


// Collect all folders of fast5 files, and send into Channels for pipelines
if (params.input.endsWith(".filelist.txt")) {
	// list of files in filelist.txt
	Channel.fromPath( params.input, checkIfExists: true )
		.splitCsv(header: false)
		.map {
			if (!file(it[0]).exists())  {
				log.warn "File not exists: ${it[0]}, check file list: ${params.input}"
			} else {
				return file(it[0])
			}
		}
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

if (genome_map[params.genome] != null) { summary['genome'] = "${params.genome} - [${genome_path}]" }
else { summary['genome'] = params.genome }

summary['\nRunning settings']         = "--------"
summary['processors'] 		= params.processors
summary['chrSet'] 			= chrSet.split(' ').join(',')
summary['dataType'] 		= dataType

if (params.runBasecall) summary['runBasecall'] = 'Yes'
if (params.skipBasecall) summary['skipBasecall'] = 'Yes'

if (params.runMethcall) {
	if (params.runNanopolish) summary['runNanopolish'] = 'Yes'
	if (params.runMegalodon) summary['runMegalodon'] = 'Yes'
	if (params.runDeepSignal) summary['runDeepSignal'] = 'Yes'
	if (params.runGuppy) summary['runGuppy'] = 'Yes'
	if (params.runTombo) summary['runTombo'] = 'Yes'
	if (params.runMETEORE) summary['runMETEORE'] = 'Yes'
	if (params.runDeepMod) summary['runDeepMod'] = 'Yes'

	if (params.runDeepMod) {
		summary['runDeepMod'] = 'Yes'
		if (params.moveOption)  summary['runDeepMod'] = summary['runDeepMod'] + ' + (move table)'
		if (isDeepModCluster)  {
			summary['runDeepMod'] = summary['runDeepMod'] + ' + (cluster model)'
		}
	}

	if (params.runNANOME) summary['runNANOME'] = 'Yes'
}

if (params.cleanAnalyses) summary['cleanAnalyses'] = 'Yes'
if (params.deepsignalDir) { summary['deepsignalDir'] = params.deepsignalDir }
if (params.rerioDir) {
	summary['rerioDir'] = params.rerioDir
	summary['MEGALODON_MODEL'] = params.MEGALODON_MODEL
}
if (params.METEOREDir) { summary['METEOREDir'] = params.METEOREDir }
if (params.guppyDir) { summary['guppyDir'] 	= params.guppyDir }
if (params.tomboResquiggleOptions) { summary['tomboResquiggleOptions'] 	= params.tomboResquiggleOptions }

if (params.outputBam) { summary['outputBam'] 	= params.outputBam }
if (params.outputONTCoverage) { summary['outputONTCoverage'] 	= params.outputONTCoverage }
if (params.outputIntermediate) { summary['outputIntermediate'] 	= params.outputIntermediate }
if (params.outputRaw) { summary['outputRaw'] 	= params.outputRaw }
if (params.outputGenomeBrowser) { summary['outputGenomeBrowser'] 	= params.outputGenomeBrowser }
if (params.deduplicate) { summary['deduplicate'] 	= params.deduplicate }
if (params.sort) { summary['sort'] 	= params.sort }

summary['\nModel summary']         = "--------"
if (params.runBasecall && !params.skipBasecall) summary['GUPPY_BASECALL_MODEL'] 	= params.GUPPY_BASECALL_MODEL
if (params.runMethcall && params.runMegalodon) summary['MEGALODON_MODEL'] 	= params.MEGALODON_MODEL
if (params.runMethcall && params.runDeepSignal) summary['DEEPSIGNAL_MODEL_DIR/DEEPSIGNAL_MODEL'] 	= params.DEEPSIGNAL_MODEL_DIR + "/" + params.DEEPSIGNAL_MODEL
if (params.runMethcall && params.runGuppy) summary['GUPPY_METHCALL_MODEL'] 	= params.GUPPY_METHCALL_MODEL
if (params.runMethcall && params.runDeepMod) {
	if (isDeepModCluster) {
		summary['DEEPMOD_RNN_MODEL;DEEPMOD_CLUSTER_MODEL'] = "${params.DEEPMOD_RNN_MODEL};${params.DEEPMOD_CLUSTER_MODEL}"
		summary['DEEPMOD_CFILE'] = params.DEEPMOD_CFILE
	} else {
		summary['DEEPMOD_RNN_MODEL'] = "${params.DEEPMOD_RNN_MODEL}"
	}
}

summary['\nPipeline settings']         = "--------"
summary['Working dir'] 		= workflow.workDir
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Profile']          = workflow.profile
summary['Config files'] 	= workflow.configFiles.join(',')
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['errorStrategy']    = params.errorStrategy
summary['maxRetries']       = params.maxRetries
if (params.echo)  summary['echo'] = params.echo
if (params.cleanup)   summary['cleanup'] = params.cleanup

if (workflow.profile.contains('hpc') || workflow.profile.contains('winter') || workflow.profile.contains('sumner') ) {
	summary['\nHPC settings']         = "--------"
    summary['queue']        = params.queue
    summary['qos']          = params.qos
    summary['memory']       = params.memory
    summary['time']         = params.time
    if (params.gresOptions) {summary['gresOptions'] = params.gresOptions }
}
if (workflow.profile.contains('google') || (params.config && params.config.contains('lifebit'))) {
	summary['\nGCP settings']         = "--------"
	if (!params.config || !params.config.contains('lifebit')) {
		summary['googleProjectName']    = params.googleProjectName
	} else { // lifebit specific settings
		summary['config']       		= params.config
		summary['networkLifebit']       = params.networkLifebit
		summary['subnetworkLifebit']	= params.subnetworkLifebit
		summary['zoneCloud']       		= params.zoneCloud
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

	script:
	"""
	date; hostname; pwd
	echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"

	## Validate nanome container/environment is correct
	bash utils/validate_nanome_container.sh  tools_version_table.tsv

	## Untar and prepare megalodon model
	if [[ ${params.runMegalodon} == true && ${params.runMethcall} == true ]]; then
		if [[ ${rerioDir} == null* ]] ; then
			# Obtain and run R9.4.1, MinION, 5mC CpG model from Rerio
			git clone ${params.rerioGithub}
			rerio/download_model.py rerio/basecall_models/${params.MEGALODON_MODEL.replace('.cfg', '')}
		else
			if [[ ${rerioDir} != rerio ]] ; then
				## rename it to rerio for output channel
				cp  -a ${rerioDir}  rerio
			fi
		fi

		ls -lh rerio/
	fi

	## Untar and prepare deepsignal model
	if [[ ${params.runDeepSignal} == true && ${params.runMethcall} == true ]]; then
		if [[ ${deepsignalDir} == *.tar.gz ]] ; then
			## Get DeepSignal Model online
			tar -xzf ${deepsignalDir}
		else
		 	if [[ ${deepsignalDir} != ${params.DEEPSIGNAL_MODEL_DIR} ]] ; then
				## rename it to deepsignal default dir name
				cp  -a ${deepsignalDir}  ${params.DEEPSIGNAL_MODEL_DIR}
			fi
		fi

		## Check DeepSignal model
		ls -lh ${params.DEEPSIGNAL_MODEL_DIR}/
	fi

	if [[ ${params.runBasecall} == true || ${params.runMethcall} == true ]]; then
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

		ls -lh reference_genome/
	fi

	echo "### Check reference genome and chrSet"
	echo "referenceGenome=${referenceGenome}"
	echo "chromSizesFile=${chromSizesFile}"
	echo "chrSet=[${chrSet}]"
	echo "dataType=${dataType}"
	echo "cpus=$task.cpus"
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

	script:
	cores = task.cpus * params.highProcTimes
	if (!params.skipBasecall) { // perform basecall
		"""
		date; hostname; pwd
		echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}"

		## Extract input files tar/tar.gz/folder
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
			find ${fast5_tar}/ -name '*.fast5' | \
				parallel -j$cores  cp {} untarTempDir/
		else
			echo "### Untar error for input=${fast5_tar}"
		fi

		## Move fast5 raw/basecalled files into XXX.untar folder
		mkdir -p ${fast5_tar.baseName}.untar

		find untarTempDir -name "*.fast5" -type f | \
			parallel -j$cores  mv {}  ${fast5_tar.baseName}.untar/

		## Clean temp files
		rm -rf untarTempDir

		## Clean old basecalled analyses in input fast5 files
		if [[ "${params.cleanAnalyses}" == true ]] ; then
			echo "### Start cleaning old analysis"
			## python -c 'import h5py; print(h5py.version.info)'
			clean_old_basecall_in_fast5.py \
				-i ${fast5_tar.baseName}.untar --is-indir --verbose\
				--processor $cores
		fi

		totalFiles=\$( find ${fast5_tar.baseName}.untar -name "*.fast5" -type f | wc -l )
		echo "### Total fast5 input files:\${totalFiles}"
		if (( totalFiles==0 )); then
			echo "### no fast5 files at ${fast5_tar.baseName}.untar, skip this job"
			rm -rf ${fast5_tar.baseName}.untar
		fi
		echo "### Untar DONE"
		"""
	} else {
		"""
		date; hostname; pwd

		## Extract input files tar/tar.gz/folder
		infn="${fast5_tar}"
		mkdir -p untarTempDir
		if [ "\${infn##*.}" == "tar" ]; then
			### deal with tar
			tar -xf \${infn} -C untarTempDir
		elif [ "\${infn##*.}" == "gz" ]; then
			### deal with tar.gz
			tar -xzf \${infn} -C untarTempDir
		elif [[ -d ${fast5_tar} ]]; then
			## user provide basecalled input dir, just cp them
			mkdir -p untarTempDir/test
			cp -rf ${fast5_tar}/*   untarTempDir/test/
		else
			echo "### Untar error for input=${fast5_tar}"
		fi

		## Move fast5 raw/basecalled files into XXX.untar folder
		mkdir -p ${fast5_tar.baseName}.untar
		## Keep the directory structure for basecalled input
		mv untarTempDir/*/*   ${fast5_tar.baseName}.untar/

		## Clean temp files
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

	if [[ ${params.skipBasecall} == false ]] ; then
		if [[ \${commandType} == "cpu" ]]; then
			## CPU version command
			guppy_basecaller --input_path ${fast5_dir} \
				--save_path "${fast5_dir.baseName}.basecalled" \
				--config ${params.GUPPY_BASECALL_MODEL} \
				--num_callers ${task.cpus} \
				--fast5_out --compress_fastq\
				--verbose_logs  &>> Basecall.run.log
		elif [[ \${commandType} == "gpu" ]]; then
			## GPU version command
			guppy_basecaller --input_path ${fast5_dir} \
				--save_path "${fast5_dir.baseName}.basecalled" \
				--config ${params.GUPPY_BASECALL_MODEL} \
				--num_callers ${task.cpus} \
				--fast5_out --compress_fastq\
				--verbose_logs \
				-x auto  &>> Basecall.run.log
		else
			echo "### error value for commandType=\${commandType}"
			exit 255
		fi
	else
		## Just use user's basecalled input
		cp -rf ${fast5_dir}/*   ${fast5_dir.baseName}.basecalled/
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
	 	parallel -j${task.cpus * params.highProcTimes} 'rm -f {}'

	## After basecall, rename and publish summary filenames, summary may also be used by resquiggle
	mv ${fast5_dir.baseName}.basecalled/sequencing_summary.txt \
		${fast5_dir.baseName}.basecalled/${fast5_dir.baseName}-sequencing_summary.txt

    ## Clean
    if [[ ${params.cleanStep} == "true" ]]; then
    	echo "### No need to clean"
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
	path "${params.dsname}_basecall_report.html",	optional: true, emit: qc_html
	path "${params.dsname}_QCReport",				emit: qc_report
	path "${params.dsname}_bam_data",		optional: true,	 emit: bam_data

	script:
	cores = task.cpus * params.highProcTimes
	samtools_cores = task.cpus * params.mediumProcTimes
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

	mkdir -p ${params.dsname}_QCReport
	if [[ ${params.skipQC} == false ]]; then
		## Perform QC report by NanoComp
		NanoComp --summary ${params.dsname}_combine_sequencing_summary.txt.gz  \
			--names ${params.dsname} --outdir ${params.dsname}_QCReport -t $cores \
			--raw  -f pdf -p ${params.dsname}_   &>> QCReport.run.log
	fi

	if [[ ${params.outputBam} == true  || ${params.outputONTCoverage} == true ]]; then
		## Combine all batch fq.gz
		> merge_all_fq.fq.gz
		cat *.basecalled/batch_basecall_combine_fq_*.fq.gz > merge_all_fq.fq.gz
		echo "### Fastq merge from all batches done!"

		## After basecall, we align results to merged, sorted bam, can be for ONT coverage analyses/output bam
		# align FASTQ files to reference genome, write sorted alignments to a BAM file
		minimap2 -t ${samtools_cores} -a  -x map-ont \
			${referenceGenome} \
			merge_all_fq.fq.gz | \
			samtools sort -@ ${samtools_cores} -T tmp -o \
				${params.dsname}_merge_all_bam.bam &&\
			samtools index -@ ${samtools_cores}  ${params.dsname}_merge_all_bam.bam
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

	[ -f ${params.dsname}_combine_sequencing_summary.txt.gz ] && \
		mv -f ${params.dsname}_combine_sequencing_summary.txt.gz ${params.dsname}_QCReport/
	[ -f ${params.dsname}_QCReport/${params.dsname}_NanoComp-report.html ] && \
		mv -f ${params.dsname}_QCReport/${params.dsname}_NanoComp-report.html ${params.dsname}_basecall_report.html

	## Clean
	if [[ ${params.cleanStep} == "true" ]]; then
		rm -f ${params.dsname}.coverage.positivestrand.bed.gz ${params.dsname}.coverage.negativestrand.bed.gz
		rm -f merge_all_fq.fq.gz
		if [[ ${params.outputBam} == false ]]; then
			rm -f ${params.dsname}_merge_all_bam.bam*
		else
			mkdir -p ${params.dsname}_bam_data
			mv  ${params.dsname}_merge_all_bam.bam*  ${params.dsname}_bam_data/
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

	script:
	cores = task.cpus * params.highProcTimes
	samtools_cores = task.cpus * params.mediumProcTimes
	resquiggle_cores = (task.cpus*params.reduceProcTimes).intValue()
	"""
	### copy basecall workspace files, due to tombo resquiggle modify base folder
	rm -rf ${basecallIndir.baseName}.resquiggle
	mkdir -p ${basecallIndir.baseName}.resquiggle/workspace

	### original basecalled results will be parrallelly used by other processes
	cp -f ${basecallIndir}/batch_basecall_combine_fq_*.fq.gz  ${basecallIndir.baseName}.resquiggle/

	## cp -rf ${basecallIndir}/workspace  ${basecallIndir.baseName}.resquiggle/
	find ${basecallIndir}/workspace -name '*.fast5' -type f| \
		parallel -j${task.cpus * params.highProcTimes}  \
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
		--overwrite --processes  ${samtools_cores} \
		&>> Resquiggle.run.log
	echo "### tombo preprocess DONE"

	### Need to check Tombo resquiggle bugs, lots of users report long runtime and hang at nearly completion for large data
	### ref: https://github.com/nanoporetech/tombo/issues/139, https://github.com/nanoporetech/tombo/issues/111
	### ref: https://github.com/nanoporetech/tombo/issues/365, https://github.com/nanoporetech/tombo/issues/167
	### ref: https://nanoporetech.github.io/tombo/resquiggle.html?highlight=processes
	### Out of memory solution for large data: --tomboResquiggleOptions '--signal-length-range 0 500000  --sequence-length-range 0 50000'
	tombo resquiggle\
		--processes ${resquiggle_cores} \
		--corrected-group ${params.ResquiggleCorrectedGroup} \
		--basecall-group ${params.BasecallGroupName} \
		--basecall-subgroup ${params.BasecallSubGroupName}\
		--ignore-read-locks ${params.tomboResquiggleOptions ? params.tomboResquiggleOptions : ''}\
		--overwrite \
		${basecallIndir.baseName}.resquiggle/workspace \
		${referenceGenome} &>> Resquiggle.run.log

	echo "### tombo resquiggle DONE"
	"""
}


// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Nanopolish {
	tag "${basecallDir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/nanopolish",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path basecallDir
	each path(reference_genome)

	output:
	path "${params.dsname}_nanopolish_batch_${basecallDir.baseName}.*.gz", 	emit: nanopolish_out

	when:
	params.runMethcall && params.runNanopolish

	script:
	samtools_cores = task.cpus * params.mediumProcTimes
	nanopolish_cores = (task.cpus*params.reduceProcTimes).intValue()

	"""
	## Put all fq and bam files into working dir, DO NOT affect the basecall dir
	bamFileName="${params.dsname}.batch_${basecallDir.baseName}.sorted.bam"

	## Do alignment firstly, find the combined fastq file
	fastqFile=\$(find ${basecallDir}/ -name 'batch_basecall_combine_fq_*.fq.gz' -type f)

	## Index, ref: https://github.com/jts/nanopolish#data-preprocessing
	## Index the raw read with fastq, we do not index in basecalled dir, in case of cache can be work
	ln -s \${fastqFile}  \${fastqFile##*/}
	nanopolish index -d ${basecallDir}/workspace -s ${basecallDir}/${basecallDir.baseName}-sequencing_summary.txt  \${fastqFile##*/}

	## Aligning reads to the reference genome, ref: https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html#aligning-reads-to-the-reference-genome
	minimap2 -t ${samtools_cores}  -a -x map-ont ${referenceGenome} \${fastqFile##*/} | \
		samtools sort -@ ${samtools_cores} -T tmp -o \${bamFileName} &&\
		samtools index -@ ${samtools_cores}  \${bamFileName}
	echo "### Alignment step: minimap2 and samtools DONE"

	## Calling methylation, ref: https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html#calling-methylation
	## there are segment fault issues, if set -t to a large number or use low memory,
	## ref: https://github.com/jts/nanopolish/issues/872
	## ref: https://github.com/jts/nanopolish/issues/683, https://github.com/jts/nanopolish/issues/580
	nanopolish call-methylation \
		-t ${nanopolish_cores}\
	 	-r \${fastqFile##*/} \
		-b \${bamFileName} -g ${referenceGenome} -q cpg | \
		awk 'NR>1' | \
		gzip -f > ${params.dsname}_nanopolish_batch_${basecallDir.baseName}.tsv.gz

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
		pattern: "${params.dsname}_nanopolish_per_read_combine.*.gz",
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
	path "${params.dsname}_nanopolish_per_read_combine.tsv.gz",	emit: nanopolish_combine_out_ch
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify

	when:
	x.size() >= 1 && params.runCombine

	"""
	> ${params.dsname}_nanopolish_per_read_combine.tsv.gz
	cat ${x} > ${params.dsname}_nanopolish_per_read_combine.tsv.gz

	if [[ ${params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		zcat ${params.dsname}_nanopolish_per_read_combine.tsv.gz |\
			sort -V -u -k1,1 -k3,3n -k4,4n -k5,5 -k2,2 |\
			gzip -f > ${params.dsname}_nanopolish_per_read_combine.sort.tsv.gz
		rm ${params.dsname}_nanopolish_per_read_combine.tsv.gz &&\
			mv ${params.dsname}_nanopolish_per_read_combine.sort.tsv.gz ${params.dsname}_nanopolish_per_read_combine.tsv.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  Nanopolish Nanopolish \
		${params.dsname}_nanopolish_per_read_combine.tsv.gz \
		.  ${task.cpus}  12 ${params.sort  ? true : false}   "${chrSet}"

	echo "### Nanopolish combine DONE"
	"""
}


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${fast5_dir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/megalodon",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path fast5_dir
	each path(reference_genome)
	each path(rerio_dir)

	output:
	path "${params.dsname}_megalodon_batch_${fast5_dir.baseName}.*.gz", emit: megalodon_out

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
	else
		echo "Detect GPU, using GPU commandType"
		commandType='gpu'
	fi

	echo "### Guppy dir:"
	which guppy_basecall_server

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		## Ref: https://github.com/nanoporetech/megalodon
		## CPU issues: https://github.com/nanoporetech/megalodon/issues/172
		megalodon \
			${fast5_dir}   --overwrite \
			--mod-motif m CG 0 \
			--outputs per_read_mods mods per_read_refs \
			--mod-output-formats bedmethyl wiggle \
			--write-mods-text --write-mod-log-probs \
			--guppy-server-path \$(which guppy_basecall_server) \
			--guppy-config ${params.MEGALODON_MODEL} \
			--guppy-params "-d ./${rerio_dir}/basecall_models/" \
			--guppy-timeout ${params.GUPPY_TIMEOUT} \
			--reference ${referenceGenome} \
			--processes $cores\
			&>> Megalodon.run.log
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		## Ref: https://github.com/nanoporetech/megalodon
		## --guppy-params "-d megalodon_model_dir/ --num_callers 8 --ipc_threads  42"
		megalodon \
			${fast5_dir}   --overwrite \
			--mod-motif m CG 0 \
			--outputs per_read_mods mods per_read_refs \
			--mod-output-formats bedmethyl wiggle \
			--write-mods-text --write-mod-log-probs \
			--guppy-server-path \$(which guppy_basecall_server) \
			--guppy-config ${params.MEGALODON_MODEL} \
			--guppy-params "-d ./${rerio_dir}/basecall_models/" \
			--guppy-timeout ${params.GUPPY_TIMEOUT} \
			--reference ${referenceGenome} \
			--processes $cores --devices 0\
			&>> Megalodon.run.log
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi

	awk 'NR>1' megalodon_results/per_read_modified_base_calls.txt | gzip -f > \
		${params.dsname}_megalodon_batch_${fast5_dir.baseName}.tsv.gz

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
process MgldnComb {
	tag "${params.dsname}"

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
	path x
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_megalodon_per_read_combine.*.gz",	emit: megalodon_combine
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify

	when:
	x.size() >= 1  && params.runCombine

	"""
	> ${params.dsname}_megalodon_per_read_combine.bed.gz
	cat ${x} > ${params.dsname}_megalodon_per_read_combine.bed.gz

	if [[ ${params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		## sort order: Chr, Start, (End), ID, Strand
		zcat ${params.dsname}_megalodon_per_read_combine.bed.gz |\
			sort -V -u -k2,2 -k4,4n -k1,1 -k3,3 |\
			gzip -f > ${params.dsname}_megalodon_per_read_combine.sort.bed.gz
		rm ${params.dsname}_megalodon_per_read_combine.bed.gz &&\
			mv ${params.dsname}_megalodon_per_read_combine.sort.bed.gz  ${params.dsname}_megalodon_per_read_combine.bed.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  Megalodon Megalodon \
		${params.dsname}_megalodon_per_read_combine.bed.gz \
		.  ${task.cpus}  12  ${params.sort  ? true : false}  "${chrSet}"

	echo "### Megalodon combine DONE"
	"""
}


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${indir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/deepsignal",
		mode: "copy",
		enabled: params.outputIntermediate

	input:
	path indir
	each path(reference_genome)
	each path(deepsignal_model_dir)

	output:
	path "${params.dsname}_deepsignal_batch_${indir.baseName}.*.gz",	emit: deepsignal_out

	when:
	params.runMethcall && params.runDeepSignal

	script:
	cores = task.cpus * params.highProcTimes
	"""
	DeepSignalModelBaseDir="."
	commandType='gpu'
	outFile="${params.dsname}_deepsignal_batch_${indir.baseName}.tsv"

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path "\${DeepSignalModelBaseDir}/${params.DEEPSIGNAL_MODEL_DIR}/${params.DEEPSIGNAL_MODEL}" \
			--result_file \${outFile} \
			--reference_path ${referenceGenome} \
			--corrected_group ${params.ResquiggleCorrectedGroup} \
			--nproc $cores \
			--is_gpu no
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		deepsignal call_mods \
			--input_path ${indir}/workspace \
			--model_path "\${DeepSignalModelBaseDir}/${params.DEEPSIGNAL_MODEL_DIR}/${params.DEEPSIGNAL_MODEL}" \
			--result_file \${outFile} \
			--reference_path ${referenceGenome} \
			--corrected_group ${params.ResquiggleCorrectedGroup} \
			--nproc $cores \
			--is_gpu yes
	else
		echo "### error value for commandType=\${commandType}"
		exit 255
	fi

	gzip -f \${outFile}
	echo "### DeepSignal methylation DONE"
	"""
}


// Combine DeepSignal runs' all results together
process DpSigComb {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_deepsignal_per_read_combine.*.gz",
		enabled: params.outputRaw

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "Site_Level-${params.dsname}/*-perSite-cov1*.gz"

	input:
	path x
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_deepsignal_per_read_combine.*.gz",	emit: deepsignal_combine_out
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify

	when:
	x.size() >= 1 && params.runCombine

	"""
	touch ${params.dsname}_deepsignal_per_read_combine.tsv.gz
	cat ${x} > ${params.dsname}_deepsignal_per_read_combine.tsv.gz

	if [[ ${params.deduplicate} == true ]] ; then
		echo "### Deduplicate for read-level outputs"
		## sort order: Chr, Start, (End), ID, Strand
		zcat ${params.dsname}_deepsignal_per_read_combine.tsv.gz |\
			sort -V -u -k1,1 -k2,2n -k5,5 -k3,3 |\
			gzip -f > ${params.dsname}_deepsignal_per_read_combine.sort.tsv.gz
		rm ${params.dsname}_deepsignal_per_read_combine.tsv.gz &&\
			mv ${params.dsname}_deepsignal_per_read_combine.sort.tsv.gz  ${params.dsname}_deepsignal_per_read_combine.tsv.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  DeepSignal DeepSignal\
		${params.dsname}_deepsignal_per_read_combine.tsv.gz \
		.  $task.cpus  12 ${params.sort  ? true : false}  "${chrSet}"
	echo "### DeepSignal combine DONE"
	"""
}


// methylation calling for Guppy
process Guppy {
	tag "${fast5_dir.baseName}"

	publishDir "${params.outdir}/${params.dsname}_intermediate/guppy",
		mode: "copy",
		pattern: "bamfile_${params.dsname}_guppy_fast5mod_batch_${fast5_dir.baseName}_guppy2sam.bam.tar.gz",
		enabled: params.outputIntermediate

	publishDir "${params.outdir}/${params.dsname}_intermediate/guppy",
		mode: "copy",
		pattern: "${params.dsname}_guppy_gcf52ref_batch_${fast5_dir.baseName}.*.gz",
		enabled: params.outputIntermediate && params.runGuppyGcf52ref

	input:
	path fast5_dir
	each path(reference_genome)
	each path(utils)

	output:
	path "${params.dsname}_guppy_fast5mod_batch_${fast5_dir.baseName}_guppy2sam.bam*",	emit: guppy_fast5mod_bam
	path "${params.dsname}_guppy_gcf52ref_batch_${fast5_dir.baseName}.*.gz",	optional: true, emit: guppy_gcf52ref_tsv

	path "bamfile_${params.dsname}_guppy_fast5mod_batch_${fast5_dir.baseName}_guppy2sam.bam.tar.gz",	optional: true,	emit: guppy_fast5mod_bam_gz

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
	else
		echo "Detect GPU, using GPU commandType"
		commandType='gpu'
	fi

	if [[ ${params.skipBasecall} == false ]]; then
		indir=${fast5_dir}
	else
		indir=${fast5_dir}/workspace
	fi

	mkdir -p ${fast5_dir.baseName}.methcalled

	if [[ \${commandType} == "cpu" ]]; then
		## CPU version command
		guppy_basecaller --input_path \${indir} --recursive\
			--save_path ${fast5_dir.baseName}.methcalled \
			--config ${params.GUPPY_METHCALL_MODEL} \
			--num_callers $task.cpus \
			--fast5_out --compress_fastq\
			--verbose_logs  &>> Guppy.run.log
	elif [[ \${commandType} == "gpu" ]]; then
		## GPU version command
		guppy_basecaller --input_path \${indir} --recursive\
			--save_path ${fast5_dir.baseName}.methcalled \
			--config ${params.GUPPY_METHCALL_MODEL} \
			--num_callers $task.cpus \
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
		parallel -j${task.cpus * params.mediumProcTimes} 'rm -f {}'

	if [[ ${params.runGuppyGcf52ref} == true ]]; then
		## gcf52ref ways
		minimap2 -t ${task.cpus * params.mediumProcTimes} -a -x map-ont ${referenceGenome} \
			batch_combine_fq.fq.gz | \
			samtools sort -@ ${samtools_cores} \
				-T tmp -o gcf52ref.batch.${fast5_dir.baseName}.bam &&\
			samtools index -@ ${samtools_cores} \
				gcf52ref.batch.${fast5_dir.baseName}.bam
		echo "### gcf52ref minimap2 alignment is done"

		## Modified version, support dir input, not all fast5 files (too long arguments)
		extract_methylation_fast5_support_dir.py \
			-p ${samtools_cores}  ${fast5_dir.baseName}.methcalled/workspace
		echo "### gcf52ref extract to db done"

		## gcf52ref files preparation
		### git clone https://github.com/kpalin/gcf52ref.git
		tar -xzf utils/gcf52ref.tar.gz -C .
		patch gcf52ref/scripts/extract_methylation_from_rocks.py < \
			utils/gcf52ref.patch

		python gcf52ref/scripts/extract_methylation_from_rocks.py \
			-d base_mods.rocksdb \
			-a gcf52ref.batch.${fast5_dir.baseName}.bam \
			-r ${referenceGenome} \
			-o tmp.batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv
		echo "### gcf52ref extract to tsv done"

		awk 'NR>1' tmp.batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv | gzip -f > \
			${params.dsname}_guppy_gcf52ref_batch_${fast5_dir.baseName}.tsv.gz &&\
			rm -f tmp.batch_${fast5_dir.baseName}.guppy.gcf52ref_per_read.tsv
		echo "### gcf52ref extraction DONE"
	fi

	## fast5mod ways
	FAST5PATH=${fast5_dir.baseName}.methcalled/workspace
	OUTBAM=${params.dsname}_guppy_fast5mod_batch_${fast5_dir.baseName}_guppy2sam.bam

	fast5mod guppy2sam \${FAST5PATH} --reference ${referenceGenome} \
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
		samtools merge -@${task.cpus * params.mediumProcTimes}  ${params.dsname}_guppy_fast5mod_combine.bam {}

	### sort is not needed due to merge the sorted bam, ref: http://www.htslib.org/doc/samtools-merge.html
	### samtools sort -@ ${samtools_cores} total.meth.bam
	samtools index -@ ${samtools_cores}  ${params.dsname}_guppy_fast5mod_combine.bam
	echo "Samtool merge and index for fast5mod DONE"

	if [[ ${params.outputIntermediate} == true ]] ; then
		tar -czf bamfile_${params.dsname}_guppy_fast5mod_combine.bam.tar.gz ${params.dsname}_guppy_fast5mod_combine.bam*
	fi

	awk '/^>/' ${referenceGenome} | awk '{print \$1}' \
		> rf_chr_all_list.txt
	if [[ "${dataType}" == "human" ]] ; then
		echo "### For human, extract chr1-22, X and Y"
		> chr_all_list.txt
		for chrname in ${chrSet}
		do
			if cat rf_chr_all_list.txt | grep -w ">\${chrname}" -q ; then
				echo \${chrname} >> chr_all_list.txt
			fi
		done
		rm  -f rf_chr_all_list.txt
		echo "### Chomosome list"
		cat chr_all_list.txt

		## Ref: https://github.com/nanoporetech/medaka/issues/177
		parallel -j$cores  -v \
			"fast5mod call  ${params.dsname}_guppy_fast5mod_combine.bam  ${referenceGenome} \
				meth.chr_{}.tsv  --meth cpg --quiet --regions {} ; \
				gzip -f meth.chr_{}.tsv" :::: chr_all_list.txt

		touch ${params.dsname}_guppy_fast5mod_per_site_combine.tsv.gz
		for chrname in ${chrSet}
		do
			if [ -f "meth.chr_\${chrname}.tsv.gz" ]; then
				cat  meth.chr_\${chrname}.tsv.gz >> \
					${params.dsname}_guppy_fast5mod_per_site_combine.tsv.gz
			fi
		done
	else
		echo "### For other genome, chrSet=[${chrSet}]"
		for chrname in ${chrSet} ;
		do
			fast5mod call ${params.dsname}_guppy_fast5mod_combine.bam  ${referenceGenome} \
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
			mv ${params.dsname}_guppy_gcf52ref_per_read_combine.sort.tsv.gz ${params.dsname}_guppy_gcf52ref_per_read_combine.tsv.gz
	fi

	if [[ ${params.runGuppyGcf52ref} == true ]] ; then
		## Unify format output for read level
		bash utils/unify_format_for_calls.sh \
			${params.dsname}  Guppy Guppy.gcf52ref\
			 ${params.dsname}_guppy_gcf52ref_per_read_combine.tsv.gz \
			.  ${task.cpus}  1  ${params.sort  ? true : false}  "${chrSet}"
	fi

	## Unify format output for site level
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  Guppy Guppy\
		${params.dsname}_guppy_fast5mod_per_site_combine.tsv.gz \
		.  ${task.cpus}  2  ${params.sort  ? true : false}  "${chrSet}"

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
	path "${params.dsname}_tombo_batch_${resquiggleDir.baseName}.*.gz",	emit: tombo_out, optional: true

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
		${resquiggleDir.baseName}.Tombo.run.log

	retry=1
	## while grep -q "BrokenPipeError:" ${resquiggleDir.baseName}.Tombo.run.log
	while ! tail -n 1 ${resquiggleDir.baseName}.Tombo.run.log |  grep -q "100%"
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

	if [ -f "${chromSizesFile}" ]; then
		ln -s  ${chromSizesFile}  genome.chome.sizes
	else
		cut -f1,2 ${referenceGenome}.fai > genome.chome.sizes
	fi

	tombo_extract_per_read_stats.py \
		genome.chome.sizes \
		"${params.dsname}_tombo_batch_${resquiggleDir.baseName}.CpG.tombo.per_read_stats" \
		"${params.dsname}_tombo_batch_${resquiggleDir.baseName}.bed"

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
			mv ${params.dsname}_tombo_per_read_combine.sort.bed.gz  ${params.dsname}_tombo_per_read_combine.bed.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  Tombo Tombo\
		${params.dsname}_tombo_per_read_combine.bed.gz \
		.  $task.cpus  12  ${params.sort  ? true : false}  "${chrSet}"
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

	if [[ "${dataType}" == "human" ]] ; then
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
			--threads ${task.cpus * params.mediumProcTimes} \
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
	sum_chr_mod.py \
		indir/ C ${params.dsname}.deepmod ${chrSet.split(' ').join(',')}  &>> DpmodComb.run.log

	> ${params.dsname}_deepmod_c_per_site_combine.bed

	## Note: for ecoli data, no pattern for chr*, but N*
	for f in \$(ls -1 indir/${params.dsname}.deepmod.*.C.bed)
	do
	  cat \$f >> ${params.dsname}_deepmod_c_per_site_combine.bed
	done
	gzip -f ${params.dsname}_deepmod_c_per_site_combine.bed

	if [[ "${dataType}" == "human" && "${isDeepModCluster}" == "true" ]] ; then
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
			\${DeepModTrainModelDir}/${params.DEEPMOD_CLUSTER_MODEL} &>> DpmodComb.run.log

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

	if [[ "${isDeepModCluster}" == "true" ]] ; then
		callfn=${params.dsname}_deepmod_clusterCpG_per_site_combine.bed.gz
	else
		callfn=${params.dsname}_deepmod_c_per_site_combine.bed.gz
	fi

	## Unify format output
	bash utils/unify_format_for_calls.sh \
		${params.dsname}  DeepMod DeepMod\
		\${callfn} \
		.  $task.cpus  2  ${params.sort ? true : false} "${chrSet}"  &>> DpmodComb.run.log

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
	## pip install -U scikit-learn==0.21.3
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
		&>> METEORE.run.log

	# Use default parameters from the sklearn library (n_estimator = 100 and max_dep = None)
	##\${combineScript}\
	##	-i \${modelContentFileName} -m default -b \${METEOREDIR} \
	##	-o ${params.dsname}.meteore.megalodon_deepsignal_default_rf_model_per_read.combine.tsv.gz\
	##	&>> METEORE.run.log

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
			.  $task.cpus  12   ${params.sort ? true : false}  "${chrSet}"\
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

	publishDir "${params.outdir}/MultiQC",
		mode: "copy", pattern: "multiqc_report.html"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "GenomeBrowser-${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Site_Level-${params.dsname}/*-perSite-cov*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_nanome_*_per_read_combine.*.gz",
		enabled: params.outputRaw

	input:
	path site_fileList
	path read_fileList
	path tools_version_tsv
	path qc_report
	path reference_genome
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_nanome_report.html",	emit:	report_out_ch
	path "README_${params.dsname}.txt",	emit: 	readme_out_ch
	path "multiqc_report.html",	emit: 	lbt_report_ch
	path "GenomeBrowser-${params.dsname}", emit:  genome_browser_ch, optional: true
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify, optional: true
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify, optional: true
	path "${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.*.gz", emit: nanome_combine_out, optional: true

	"""
	if [[ ${params.runNANOME} == true ]] ; then
		## NANOME XGBoost method
		modelContentTSVFileName=${params.dsname}_nanome_${params.NANOME_MODEL}_model_content.tsv
		> \$modelContentTSVFileName
		passModelTsv=false
		if [[ "${params.NANOME_CONSENSUS_TOOLS}" == *"Nanopolish"* ]]; then
			NanopolishReadReport=\$(find . -maxdepth 1 -name '*Nanopolish-perRead-score.tsv.gz')
			if [[ -z \$NanopolishReadReport ]] ; then
				echo "### Not found Nanopolish read-level outputs"
				NanopolishReadReport="None"
			else
				passModelTsv=true
			fi
			printf '%s\t%s\n' nanopolish \${NanopolishReadReport} >> \$modelContentTSVFileName
		fi

		if [[ "${params.NANOME_CONSENSUS_TOOLS}" == *"Megalodon"* ]]; then
			MegalodonReadReport=\$(find . -maxdepth 1 -name '*Megalodon-perRead-score.tsv.gz')
			if [[ -z \$MegalodonReadReport ]] ; then
				echo "### Not found Megalodon read-level outputs"
				MegalodonReadReport="None"
			else
				passModelTsv=true
			fi
			printf '%s\t%s\n' megalodon \${MegalodonReadReport} >> \$modelContentTSVFileName
		fi

		if [[ "${params.NANOME_CONSENSUS_TOOLS}" == *"DeepSignal"* ]]; then
			DeepSignalReadReport=\$(find . -maxdepth 1 -name '*DeepSignal-perRead-score.tsv.gz')
			if [[ -z \$DeepSignalReadReport ]] ; then
				echo "### Not found DeepSignal read-level outputs"
				DeepSignalReadReport="None"
			else
				passModelTsv=true
			fi
			printf '%s\t%s\n' deepsignal \${DeepSignalReadReport} >> \$modelContentTSVFileName
		fi

		if [[ "\$passModelTsv" == true ]] ; then
			## NANOME XGBoost model results, if there are model results exists
			echo "### NANOME XGBoost predictions"
			PYTHONPATH=src python src/nanocompare/xgboost/xgboost_predict.py \
				--contain-na --tsv-input\
				--dsname ${params.dsname} -i \${modelContentTSVFileName}\
				-m ${params.NANOME_MODEL}  -o ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz &>> Report.run.log  || true

			if [[ ${params.deduplicate} == true ]] ; then
				echo "### Deduplicate for read-level outputs"
				## sort order: Chr, Start, (End), ID, Strand
				zcat ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz |\
					sort -V -u -k2,2 -k3,3n -k1,1 -k4,4 |\
					gzip -f > ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.sort.tsv.gz
				rm ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz &&\
					mv ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.sort.tsv.gz\
						${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz
			fi

			## Unify format output
			echo "### NANOME read/site level results"
			bash utils/unify_format_for_calls.sh \
				${params.dsname}  NANOME NANOME\
				${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz \
				.  $task.cpus  12  ${params.sort ? true : false}  "${chrSet}"
			ln -s Site_Level-${params.dsname}/${params.dsname}_NANOME-perSite-cov1.sort.bed.gz\
				${params.dsname}_NANOME-perSite-cov1.sort.bed.gz
		fi
	fi

	## Generate NF pipeline running information tsv
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

	if [[ -z "\${basecallOutputFile}" ]] ; then
		basecallOutputFile=None
	fi

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

	## Generate html NANOME report
	PYTHONPATH=src python src/nanocompare/report/gen_html_report.py\
		${params.dsname} \
		running_information.tsv \
		\${basecallOutputFile} \
		. \
		${params.dsname}_NANOME_report \
		./src/nanocompare/report\
		${tools_version_tsv}  &>> Report.run.log

	## Combine a single html report
	## No tty usage, ref: https://github.com/remy/inliner/issues/151
	script -qec "inliner ${params.dsname}_NANOME_report/nanome_report.html" /dev/null  > ${params.dsname}_nanome_report.html

	## Used for lifebit rendering feature
	cp ${params.dsname}_nanome_report.html   multiqc_report.html

	## Generate readme.txt
	PYTHONPATH=src PYTHONIOENCODING=UTF-8 python src/nanocompare/report/gen_txt_readme.py\
		src/nanocompare/report/readme.txt.template ${params.dsname} ${params.outdir}\
		${workflow.projectDir} ${workflow.workDir} "${workflow.commandLine}"\
		${workflow.runName} "${workflow.start}"\
		> README_${params.dsname}.txt   2>> Report.run.log

	## Output BigWig format for IGV
	if [[ ${params.outputGenomeBrowser} == true ]] ; then
		if command -v bedGraphToBigWig  &> /dev/null ; then
			mkdir -p GenomeBrowser-${params.dsname}
			find . -maxdepth 1 -name '*-perSite-cov1.sort.bed.gz' -print0 | \
				while IFS= read -r -d '' infn ; do
					echo "### processing infn=\$infn"
					## methfreq bw generation
					zcat \${infn} | awk '{printf "%s\\t%d\\t%d\\t%2.5f\\n" , \$1,\$2,\$3,\$7}' > \
						GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.bedgraph}
					LC_COLLATE=C sort -u -k1,1 -k2,2n \
						GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.bedgraph} > \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.sorted.bedgraph}
					bedGraphToBigWig GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.sorted.bedgraph} \
						reference_genome/chrom.sizes   GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.bw}
					rm -f GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.bedgraph}  \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.sorted.bedgraph}

					## coverage bw generation
					zcat \${infn} | \
						awk '{printf "%s\\t%d\\t%d\\t%d\\n" , \$1,\$2,\$3,\$8}' > \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.bedgraph}
					LC_COLLATE=C sort -u -k1,1 -k2,2n \
						GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.bedgraph} > \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.sorted.bedgraph}
					bedGraphToBigWig GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.sorted.bedgraph} \
						reference_genome/chrom.sizes   GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.bw}
					rm -f GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.bedgraph}  \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.sorted.bedgraph}
				done
		else
			echo "### ERROR: No bedGraphToBigWig in PATH, please install it"
		fi
	fi
	echo "### report html DONE"
	"""
}


workflow {
	if ( !file(genome_path.toString()).exists() )   exit 1, "genome reference file does not exist, check params: --genome ${params.genome}"
	genome_ch = Channel.fromPath(genome_path, type: 'any', checkIfExists: true)

	if (!params.rerioDir) { // default if null, will online downloading
		// This is only a place holder for input
		rerioDir = Channel.fromPath("${projectDir}/utils/null1", type: 'any', checkIfExists: false)
	} else {
		// User provide the dir
		if ( !file(params.rerioDir.toString()).exists() )   exit 1, "rerioDir does not exist, check params: --rerioDir ${params.rerioDir}"
		rerioDir = Channel.fromPath(params.rerioDir, type: 'any', checkIfExists: true)
	}

	if (!params.deepsignalDir) {
		// default if null, will online downloading
		deepsignalDir = Channel.fromPath(params.DEEPSIGNAL_MODEL_ONLINE, type: 'any', checkIfExists: true)
	} else {
		// User provide the dir
		if ( !file(params.deepsignalDir.toString()).exists() )   exit 1, "deepsignalDir does not exist, check params: --deepsignalDir ${params.deepsignalDir}"
		deepsignalDir = Channel.fromPath(params.deepsignalDir, type: 'any', checkIfExists: true)
	}

	EnvCheck(genome_ch, ch_utils, rerioDir, deepsignalDir)
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
		comb_nanopolish = NplshComb(Nanopolish.out.nanopolish_out.collect(), ch_src, ch_utils)
		s1 = comb_nanopolish.site_unify
		r1 = comb_nanopolish.read_unify
	} else {
		s1 = Channel.empty()
		r1 = Channel.empty()
	}

	if (params.runMegalodon && params.runMethcall) {
		Megalodon(Untar.out.untar, EnvCheck.out.reference_genome, EnvCheck.out.rerio)
		comb_megalodon = MgldnComb(Megalodon.out.megalodon_out.collect(), ch_src, ch_utils)
		s2 = comb_megalodon.site_unify
		r2 = comb_megalodon.read_unify
	} else {
		s2 = Channel.empty()
		r2 = Channel.empty()
	}

	if (params.runDeepSignal && params.runMethcall) {
		DeepSignal(Resquiggle.out.resquiggle, EnvCheck.out.reference_genome, EnvCheck.out.deepsignal_model)
		comb_deepsignal = DpSigComb(DeepSignal.out.deepsignal_out.collect(), ch_src, ch_utils)
		s3 = comb_deepsignal.site_unify
		r3 = comb_deepsignal.read_unify
	} else {
		s3 = Channel.empty()
		r3 = Channel.empty()
	}

	if (params.runGuppy && params.runMethcall) {
		Guppy(Untar.out.untar, EnvCheck.out.reference_genome, ch_utils)

		gcf52ref_ch = Channel.fromPath("${projectDir}/utils/null1").concat(Guppy.out.guppy_gcf52ref_tsv.collect())

		comb_guppy = GuppyComb(Guppy.out.guppy_fast5mod_bam.collect(),
								gcf52ref_ch,
								EnvCheck.out.reference_genome,
								ch_src, ch_utils)
		s4 = comb_guppy.site_unify
		r4 = comb_guppy.read_unify
	} else {
		s4 = Channel.empty()
		r4 = Channel.empty()
	}

	if (params.runTombo && params.runMethcall) {
		Tombo(Resquiggle.out.resquiggle, EnvCheck.out.reference_genome)
		comb_tombo = TomboComb(Tombo.out.tombo_out.collect(), ch_src, ch_utils)
		s5 = comb_tombo.site_unify
		r5 = comb_tombo.read_unify
	} else {
		s5 = Channel.empty()
		r5 = Channel.empty()
	}

	if (params.runDeepMod && params.runMethcall) {
		if (!isDeepModCluster) {
			// not use cluster model, only a place holder here
			ch_ctar = Channel.fromPath("${projectDir}/utils/null1", type:'any', checkIfExists: false)
		} else {
			if ( !file(params.DEEPMOD_CFILE.toString()).exists() )   exit 1, "DEEPMOD_CFILE does not exist, check params: --DEEPMOD_CFILE ${params.DEEPMOD_CFILE}"
			ch_ctar = Channel.fromPath(params.DEEPMOD_CFILE, type:'any', checkIfExists: true)
		}
		DeepMod(Basecall.out.basecall, EnvCheck.out.reference_genome)
		comb_deepmod = DpmodComb(DeepMod.out.deepmod_out.collect(), ch_ctar, ch_src, ch_utils)
		s6 = comb_deepmod.site_unify
	} else {
		s6 = Channel.empty()
	}

	if (params.runMETEORE && params.runMethcall) {
		// Read level combine a list for top3 used by METEORE
		if (!params.METEOREDir) {
			METEOREDir_ch = Channel.fromPath(params.METEORE_GITHUB_ONLINE, type: 'any', checkIfExists: true)
		} else {
			if ( !file(params.METEOREDir.toString()).exists() )   exit 1, "METEOREDir does not exist, check params: --METEOREDir ${params.METEOREDir}"
			METEOREDir_ch = Channel.fromPath(params.METEOREDir, type: 'any', checkIfExists: true)
		}
		METEORE(r1, r2, r3, ch_src, ch_utils, METEOREDir_ch)
		s7 = METEORE.out.site_unify
		r7 = METEORE.out.read_unify
	} else {
		s7 = Channel.empty()
		r7 = Channel.empty()
	}

	// Site level combine a list
	Channel.fromPath("${projectDir}/utils/null1").concat(
		s1, s2, s3, s4, s5, s6, s7
		).toList().set { tools_site_unify }

	Channel.fromPath("${projectDir}/utils/null2").concat(
		r1, r2, r3
		).toList().set { tools_read_unify }

	Report(tools_site_unify, tools_read_unify, EnvCheck.out.tools_version_tsv, QCExport.out.qc_report, EnvCheck.out.reference_genome, ch_src, ch_utils)
}
