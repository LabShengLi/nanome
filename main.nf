#!/usr/bin/env nextflow
/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : main.nf
 @Software : NANOME project
 @Organization : JAX Li Lab
----------------------------------------------------------------------------------------
*/
// We now support both latest and lower versions, due to Lifebit CloudOS is only support 20.04
// Note: NXF_VER=20.04.1 nextflow run main.nf -profile test,singularity
if( nextflow.version.matches(">= 20.07.1") ){
	nextflow.enable.dsl = 2
} else {
	// Support lower version of nextflow
	nextflow.preview.dsl = 2
}

include {helpMessage} from './modules/HELP'

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Check mandatory params
if (! params.dsname)  exit 1, "Missing --dsname option for dataset name, check command help use --help"
if (! params.input)  exit 1, "Missing --input option for input data, check command help use --help"
//if ( !file(params.input.toString()).exists() )   exit 1, "input does not exist, check params: --input ${params.input}"

// Parse genome params
gbl_genome_map = params.genome_map

gbl_genome_path = gbl_genome_map[params.genome] ? gbl_genome_map[params.genome] : params.genome

// infer dataType, chrSet based on reference genome name, hg - human, ecoli - ecoli, otherwise is other reference genome
humanChrSet = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY'

genome_basefn = (new File(params.genome)).name
if (genome_basefn.startsWith('hg') || genome_basefn.startsWith('chm13') ||
		(params.dataType && params.dataType == 'human')) {
	dataType = params.dataType ? params.dataType : "human"
	// default for human chr
	chrSet = params.chrSet ? params.chrSet : humanChrSet
} else if (genome_basefn.startsWith('mm') || (params.dataType && params.dataType == 'mouse') ){
	dataType = params.dataType ? params.dataType : "mouse"
	// default for mouse chr
	chrSet = params.chrSet ? params.chrSet : humanChrSet
} else if (genome_basefn.startsWith('ecoli') || (params.dataType && params.dataType == 'ecoli')) {
	dataType = params.dataType ? params.dataType : "ecoli"
	chrSet = params.chrSet ? params.chrSet : 'NC_000913.3'
} else {
	// if not infer data type, use other
	dataType = params.dataType ? params.dataType : "other"

	if (params.chrSet)  chrSet = params.chrSet
	else exit 1, "Missing --chrSet option for other reference genome, please specify chromosomes used in reference genome [${params.genome}]"
}

// chrSet1 and dataType1 is the infered params, defined from chrSet and dataType (not in scope of params), will be used in every modules
params.chrSet1 = chrSet
params.dataType1 = dataType

// Get src and utils dir
projectDir = workflow.projectDir
ch_utils = Channel.fromPath("${projectDir}/utils",  type: 'dir', followLinks: false)
ch_src   = Channel.fromPath("${projectDir}/src",  type: 'dir', followLinks: false)

// Reference genome, chom size file name, will be used in every modules
params.referenceGenome = "${params.GENOME_DIR}/${params.GENOME_FN}"
params.chromSizesFile = "${params.GENOME_DIR}/${params.CHROM_SIZE_FN}"

if (dataType == 'human') { isDeepModCluster = params.useDeepModCluster }
else { isDeepModCluster = false }
params.isDeepModCluster = isDeepModCluster


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
		.set{ ch_inputs }
} else if (params.input.contains('*') || params.input.contains('?')) {
	// match all files in the folder, note: input must use quote string '', prevent expand in advance
	// such as --input '/fastscratch/liuya/nanome/NA12878/NA12878_CHR22/input_chr22/*'
	Channel.fromPath(params.input, type: 'any', checkIfExists: true)
		.set{ ch_inputs }
} else {
	// For single file/wildcard matched files
	Channel.fromPath( params.input, checkIfExists: true ).set{ ch_inputs }
}

// Header log info
def summary = [:]
summary['dsname'] 			= params.dsname
summary['input'] 			= params.input

if (gbl_genome_map[params.genome]) { summary['genome'] = "${params.genome} - [${gbl_genome_path}]" }
else { summary['genome'] = params.genome }

summary['\nRunning settings']         = "--------"
summary['processors'] 		= params.processors
summary['chrSet'] 			= chrSet   // .split(' ').join(',')
summary['dataType'] 		= dataType

if (params.runBasecall) summary['runBasecall'] = 'Yes'
if (params.skipBasecall) summary['skipBasecall'] = 'Yes'
if (params.runResquiggle) summary['runResquiggle'] = 'Yes'

if (params.runMethcall) {
	if (params.runNanopolish) summary['runNanopolish'] = 'Yes'
	if (params.runMegalodon) summary['runMegalodon'] = 'Yes'
	if (params.runDeepSignal2) summary['runDeepSignal2'] = 'Yes'
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

	if (params.runNewTool && params.newModuleConfigs)
		summary['runNewTool'] = params.newModuleConfigs.collect{it.name}.join(',')
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
if (params.multi_to_single_fast5) { summary['multi_to_single_fast5'] 	= params.multi_to_single_fast5 }
if (params.phasing) { summary['phasing'] 	= params.phasing }
if (params.hmc) { summary['hmc'] 	= params.hmc }
if (params.ctg_name) { summary['ctg_name'] 	= params.ctg_name }


summary['\nModel summary']         = "--------"
if (params.runBasecall && !params.skipBasecall) summary['GUPPY_BASECALL_MODEL'] 	= params.GUPPY_BASECALL_MODEL
if (params.runMethcall && params.runMegalodon)
	summary['MEGALODON_MODEL'] 	= params.rerio? 'Rerio:' + params.MEGALODON_MODEL : 'Remora:' + params.remoraModel
if (params.runMethcall && params.runDeepSignal) summary['DEEPSIGNAL_MODEL_DIR/DEEPSIGNAL_MODEL'] =\
 	params.DEEPSIGNAL_MODEL_DIR + "/" + params.DEEPSIGNAL_MODEL
if (params.runMethcall && params.runGuppy) summary['GUPPY_METHCALL_MODEL'] 	= params.GUPPY_METHCALL_MODEL
if (params.runMethcall && params.runDeepMod) {
	if (isDeepModCluster) {
		summary['DEEPMOD_RNN_MODEL;DEEPMOD_CLUSTER_MODEL'] = \
			"${params.DEEPMOD_RNN_MODEL};${params.DEEPMOD_CLUSTER_MODEL}"
		summary['DEEPMOD_CFILE'] = params.DEEPMOD_CFILE
	} else {
		summary['DEEPMOD_RNN_MODEL'] = "${params.DEEPMOD_RNN_MODEL}"
	}
}
if (params.runNANOME) {
	summary['NANOME_MODEL'] = "${params.NANOME_MODEL}"
	summary['CS_MODEL_FILE'] = "${params.CS_MODEL_FILE}"
	// summary['CS_MODEL_SPEC'] = "${params.CS_MODEL_SPEC}"
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
if (params.echo)  		summary['echo'] = params.echo
if (params.cleanup)   	summary['cleanup'] = params.cleanup

if (workflow.profile.contains('hpc') || workflow.profile.contains('winter') ||\
 	workflow.profile.contains('sumner') ) {
	summary['\nHPC settings']         = "--------"
    summary['queue']        = params.queue
    summary['qos']          = params.qos
    summary['memory']       = params.memory
    summary['time']         = params.time
    summary['queueSize']    = params.queueSize
    if (params.gresOptions) {summary['gresOptions'] = params.gresOptions }
}
if (workflow.profile.contains('google') || (params.config && params.config.contains('lifebit'))) {
	summary['\nGCP settings']         = "--------"
	if (params.projectCloud) {
		summary['projectCloud']    = params.projectCloud
	}
	if (params.config) { // lifebit specific settings
		summary['config']       		= params.config
	}
	summary['networkCloud']       = params.networkCloud
	summary['subnetworkCloud']	= params.subnetworkCloud

    summary['locationCloud']          = params.locationCloud
    summary['regionCloud']            = params.regionCloud
    summary['zoneCloud']       		= params.zoneCloud

    summary['bootDiskSizeCloud']       = params.bootDiskSizeCloud

	if (params.machineType)   	summary['machineType'] = params.machineType
	else {
		summary['machineType:cpus']         	= params.processors
		summary['machineType:memory']         	= params.memory
	}
	summary['gpuType']         	= params.gpuType
	summary['gpuNumber']        = params.gpuNumber

	// summary['lowDiskSize']      = params.lowDiskSize
	summary['midDiskSize']      = params.midDiskSize
	summary['highDiskSize']     = params.highDiskSize
}

log.info """\
NANOME - NF PIPELINE (v$workflow.manifest.version)
by Li Lab at The Jackson Laboratory
https://github.com/LabShengLi/nanome
================================="""
.stripIndent()

log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "================================="


include { ENVCHECK } from './modules/ENVCHECK'  // addParams(chrSet1: "${chrSet}", dataType1:"${dataType}")

include { UNTAR } from './modules/UNTAR'

include { BASECALL } from './modules/BASECALL'

include { ALIGNMENT } from './modules/ALIGNMENT'

include { QCEXPORT } from './modules/QCEXPORT'

include { RESQUIGGLE } from './modules/RESQUIGGLE'

include { NANOPOLISH; NPLSHCOMB } from './modules/NANOPOLISH'

include { MEGALODON; MGLDNCOMB } from './modules/MEGALODON'

include { DEEPSIGNAL; DPSIGCOMB } from './modules/DEEPSIGNAL'

include { DEEPSIGNAL2; DEEPSIGNAL2COMB } from './modules/DEEPSIGNAL2'

include { Guppy; GuppyComb; Tombo; TomboComb; DeepMod; DpmodComb; METEORE } from './modules/OLDTOOLS'

include { Guppy6; Guppy6Comb } from './modules/GUPPY6'

include { NewTool; NewToolComb } from './modules/NEWTOOLS'

include { CLAIR3; PHASING } from './modules/PHASING'

include { CONSENSUS } from './modules/CONSENSUS'

include { EVAL } from './modules/EVAL'

include { REPORT } from './modules/REPORT'

// place holder channel, used for empty file of a channel
null1 = Channel.fromPath("${projectDir}/utils/null1")
null2 = Channel.fromPath("${projectDir}/utils/null2")
null3 = Channel.fromPath("${projectDir}/utils/null3")

workflow {
	if ( !file(gbl_genome_path.toString()).exists() )
		exit 1, "genome reference path does not exist, check params: --genome ${params.genome}"

	ch_genome = Channel.fromPath(gbl_genome_path, type: 'any', checkIfExists: true)

	// rerio model dir will be download in ENVCHECK if needed
	ch_rerio_dir = (params.rerio && params.rerioDir) ? Channel.fromPath(params.rerioDir, type: 'any', checkIfExists: true) :
							null1

	// deepsignal model dir will be downloaded in ENVCHECK if needed
	if (params.runDeepSignal) {
		ch_deepsignal_dir = params.deepsignalDir ?
				Channel.fromPath(params.deepsignalDir, type: 'any', checkIfExists: true) :
				Channel.fromPath(params.DEEPSIGNAL_MODEL_ONLINE, type: 'any', checkIfExists: true)
	} else {
		// use null placeholder
		ch_deepsignal_dir = null2
	}

	ENVCHECK(ch_genome, ch_utils, ch_rerio_dir, ch_deepsignal_dir)
	UNTAR(ch_inputs)

	if (params.runBasecall) {
		BASECALL(UNTAR.out.untar)
		ALIGNMENT(BASECALL.out.basecall, ENVCHECK.out.reference_genome)
		QCEXPORT(BASECALL.out.basecall.collect(),
					ALIGNMENT.out.alignment.collect(),
					ENVCHECK.out.reference_genome)
	}

	// Resquiggle running if use Tombo or DeepSignal
	if (((params.runDeepSignal || params.runTombo || params.runDeepSignal2) && params.runMethcall)
		|| params.runResquiggle) {
		resquiggle = RESQUIGGLE(UNTAR.out.untar_tuple.join(BASECALL.out.basecall_tuple), ENVCHECK.out.reference_genome)
		f1 = params.feature_extract ? resquiggle.feature_extract : Channel.empty()
	} else {
		f1 = Channel.empty()
	}

	if (params.runNanopolish && params.runMethcall) {
		NANOPOLISH(UNTAR.out.untar_tuple.join(BASECALL.out.basecall_tuple).join(ALIGNMENT.out.alignment_tuple),
			ENVCHECK.out.reference_genome)
		comb_nanopolish = NPLSHCOMB(NANOPOLISH.out.nanopolish_tsv.collect(), ch_src, ch_utils)
		s1 = comb_nanopolish.site_unify
		r1 = comb_nanopolish.read_unify
	} else {
		s1 = Channel.empty()
		r1 = Channel.empty()
	}

	if (params.runMegalodon && params.runMethcall) {
		MEGALODON(UNTAR.out.untar, ENVCHECK.out.reference_genome, ENVCHECK.out.rerio)
		comb_megalodon = MGLDNCOMB(MEGALODON.out.megalodon_tsv.collect(),
							MEGALODON.out.megalodon_mod_mappings.collect(),
							ch_src, ch_utils)
		s2 = comb_megalodon.site_unify
		r2 = comb_megalodon.read_unify
	} else {
		s2 = Channel.empty()
		r2 = Channel.empty()
	}

	if (params.runDeepSignal && params.runMethcall) {
		DEEPSIGNAL(RESQUIGGLE.out.resquiggle, ENVCHECK.out.reference_genome,
					ENVCHECK.out.deepsignal_model)
		comb_deepsignal = DPSIGCOMB(DEEPSIGNAL.out.deepsignal_tsv.collect(), ch_src, ch_utils)
		s3 = comb_deepsignal.site_unify
		r3 = comb_deepsignal.read_unify
	} else {
		s3 = Channel.empty()
		r3 = Channel.empty()
	}

	if (params.runDeepSignal2 && params.runMethcall) {
		deepsignal2_model_file = Channel.fromPath(params.DEEPSIGNAL2_MODEL_FILE, type: 'any', checkIfExists: true)
		deepsignal2 = DEEPSIGNAL2(RESQUIGGLE.out.resquiggle,
					ENVCHECK.out.reference_genome,
					ch_src, ch_utils, deepsignal2_model_file)
		comb_deepsignal2 = DEEPSIGNAL2COMB(DEEPSIGNAL2.out.deepsignal2_batch_per_read.collect(),
						DEEPSIGNAL2.out.deepsignal2_batch_feature.collect(),
						ch_src, ch_utils
						)
		f2 = Channel.empty() //deepsignal2.deepsignal2_feature_out
		s3_1 = comb_deepsignal2.site_unify
		r3_1 = comb_deepsignal2.read_unify
	} else {
		f2 = Channel.empty()
		s3_1 = Channel.empty()
		r3_1 = Channel.empty()
	}

	if (params.runGuppy && params.runMethcall) {
		Guppy6(UNTAR.out.untar, ENVCHECK.out.reference_genome, ch_utils)

		comb_guppy6 = Guppy6Comb(Guppy6.out.guppy_batch_bam_out.collect(),
								ENVCHECK.out.reference_genome,
								ch_src, ch_utils)

		s4 = Channel.empty()
		r4 = Channel.empty()
//		s4 = comb_guppy.site_unify
//		r4 = comb_guppy.read_unify
	} else {
		s4 = Channel.empty()
		r4 = Channel.empty()
	}

	if (params.runTombo && params.runMethcall) {
		Tombo(RESQUIGGLE.out.resquiggle, ENVCHECK.out.reference_genome)
		comb_tombo = TomboComb(Tombo.out.tombo_tsv.collect(), ch_src, ch_utils)
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
			if ( !file(params.DEEPMOD_CFILE.toString()).exists() )
				exit 1, "DEEPMOD_CFILE does not exist, check params: --DEEPMOD_CFILE ${params.DEEPMOD_CFILE}"
			ch_ctar = Channel.fromPath(params.DEEPMOD_CFILE, type:'any', checkIfExists: true)
		}
		DeepMod(BASECALL.out.basecall, ENVCHECK.out.reference_genome)
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
			if ( !file(params.METEOREDir.toString()).exists() )
				exit 1, "METEOREDir does not exist, check params: --METEOREDir ${params.METEOREDir}"
			METEOREDir_ch = Channel.fromPath(params.METEOREDir, type: 'any', checkIfExists: true)
		}
		METEORE(r1, r2, r3, ch_src, ch_utils, METEOREDir_ch)
		s7 = METEORE.out.site_unify
		r7 = METEORE.out.read_unify
	} else {
		s7 = Channel.empty()
		r7 = Channel.empty()
	}

	if (params.runNewTool && params.newModuleConfigs) {
		newModuleCh = Channel.of( params.newModuleConfigs ).flatten()
		// ref: https://www.nextflow.io/docs/latest/operator.html#combine
		NewTool(newModuleCh.combine(BASECALL.out.basecall), ENVCHECK.out.reference_genome, params.referenceGenome)
		NewToolComb(NewTool.out.batch_out.collect(), newModuleCh, ch_src)

		s_new = NewToolComb.out.site_unify
		r_new = NewToolComb.out.read_unify
	} else {
		s_new = Channel.empty()
		r_new = Channel.empty()
	}

	null2.concat(
		r1, r2, r3, r3_1, f1, f2
		).toList().set { top3_tools_read_unify }

	if (params.runNANOME) {
		consensus = CONSENSUS(top3_tools_read_unify, ch_src, ch_utils)
		s8 = consensus.site_unify
		r8 = consensus.read_unify
	} else {
		s8 = Channel.empty()
		r8 = Channel.empty()
	}

	null2.concat(
		r1, r2, r3, r8, f1, f2
		).toList().set { tools_read_unify }

	// perform evaluation of tools' methylation results
	if (params.runEval) {
		bg1 = params.bg1 ? Channel.fromPath(params.bg1) : Channel.empty()
		bg2 = params.bg2 ? Channel.fromPath(params.bg2) : Channel.empty()

		null1.concat(
			bg1, bg2
		).toList().set { bg_list }

		if (params.genome_annotation_dir) {
			genome_annotation_ch = Channel.fromPath(params.genome_annotation_dir)
		} else {
			genome_annotation_ch = null3
		}

		EVAL(tools_read_unify, bg_list, ch_src, ch_utils, genome_annotation_ch)
	}

	// Site level combine a list
	null1.concat(
		s1, s2, s3, s3_1, s4, s5, s6, s7, s_new, s8
		).toList().set { tools_site_unify }

	REPORT(tools_site_unify, top3_tools_read_unify,
			ENVCHECK.out.tools_version_tsv, QCEXPORT.out.qc_report,
			ENVCHECK.out.reference_genome, ch_src, ch_utils)

	if (params.phasing) {
		CLAIR3(QCEXPORT.out.bam_data, ENVCHECK.out.reference_genome)
		null1.concat(
			MGLDNCOMB.out.megalodon_combine,
			MGLDNCOMB.out.read_unify,
			CONSENSUS.out.nanome_combine_out,
			CONSENSUS.out.read_unify,
			NPLSHCOMB.out.nanopolish_combine_out_ch
			).toList().set { mega_and_nanome_ch }
		PHASING(mega_and_nanome_ch, CLAIR3.out.clair3_out_ch,
				ch_src, QCEXPORT.out.bam_data, ENVCHECK.out.reference_genome)
	}
}
