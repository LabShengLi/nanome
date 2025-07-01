/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : HELP.nf
 @Software : NANOME project
 @Organization : Sheng Li Lab
----------------------------------------------------------------------------------------
*/
def helpMessage() {
	log.info"""
	NANOME - Nextflow PIPELINE (v$workflow.manifest.version)
	by Sheng Li Lab
	https://github.com/LabShengLi/nanome
	=================================
	Usage:
	The typical command is as follows:

	nextflow run LabShengLi/nanome -profile test,docker
	nextflow run LabShengLi/nanome -profile test,singularity
	nextflow run LabShengLi/nanome -profile [docker/singularity] \\
		--dsname DSNAME --input INPUT --genome GENOME

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
	  --queue		SLURM job submission queue name, e.g., 'gpu'
	  --qos			SLURM job submission QOS name, e.g., 'inference'
	  --gresOptions		SLURM job submission GPU allocation option, e.g., 'gpu:v100:1'
	  --time		SLURM job submission running time, e.g., '2h', '1d'
	  --memory		SLURM job submission memory, e.g., '32GB'

	  --projectCloud	Google Cloud Platform (GCP) project name for google-lifesciences
	  --config		Lifebit CloudOS config file, e.g., 'conf/executors/lifebit.config'

	-profile options:
	  Use this parameter to choose a predefined configuration profile. Profiles can give configuration presets for different compute environments.

	  test		A bundle of input params for ecoli test
	  test_human	A bundle of input params for human test
	  docker 	A generic configuration profile to be used with Docker, pulls software from Docker Hub: liuyangzzu/nanome:latest
	  singularity	A generic configuration profile to be used with Singularity, pulls software from: docker://liuyangzzu/nanome:latest
	  conda		Please only use conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity. Check our GitHub for how to install local conda environment
	  hpc		A generic configuration profile to be used on HPC cluster with SLURM
	  google	A generic configuration profile to be used on Google Cloud platform with 'google-lifesciences'

	Contact to https://github.com/LabShengLi/nanome/issues for bug report.
	""".stripIndent()
}

