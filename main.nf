#!/usr/bin/env nextflow

log.info """\
	NANOME - NF PIPELINE (v1.0)
	by Li Lab at The Jackon Laboratory
	http://nanome.jax.org
	=================================
	dsname              :${params.dsname}
	input               :${params.input}
	reference_genome    :${params.referenceGenome}
	chromSizesFile      :${params.chromSizesFile}
	runBasecall         :${params.runBasecall}
	runMethcall         :${params.runMethcall}
	eval                :${params.eval}
	=================================
	"""
	.stripIndent()

//	benchmarking        :${params.benchmarking}
//	debug               :${params.debug}

if (params.input.endsWith("filelist.txt")) { // filelist
	println("Detect filelist:");
	Channel.fromPath( params.input )
	    .splitCsv(header: false)
	    .map { file(it[0]) }
	    .toList()
	    .set{ fast5_tar_ch }

	fast5_tar_ch.flatten().into{fast5_tar_ch1; fast5_tar_ch2}
	fast5_tar_ch2.view()
} else { // single file
	println("Detect one file:");
	Channel.fromPath( params.input ).set{fast5_tar_ch1}
}


// Inputs for methcalling pipelines (such as reference genome, deepsignal model, DeepMod Clustering data, Megalodon model, etc.)
reference_genome_tar_ch = Channel.fromPath(params.reference_genome_tar)  //reference_genome_tar
deepmod_ctar_ch = Channel.fromPath(params.deepmod_ctar)
deepsignel_model_tar_ch = Channel.fromPath(params.deepsignel_model_tar)
megalodon_model_tar_ch = Channel.fromPath(params.megalodon_model_tar)


// Check all tools work well on the platform
process EnvCheck {

	tag 'EnvCheck'
	// teminate all later processes if this process is not passed
	errorStrategy 'terminate'
	label 'with_gpus'

	"""
	set -x

	which guppy_basecaller
    guppy_basecaller -v

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
    """
}


// Get input (model, reference_genome, etc.) from cloud drive zenodo, then output to channel for each tool
process GetRefGenome{
	tag 'GetRefGenome'
	cache  'lenient'

	when:
	params.runMethcall

	input:
	file x2 from deepmod_ctar_ch
	file x3 from deepsignel_model_tar_ch
	file x4 from megalodon_model_tar_ch
	file x5 from reference_genome_tar_ch

	output:
	file "reference_genome" into reference_genome_ch //hg38_fa_ch
	file "deepmod_c" into deepmod_c_ch
	file "DeepMod" into deepmod_prj_ch //all model online
	file "deepsignal_model" into deepsignal_model_ch
	file "megalodon_model" into megalodon_model_ch

	"""
	set -x

	mkdir -p deepmod_c
	tar xzf ${x2} -C deepmod_c

	mkdir -p deepsignal_model
	tar xzf ${x3} -C deepsignal_model

	tar xzf ${x4} -C .

	tar xzf ${x5} -C .

	## TODO: The Docker container do not have git,
	## how about not use container for this process????
	git clone https://github.com/WGLab/DeepMod.git
    """
}


// The reference genome will be used later for: Resquiggle, Nanopolish, Megalodon, etc.
reference_genome_ch.into{hg38_fa_ch1; hg38_fa_ch2;hg38_fa_ch3;hg38_fa_ch4;hg38_fa_ch5;hg38_fa_ch6}


// Get input fast5.tar files from local or online, output to Basecall/Megalodon channel
process GetFast5Files{
	tag 'GetFast5Files'
	cache  'lenient'

	input:
	/* Note: this channel may have three conditions (containing fast5 files):
				1. tar.gz
				2. tar
				3. folder
	*/
	file x from fast5_tar_ch1 // flattened, emit 1 at a time

	//TODO: why M_test_dir is ok in cloud, but M_*_dir is not ok ('// problem')? We want send multiple file to output here????
	//I checked file or path still not working in cloud, if I using *, failed, if I using M_${x}_dir, still failed.
	output:
	//file "M_*_dir" into fast5_dir_out_ch
	file "M_test_dir" into fast5_dir_out_ch

	"""
	set -x

	infn=${x}

	if [ "\${infn##*.}" == "tar" ]; then ### deal with tar
		mkdir -p untarDir
		tar -xf ${x} -C untarDir
	elif [ "\${infn##*.}" == "gz" ]; then ### deal with tar.gz
		mkdir -p untarDir
		tar -xzf ${x} -C untarDir
	else ### deal with ready folder
		cp -d ${x} M_${x}_dir
		exit 0
	fi

	## convert file name, replace . with _, no tar.gz suffix now,
	## example: demo.fast5.tar.gz -> demo_fast5
	newx=${x}
	newx=\${newx%".tar.gz"}
	newx=\${newx%".tar"}
	newx="\${newx//./_}"

	echo \${newx}

	mv untarDir M_test_dir

    """
}


// We collect all folders of fast5 files, and send into Channels for pipelines
fast5_dir_out_ch.into { basecall_input_ch; megalodon_in_ch; fast5_dir_out_ch2 }

println(fast5_dir_out_ch2)

// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
	tag "${x}"
	cache  'lenient'
	label 'with_gpus'

	input: // input here is not passed
    file x from basecall_input_ch.flatten()  // we are going to solve this input failed passing problems

    output:
    file "basecall_dir/${x}_basecalled" into basecall_out_ch  // try to fix the christina proposed problems
    file "basecall_dir/${x}_basecalled/*-sequencing_summary.txt" into qc_ch

    when:
    params.runBasecall

    """
    set -x

    mkdir -p basecall_dir/${x}_basecalled

    guppy_basecaller --input_path ${x} \
        --save_path "basecall_dir/${x}_basecalled" \
        --config ${params.GUPPY_BASECALL_MODEL} \
        --num_callers ${params.GuppyNumCallers} \
        --fast5_out \
        --recursive \
        --verbose_logs ${params.GuppyGPUOptions}

    ## After basecall, rename summary file names
	mv basecall_dir/${x}_basecalled/sequencing_summary.txt basecall_dir/${x}_basecalled/${x}-sequencing_summary.txt

    """
}


// Collect and output QC results for basecall
process QCExport {
	cache  'lenient'

	publishDir "${params.outputDir}" , mode: "copy"

	input:
    file flist from qc_ch.collect()

    output:
    file "${params.dsname}-qc-report/*-sequencing_summary.txt" into qc_out_ch

    """
    echo $flist


	mkdir -p ${params.dsname}-qc-report
    cp -L -f *-sequencing_summary.txt ${params.dsname}-qc-report/
    """
}


// Duplicates basecall outputs
basecall_out_ch
	.into { resquiggle_in_ch; nanopolish_in_ch; deepmod_in_ch }


// Duplicates resquiggle outputs
resquiggle_in2_ch = resquiggle_in_ch.flatten().combine(hg38_fa_ch1)


// resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
	tag "${ttt[0]}"
	cache  'lenient'

	input:
	file ttt from resquiggle_in2_ch  // TODO: how to pass [basecallDir, refFile] int two vars

    output:
    file 'resquiggle_dir/M*' into resquiggle_out_ch

    when:
    params.runMethcall

    """
    set -x

    filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}

	refGenome=\${ref}/hg38.fasta

	refGenome=${params.referenceGenome}

	mkdir -p resquiggle_dir/\${x}/workspace
	### only copy workspace files, due to nanpolish modify base folder at x
	cp -rf \${x}/workspace/* resquiggle_dir/\${x}/workspace

    tombo resquiggle --dna --processes ${params.processors} \
        --corrected-group ${params.correctedGroup} \
		--basecall-group Basecall_1D_000 --overwrite \
		resquiggle_dir/\${x}/workspace \${refGenome}
    """
}


// duplicate resquiggle results for DeepSignal and Tombo
resquiggle_out_ch.into { deepsignal_in_ch; tombo_in_ch }


deepsignal_in2_ch =
	deepsignal_in_ch.flatten().combine(hg38_fa_ch3).combine(deepsignal_model_ch)



// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${ttt[0]}"

	cache  'lenient'
	label 'with_gpus'

	input:
//    file x from deepsignal_in_ch.flatten()
// TODO: how to add more than two pair files [bdir, refGenome, DeepSignalModel]
	file ttt from deepsignal_in2_ch

    output:
    file '*.tsv' into deepsignal_out_ch

    when:
    params.runMethcall

    """
    set -x

    filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}
	refGenome=\${ref}/hg38.fasta
	refGenome=${params.referenceGenome}

	deepsignalModelDir=\${filepair[2]}


	deepsignal call_mods --input_path \${x}/workspace \
	    --model_path \${deepsignalModelDir}/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt \
		--result_file "${params.dsname}-DeepSignal.batch_\${x}.CpG.deepsignal.call_mods.tsv" \
		--reference_path \${refGenome} \
		--corrected_group ${params.correctedGroup} \
		--nproc ${params.processors} \
		--is_gpu ${params.DeepMod_isgpu}
    """
}


tombo_in_ch.flatten().combine(hg38_fa_ch6).set{tombo_in_ch1}


// Tombo runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Tombo {
	tag "${ttt[0]}"
	cache  'lenient'
//	label 'with_gpus'

	input:
    file ttt from tombo_in_ch1 // [resquiggleDir, reference_genome]

    output:
    file '*.per_read_stats.bed' into tombo_out_ch

    when:
    params.runMethcall

    """
    x=${ttt[0]}
    reference_genome_dir=${ttt[1]}

	tombo detect_modifications alternative_model \
		--fast5-basedirs \${x}/workspace \
		--dna --standard-log-likelihood-ratio \
		--statistics-file-basename \
		${params.dsname}.batch_\${x} \
		--per-read-statistics-basename ${params.dsname}.batch_\${x} \
		--alternate-bases CpG \
		--processes ${params.processors} \
		--corrected-group ${params.correctedGroup}

	python ${workflow.projectDir}/utils/Tombo_extract_per_read_stats.py ${params.chromSizesFile} \
				"${params.dsname}.batch_\${x}.CpG.tombo.per_read_stats" \
				"${params.dsname}.batch_\${x}.CpG.tombo.per_read_stats.bed"
    """
}


// Megalodon will need 3 input together for each emit: [fast5Dir, hg38ReferenceGenome, ModelDir]
megalodon_in2_ch = megalodon_in_ch.flatten().combine(hg38_fa_ch2).combine(megalodon_model_ch)


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${ttt[0]}"
	cache  'lenient'

	//errorStrategy 'ignore'
	label 'with_gpus'

	input:
	file ttt from megalodon_in2_ch  // TODO: how to pass [basecallDir, refFile] int two vars

    output:
    file 'megalodon_results/*.per_read_modified_base_calls.txt' into megalodon_out_ch

    when:
    params.runMethcall

	// TODO: how to specify Megalodon / Guppy model using full path? or put them into folder
    """
    set -x

    ### git clone https://github.com/nanoporetech/rerio.git

	echo $ttt

	## TODO: reference to array index
	echo ${ttt[0]},${ttt[1]}, ${ttt[2]}


	filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}
	refGenome=\${ref}/hg38.fasta
	refGenome=${params.referenceGenome}


    mkdir -p indir/\${x}
    cp -rf \${x}/* indir/\${x}/

    megalodon \
	    indir/\${x} \
	    --overwrite \
	    --outputs basecalls mod_basecalls mappings \
	    per_read_mods mods mod_mappings \
	    per_read_refs \
	    --guppy-server-path guppy_basecall_server \
	    --guppy-config ${params.MEGALODON_MODEL_FOR_GUPPY_CONFIG} \
	    --guppy-params "-d ./megalodon_model/ --num_callers 5 --ipc_threads 80" \
	    --reads-per-guppy-batch ${params.READS_PER_GUPPY_BATCH} \
	    --guppy-timeout ${params.GUPPY_TIMEOUT} \
	    --samtools-executable ${params.SAMTOOLS_PATH} \
	    --sort-mappings \
	    --mappings-format bam \
	    --reference \${refGenome} \
	    --mod-motif m CG 0 \
	    --mod-output-formats bedmethyl wiggle \
	    --write-mods-text \
	    --write-mod-log-probs \
	    --devices 0 \
	    --processes ${params.processors}

	mv megalodon_results/per_read_modified_base_calls.txt megalodon_results/${params.dsname}.batch_\${x}.per_read_modified_base_calls.txt
    """
}


deepmod_prj_ch.into{deepmod_prj_ch1; deepmod_prj_ch2}

deepmod_in2_ch = deepmod_in_ch.flatten().combine(hg38_fa_ch4).combine(deepmod_prj_ch1)

// DeepMod runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepMod {
	tag "${ttt[0]}"

	cache  'lenient'

	errorStrategy 'ignore'
	label 'with_gpus'

	input:
	file ttt from deepmod_in2_ch  // [basecallDir, refGenome, DeepModPrj]

    output:
    file 'mod_output/batch_*_num' into deepmod_out_ch

    when:
    params.runDeepMod && params.runMethcall

    """
    set -x

    echo $ttt
	filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}
	refGenome=\${ref}/hg38.fasta
	refGenome=${params.referenceGenome}

	DeepModProjectDir=${ttt[2]}

    DeepMod.py detect \
			--wrkBase \${x}/workspace --Ref \${refGenome} \
			--Base C --modfile \${DeepModProjectDir}/train_deepmod/${params.deepModModel} \
			--FileID batch_\${x}_num \
			--threads ${params.processors} ${params.DeepModMoveOptions}  #--move
    """
}


nanopolish_in2_ch =  nanopolish_in_ch.flatten().combine(hg38_fa_ch5)

// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Nanopolish {
	tag "${ttt[0]}"
	cache  'lenient'

	input:
	file ttt from nanopolish_in2_ch

    output:
    file '*.tsv' into nanopolish_out_ch

    when:
    params.runMethcall

    """
    set -x

    filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}
	refGenome=\${ref}/hg38.fasta
	refGenome=${params.referenceGenome}


    ### We put all fq and bam files into working dir, DO NOT affect the basecall dir
	##fastqFile=\${x}/reads.fq
	fastqFile=\${x}.reads.fq
	fastqNoDupFile="\${fastqFile}.noDups.fq"
	bamFileName="${params.dsname}.batch_\${x}.sorted.bam"


	echo \${fastqFile}
	echo \${fastqNoDupFile}

	##rm -rf \${fastqFile} \${fastqNoDupFile}
	##rm -rf \${x}/\${bamFileName} \${fastqNoDupFile}

	## Do alignment firstly
	touch \${fastqFile}
	for f in \$(ls -1 \${x}/*.fastq)
	do
		cat \$f >> \$fastqFile
		# echo "cat \$f >> \$fastqFile - COMPLETED"
	done

	python ${workflow.projectDir}/utils/nanopore_nanopolish_preindexing_checkDups.py \${fastqFile} \${fastqNoDupFile}

	# Index the raw read with fastq
	nanopolish index -d \${x}/workspace \${fastqNoDupFile}

	minimap2 -t ${params.processors} -a -x map-ont \${refGenome} \${fastqNoDupFile} | samtools sort -T tmp -o \${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ threads \${bamFileName}
	echo "### samtools finished"

	echo "### Alignment step DONE"

	nanopolish call-methylation -t ${params.processors} -r \${fastqNoDupFile} -b \${bamFileName} -g \${refGenome} > ${params.dsname}.batch_\${x}.nanopolish.methylation_calls.tsv
    """
}


// prepare combining results
nanopolish_combine_in_ch = nanopolish_out_ch.toList()
deepmod_combine_in_ch=deepmod_out_ch.toList()
megalodon_combine_in_ch = megalodon_out_ch.toList()
tombo_combine_in_ch = tombo_out_ch.toList()
deepsignal_combine_in_ch = deepsignal_out_ch.toList()


// Combine DeepSignal runs' all results together
process DpSigCombine {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"
	//label 'with_gpus'

	input:
    file x from deepsignal_combine_in_ch

    output:
    file '*.combine.tsv.gz' into deepsignal_combine_out_ch

    when:
    x.size() >= 1

    """
	echo ${x}
	touch ${params.dsname}.DeepSignal.combine.tsv
	cat ${x} > ${params.dsname}.DeepSignal.combine.tsv

	gzip ${params.dsname}.DeepSignal.combine.tsv

    """
}


// Combine Tombo runs' all results together
process TomboCombine {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"
	//label 'with_gpus'

	input:
    file x from tombo_combine_in_ch // list of tombo bed files

    output:
    file '*.combine.bed.gz' into tombo_combine_out_ch

    when:
    x.size() >= 1

    """
	touch ${params.dsname}.Tombo.combine.bed
	cat ${x} > ${params.dsname}.Tombo.combine.bed

	gzip ${params.dsname}.Tombo.combine.bed
    """
}


// Combine Megalodon runs' all results together
process MgldnCombine {
	//label 'with_gpus'

	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
    file x from megalodon_combine_in_ch

    output:
    file '*.combine.tsv.gz' into megalodon_combine_out_ch

    when:
    x.size() >= 1

    """
	> ${params.dsname}.Megalodon.combine.tsv

	for fn in $x
	do
		break
	done
	#sed -n '1p' \${fn} > ${params.dsname}.Megalodon.combine.tsv

    for fn in $x
	do
		sed '1d' \${fn} >> ${params.dsname}.Megalodon.combine.tsv
	done

	gzip ${params.dsname}.Megalodon.combine.tsv
    """
}


// Combine Nanopolish runs' all results together
process NplshCombine {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"
	//label 'with_gpus'

	input:
    file x from nanopolish_combine_in_ch

    output:
    file '*.combine.tsv.gz' into nanopolish_combine_out_ch

    when:
    x.size() >= 1

    """
    set -x

    > ${params.dsname}.Nanopolish.combine.tsv

    for fn in $x
	do
		break
	done
	sed -n '1p' \${fn} > ${params.dsname}.Nanopolish.combine.tsv

    for fn in $x
	do
		sed '1d' \${fn} >> ${params.dsname}.Nanopolish.combine.tsv
	done

	gzip ${params.dsname}.Nanopolish.combine.tsv
    """
}



// Combine DeepMod runs' all results together
process DpmodCombine {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"
	label 'with_gpus'  // cluster model need gpu

	input:
    file x from deepmod_combine_in_ch
	file deepmod_project_dir from deepmod_prj_ch2
    file deepmodCDir from deepmod_c_ch

    output:
    file '*.combine.bed.gz' into deepmod_combine_out_ch

    when:
    x.size() >= 1

    """
    set -x

    echo ${deepmodCDir}/C

    DeepMod_Cluster_CDir=${deepmodCDir}/C

	mkdir -p indir
    for dx in $x
    do
        mkdir -p indir/\$dx
        cp -rf \$dx/* indir/\$dx
    done

    python ${workflow.projectDir}/utils/sum_chr_mod.py \
        indir/ C ${params.dsname}.deepmod ${params.DeepModSumChrSet}

	## Only apply to human genome
	python ${workflow.projectDir}/utils/hm_cluster_predict.py \
		indir/${params.dsname}.deepmod \
		\${DeepMod_Cluster_CDir} \
		${deepmod_project_dir}/train_deepmod/${params.clusterDeepModModel}  || true

	> ${params.dsname}.DeepModC.combine.bed

	## Note: for ecoli data, no chr*, but N*
	for f in \$(ls -1 indir/${params.dsname}.deepmod.*.C.bed)
	do
	  cat \$f >> ${params.dsname}.DeepModC.combine.bed
	done

	> ${params.dsname}.DeepModC_clusterCpG.combine.bed

	for f in \$(ls -1 indir/${params.dsname}.deepmod_clusterCpG.*.C.bed)
	do
	  cat \$f >> ${params.dsname}.DeepModC_clusterCpG.combine.bed
	done

	gzip ${params.dsname}.DeepModC.combine.bed
	gzip ${params.dsname}.DeepModC_clusterCpG.combine.bed

	echo "###DeepMod combine DONE###"
    """
}


//TODO: how sort the list???
deepsignal_combine_out_ch.concat(tombo_combine_out_ch,megalodon_combine_out_ch, \
	nanopolish_combine_out_ch,deepmod_combine_out_ch.flatten())
	.toSortedList()
	.into { readlevel_in_ch; sitelevel_in_ch }

// Read level evaluations
process ReadLevelPerf {
	publishDir "${params.outputDir}/${params.dsname}-nanome-analysis" , mode: "copy"

	input: // TODO: I can not sort fileList by name, seems sorted by date????
    file fileList from readlevel_in_ch

    output:
    file "MethPerf-*" into readlevel_out_ch

    when:
    params.eval && (fileList.size() >= 1)

    """
    set -x

	# Sort files
    flist=(\$(ls *.combine.{tsv,bed}.gz))

    echo \${flist[@]}

    deepsignalFile=\$(find . -maxdepth 1 -name '*.DeepSignal.combine.tsv.gz')
    tomboFile=\$(find . -maxdepth 1 -name '*.Tombo.combine.bed.gz')
    nanopolishFile=\$(find . -maxdepth 1 -name '*.Nanopolish.combine.tsv.gz')
    deepmodFile=\$(find . -maxdepth 1 -name '*.DeepModC.combine.bed.gz')
    megalodonFile=\$(find . -maxdepth 1 -name '*.Megalodon.combine.tsv.gz')

    export PYTHONPATH=${workflow.projectDir}:\${PYTHONPATH}

	## Read level evaluations
	python ${workflow.projectDir}/src/nanocompare/read_level_eval.py \
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

    output:
    file "MethCorr-*" into sitelevel_out_ch

    when:
    params.eval && (fileList.size() >= 1)

    """
    set -x

	# Sort file by my self
    flist=(\$(ls *.combine.{tsv,bed}.gz))
    echo \${flist[@]}

    deepsignalFile=\$(find . -maxdepth 1 -name '*.DeepSignal.combine.tsv.gz')
    tomboFile=\$(find . -maxdepth 1 -name '*.Tombo.combine.bed.gz')
    nanopolishFile=\$(find . -maxdepth 1 -name '*.Nanopolish.combine.tsv.gz')
    deepmodFile=\$(find . -maxdepth 1 -name '*.DeepModC_clusterCpG.combine.bed.gz')
    megalodonFile=\$(find . -maxdepth 1 -name '*.Megalodon.combine.tsv.gz')

    export PYTHONPATH=${workflow.projectDir}:\${PYTHONPATH}

	## Site level evaluations
	python ${workflow.projectDir}/src/nanocompare/site_level_eval.py \
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
