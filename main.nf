#!/usr/bin/env nextflow


log.info """\
	NANOME - NF PIPELINE (v1.0)
	by The Jackon Laboratory
	Sheng Li Lab http://nanome.jax.org
	=================================
	dsname              :${params.dsname}
	input               :${params.input}
	benchmarking        :${params.benchmarking}
	eval                :${params.eval}
	debug               :${params.debug}
	online              :${params.online}
	=================================
	"""
	.stripIndent()



// Input for preprocessing (untar, seperate)
ch_input = params.benchmarking ? Channel.empty() : Channel.fromPath(params.input, checkIfExists: true)

// Benchmarking inputs
benchmarking_in_ch = params.benchmarking ? Channel
        .fromPath(params.input, type: 'file', checkIfExists: true) : Channel.empty()

fast5_tar_ch = Channel.fromPath(params.input)


// Inputs for methcalling pipelines (such as reference genome, deepsignal model, DeepMod Clustering data, Megalodon model, etc.)
hg38_tar_ch = Channel.fromPath(params.hg38_tar)
deepmod_ctar_ch = Channel.fromPath(params.deepmod_ctar)
deepsignel_model_tar_ch = Channel.fromPath(params.deepsignel_model_tar)
megalodon_model_tar_ch = Channel.fromPath(params.megalodon_model_tar)


// Get input from internet, and output to channel for later usage
process GetInputData{
	tag 'GetInputData'
	cache  'lenient'

	when:
	params.online

	input:
	file x1 from hg38_tar_ch
	file x2 from deepmod_ctar_ch
	file x3 from deepsignel_model_tar_ch
	file x4 from megalodon_model_tar_ch

	output:
	file "hg38" into hg38_fa_ch
	file "deepmod_c" into deepmod_c_ch
	file "deepsignal_model" into deepsignal_model_ch
	file "megalodon_model" into megalodon_model_ch

	"""
	set -x
	mkdir -p hg38
	tar xzf ${x1} -C hg38

	mkdir -p deepmod_c
	tar xzf ${x2} -C deepmod_c

	mkdir -p deepsignal_model
	tar xzf ${x3} -C deepsignal_model

	mkdir -p megalodon_model
	tar xzf ${x4} -C megalodon_model
    """
}


// The reference genome will be used later for: Resquiggle, Nanopolish, Megalodon, etc.
hg38_fa_ch.into{hg38_fa_ch1; hg38_fa_ch2;hg38_fa_ch3;hg38_fa_ch4;hg38_fa_ch5}


// Get the NA19240/NA12878 input fast5.tar files, output to Basecall/Megalodon channel
process GetOnlineFast5{
	tag 'GetOnlineFast5'
	cache  'lenient'

	when:
	params.online

	input:
	file x from fast5_tar_ch

	output:
	file "M_*" into online_out_ch

	"""
	set -x

	mkdir -p untar_${x}

	##tar -xf ${x} -C untar_${x}

	infn=${x}

	if [ "\${infn##*.}" = "tar" ]; then
		tar -xf ${x} -C untar_${x}
	elif [ "\${infn##*.}" = "gz" ]; then
		tar -xzf ${x} -C untar_${x}
	fi


	newx=${x}
	newx=\${newx%".tar"}

	mkdir -p M_\${newx}
	find untar_${x} -name '*.fast5' -exec mv {} M_\${newx}/ \\;
    """
}


// Check all tools work well on the platform
process UntarBenchmarking {
	// echo true
	tag 'UntarBench'
	cache  'lenient'
	// teminate all later processes if this process is not passed
	errorStrategy 'terminate'

	when:
	params.benchmarking && !params.online

	input:
	file x from benchmarking_in_ch

	output:
	file "untar_dir/BenchmarkingData/MB*" into benchmarking_out_ch

	when:
	! params.online

	"""
	echo $x

	infn=$x
	mkdir -p untar_dir
    if [ "\${infn##*.}" = "tar" ]; then
		tar -xf ${x} -C untar_dir
	elif [ "\${infn##*.}" = "gz" ]; then
		tar -xzf ${x} -C untar_dir
	fi
    """
}



// Check all tools work well on the platform
process EnvCheck {

	tag 'EnvCheck'
	// teminate all later processes if this process is not passed
	errorStrategy 'terminate'

	when:
	! params.online

	"""
	set -x

	pwd
    ls -la
    which conda
	## source activate nanocompare
	## . activate nanocompare
    conda env list
	echo ${params.genomeMotifC}
    ##ls -la ${params.genomeMotifC}

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

    ${params.GuppyDir}/bin/guppy_basecaller -v
    """
}


// untar file, seperate into N folders named 'M1', ..., 'M10', etc.
process Preprocess {
//	echo true
	errorStrategy 'ignore'
	cache  'lenient'

    input:
    file fast5_tar from ch_input

	output:
    file 'sept_dir/M*' into preprocess_out_ch

    when:
    ! params.benchmarking && ! params.online

    """
    set -x
    echo "Input=${fast5_tar}"

    infn=${fast5_tar}

    mkdir -p untar_dir
    if [ "${params.isfile}" = true ] ; then
        #mkdir -p untar_dir
        #tar xzf \${infn} -C untar_dir
        if [ "\${infn##*.}" = "tar" ]; then
			tar -xf \${infn} -C untar_dir
		elif [ "\${infn##*.}" = "gz" ]; then
			tar -xzf \${infn} -C untar_dir
		fi
    else
        echo "Multiple fast5.tar need to be untared"
		#filelist=\$(find \${infn}/ -type f \\( -iname \\*.fast5.tar -o -iname \\*.fast5.tar.gz -o -iname \\*.fast5 \\) )
		filelist=\$( find \${infn}/ -type f -name '*.tar.gz' )
		for fast5tarfn in \${filelist}; do
			echo "fn=\${fast5tarfn}"
			if [ "\${fast5tarfn##*.}" = "tar" ]; then
				tar -xf \${fast5tarfn} -C untar_dir
	        elif [ "\${fast5tarfn##*.}" = "gz" ]; then
				tar -xzf \${fast5tarfn} -C untar_dir
	        elif [ "\${fast5tarfn##*.}" = "fast5" ]; then
				cp \${fast5tarfn} untar_dir
	        fi
		done
    fi
    echo "### Untar input done. ###"

    mkdir -p sept_dir
    python ${workflow.projectDir}/model_params/FilesSeparatorNew.py untar_dir ${params.ntask} sept_dir
    echo "### Seperation fast5 files done. ###"
    """
}



// We collect all folders of fast5  files, and send into Channels for pipelines
preprocess_out_ch
	.mix(benchmarking_out_ch, online_out_ch)
	.into { basecall_input_ch; megalodon_in_ch }



// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
	tag "${x}"
	cache  'lenient'

	input:
    file x from basecall_input_ch.flatten()
//    file hg38_tar from hg38_tar_ch

    output:
    file 'basecall_dir/M*' into basecall_out_ch
    file "basecall_dir/${x}/*-sequencing_summary.txt" into qc_ch

    when:
    ! params.debug

    """
    set -x
    mkdir -p basecall_dir/${x}
    guppy_basecaller --input_path $x \
        --save_path "basecall_dir/${x}" \
        --config ${params.GUPPY_BASECALL_MODEL} \
        --num_callers 3 \
        --fast5_out \
        --recursive \
        --verbose_logs \
        ${params.GuppyGPUOptions}

        #--gpu_runners_per_device ${params.processors} \
        #--device auto

    ## After basecall, rename summary file names
	mv basecall_dir/${x}/sequencing_summary.txt basecall_dir/${x}/${x}-sequencing_summary.txt

    """
}

process QCStep {
//	tag "${x}"
	cache  'lenient'

	publishDir "outputs/${params.dsname}-qc-report" //, mode: "copy"

	input:
    file flist from qc_ch.collect()

    output:
    file "summary/*-sequencing_summary.txt" into qc_out_ch

    """
    echo $flist


	mkdir -p summary
    cp -L -f *-sequencing_summary.txt summary/
    """
}

basecall_out_ch
	.into { resquiggle_in_ch; nanopolish_in_ch; deepmod_in_ch }


resquiggle_in2_ch = resquiggle_in_ch.flatten().combine(hg38_fa_ch1)


// resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
	tag "${ttt[0]}"
	cache  'lenient'

	//clusterOptions '-q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem=50G'

	input:
//    file x from resquiggle_in_ch.flatten()
	file ttt from resquiggle_in2_ch  // TODO: how to pass [basecallDir, refFile] int two vars


    output:
    file 'resquiggle_dir/M*' into resquiggle_out_ch

    """
    set -x

    filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}

	mkdir -p resquiggle_dir/\${x}/workspace
	### only copy workspace files, due to nanpolish modify base folder at x
	cp -rf \${x}/workspace/* resquiggle_dir/\${x}/workspace

    tombo resquiggle --dna --processes ${params.processors} \
        --corrected-group ${params.correctedGroup} \
		--basecall-group Basecall_1D_000 --overwrite \
		resquiggle_dir/\${x}/workspace \${ref}/hg38.fa
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

	input:
//    file x from deepsignal_in_ch.flatten()
// TODO: how to add more than two pair files [bdir, refGenome, DeepSignalModel]
	file ttt from deepsignal_in2_ch

    output:
    file '*.tsv' into deepsignal_out_ch

    """
    set -x

    filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}
	refGenome=\${ref}/hg38.fa
	deepsignalModelDir=\${filepair[2]}


	deepsignal call_mods --input_path \${x}/workspace \
	    --model_path \${deepsignalModelDir}/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt \
		--result_file "${params.dsname}-N${params.ntask}-DeepSignal.batch_\${x}.CpG.deepsignal.call_mods.tsv" \
		--reference_path \${refGenome} \
		--corrected_group ${params.correctedGroup} \
		--nproc ${params.processors} \
		--is_gpu ${params.isgpu}
    """
}


// Tombo runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Tombo {
	tag "${x}"
	cache  'lenient'

	input:
    file x from tombo_in_ch.flatten()

    output:
    file '*.per_read_stats.bed' into tombo_out_ch

    when:
    ! params.top3

    """
	tombo detect_modifications alternative_model \
		--fast5-basedirs ${x}/workspace \
		--dna --standard-log-likelihood-ratio \
		--statistics-file-basename \
		${params.analysisPrefix}.batch_${x} \
		--per-read-statistics-basename ${params.analysisPrefix}.batch_${x} \
		--alternate-bases CpG \
		--processes ${params.processors} \
		--corrected-group ${params.correctedGroup}

	python ${workflow.projectDir}/model_params/Tombo_extract_per_read_stats.py ${workflow.projectDir}/model_params/${params.chromSizesFile} \
				"${params.analysisPrefix}.batch_${x}.CpG.tombo.per_read_stats" \
				"${params.analysisPrefix}.batch_${x}.CpG.tombo.per_read_stats.bed"
    """
}


megalodon_in2_ch = megalodon_in_ch.flatten().combine(hg38_fa_ch2)

// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${ttt[0]}"
	cache  'lenient'

	input:
//    file x from megalodon_in_ch.flatten()
//    file ref from hg38_fa_ch
	file ttt from megalodon_in2_ch  // TODO: how to pass [basecallDir, refFile] int two vars
	val x from 'abc'
	val ref from 'ref'

    output:
    file 'megalodon_results/*.per_read_modified_base_calls.txt' into megalodon_out_ch

    when:
    ! params.debug

    """
    set -x

	echo $ttt

	## TODO: reference to array index
	echo ${ttt[0]},${ttt[1]}


	filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}
	refGenome=\${ref}/hg38.fa


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
	    --guppy-params "--num_callers 5 --ipc_threads 80" \
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

	mv megalodon_results/per_read_modified_base_calls.txt megalodon_results/${params.analysisPrefix}.batch_\${x}.per_read_modified_base_calls.txt
    """
}


deepmod_in2_ch = deepmod_in_ch.flatten().combine(hg38_fa_ch4)


// DeepMod runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepMod {
	tag "${ttt[0]}"

	cache  'lenient'

	errorStrategy 'ignore'

	input:
	file ttt from deepmod_in2_ch
//    file x from deepmod_in_ch.flatten()

    output:
    file 'mod_output/batch_*_num' into deepmod_out_ch

    when:
    ! params.top3 && params.runDeepMod

    """
    set -x

    echo $ttt
	filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}
	refGenome=\${ref}/hg38.fa

    DeepMod.py detect \
			--wrkBase \${x}/workspace --Ref \${refGenome} \
			--Base C --modfile ${workflow.projectDir}/model_params/${params.deepModModel} \
			--FileID batch_\${x}_num \
			--threads ${params.processors} --move
    """
}


nanopolish_in2_ch =  nanopolish_in_ch.flatten().combine(hg38_fa_ch5)

// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Nanopolish {
	tag "${ttt[0]}"
	cache  'lenient'

	input:
//    file x from nanopolish_in_ch.flatten()
	file ttt from nanopolish_in2_ch

    output:
    file '*.tsv' into nanopolish_out_ch

    """
    set -x

    filepair=(${ttt})
	x=\${filepair[0]}
	ref=\${filepair[1]}
	refGenome=\${ref}/hg38.fa


    ### We put all fq and bam files into working dir, DO NOT affect the basecall dir
	##fastqFile=\${x}/reads.fq
	fastqFile=\${x}.reads.fq
	fastqNoDupFile="\${fastqFile}.noDups.fq"
	bamFileName="${params.analysisPrefix}.batch_\${x}.sorted.bam"


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

	python ${workflow.projectDir}/model_params/nanopore_nanopolish_preindexing_checkDups.py \${fastqFile} \${fastqNoDupFile}

	# Index the raw read with fastq
	nanopolish index -d \${x}/workspace \${fastqNoDupFile}

	minimap2 -t ${params.processors} -a -x map-ont \${refGenome} \${fastqNoDupFile} | samtools sort -T tmp -o \${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ threads \${bamFileName}
	echo "### samtools finished"

	echo "### Alignment step DONE"

	nanopolish call-methylation -t ${params.processors} -r \${fastqNoDupFile} -b \${bamFileName} -g \${refGenome} > ${params.analysisPrefix}.batch_\${x}.nanopolish.methylation_calls.tsv
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
	publishDir "outputs/${params.dsname}-methylation-callings" //, mode: "copy"

	input:
    file x from deepsignal_combine_in_ch

    output:
    file '*.combine.tsv' into deepsignal_combine_out_ch

    when:
    x.size() >= 1

    """
	echo ${x}
	touch ${params.dsname}.DeepSignal.combine.tsv
	cat ${x} > ${params.dsname}.DeepSignal.combine.tsv
    """
}


// Combine Tombo runs' all results together
process TomboCombine {
	publishDir "outputs/${params.dsname}-methylation-callings" //, mode: "copy"

	input:
    file x from tombo_combine_in_ch

    output:
    file '*.combine.tsv' into tombo_combine_out_ch

    when:
    x.size() >= 1

    """
	touch ${params.dsname}.Tombo.combine.tsv
	cat ${x} > ${params.dsname}.Tombo.combine.tsv
    """
}


// Combine Megalodon runs' all results together
process MgldnCombine {

	publishDir "outputs/${params.dsname}-methylation-callings" //, mode: "copy"

	input:
    file x from megalodon_combine_in_ch

    output:
    file '*.combine.tsv' into megalodon_combine_out_ch

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
    """
}


// Combine Nanopolish runs' all results together
process NplshCombine {
	publishDir "outputs/${params.dsname}-methylation-callings" //, mode: "copy"

	input:
    file x from nanopolish_combine_in_ch

    output:
    file '*.combine.tsv' into nanopolish_combine_out_ch

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
    """
}



// Combine DeepMod runs' all results together
process DpmodCombine {
	publishDir "outputs/${params.dsname}-methylation-callings" //, mode: "copy"

	input:
    file x from deepmod_combine_in_ch

    file deepmodCDir from deepmod_c_ch

    output:
    file '*.combine.tsv' into deepmod_combine_out_ch

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

    python ${workflow.projectDir}/model_params/sum_chr_mod.py \
        indir/ C ${params.dsname}.deepmod

	python ${workflow.projectDir}/model_params/hm_cluster_predict.py \
		indir/${params.dsname}.deepmod \
		\${DeepMod_Cluster_CDir} \
		${workflow.projectDir}/model_params/${params.clusterDeepModModel}  || true

	> ${params.dsname}.DeepModC.combine.tsv

	for f in \$(ls -1 indir/${params.dsname}.deepmod.chr*.C.bed)
	do
	  cat \$f >> ${params.dsname}.DeepModC.combine.tsv
	done

	> ${params.dsname}.DeepModC_clusterCpG.combine.tsv

	for f in \$(ls -1 indir/${params.dsname}.deepmod_clusterCpG.chr*.C.bed)
	do
	  cat \$f >> ${params.dsname}.DeepModC_clusterCpG.combine.tsv
	done
    """
}

//TODO: how sort the list???
deepsignal_combine_out_ch.concat(tombo_combine_out_ch,megalodon_combine_out_ch, \
	nanopolish_combine_out_ch,deepmod_combine_out_ch.flatten())
	.toSortedList()
	.into { readlevel_in_ch; sitelevel_in_ch }

//println("test_out_ch=" + test_out_ch)
//test_out_ch.view()


// Read level evaluations
process ReadLevelPerf {
	publishDir "outputs/${params.dsname}-nanome-analysis" //, mode: "copy"

	input: // TODO: I can not sort fileList by name, seems sorted by date????
    file fileList from readlevel_in_ch

    output:
    file "MethPerf-*" into readlevel_out_ch

    when:
    params.eval && (fileList.size() >= 1)

    """
    set -x

    echo ${fileList}

	# Sort file by my self
    flist=(\$(ls *.combine.tsv))

    echo \${flist[@]}

    for fn in \${flist[@]}
    do
        echo File: \$fn
        # head -n 3 \$fn
    done

    deepsignalFile=\$(find . -maxdepth 1 -name '*.DeepSignal.combine.tsv')
    tomboFile=\$(find . -maxdepth 1 -name '*.Tombo.combine.tsv')
    nanopolishFile=\$(find . -maxdepth 1 -name '*.Nanopolish.combine.tsv')
    deepmodFile=\$(find . -maxdepth 1 -name '*.DeepModC.combine.tsv')
    megalodonFile=\$(find . -maxdepth 1 -name '*.Megalodon.combine.tsv')

    export PYTHONPATH=${workflow.projectDir}:\${PYTHONPATH}

    ## python ${workflow.projectDir}/src/nanocompare/read_level_eval.py --help

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
		--report-joined -mpi -o . ###--enable-cache --using-cache

	echo "### Read level analysis DONE"
    """
}

// Site level correlation analysis
process SiteLevelCorr {
	publishDir "outputs/${params.dsname}-nanome-analysis" //, mode: "copy"

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
    flist=(\$(ls *.combine.tsv))
    echo \${flist[@]}

    deepsignalFile=\$(find . -maxdepth 1 -name '*.DeepSignal.combine.tsv')
    tomboFile=\$(find . -maxdepth 1 -name '*.Tombo.combine.tsv')
    nanopolishFile=\$(find . -maxdepth 1 -name '*.Nanopolish.combine.tsv')
    deepmodFile=\$(find . -maxdepth 1 -name '*.DeepModC_clusterCpG.combine.tsv')
    megalodonFile=\$(find . -maxdepth 1 -name '*.Megalodon.combine.tsv')

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
