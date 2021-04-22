#!/usr/bin/env nextflow


log.info """\
	NANOME - NF PIPELINE
	by The Jackon Laboratory
	Sheng Li Lab http://nanome.jax.org
	=================================
	dsname      :${params.dsname}
	input       :${params.input}
	benchmarking:${params.benchmarking}
	eval        :${params.eval}
	debug       :${params.debug}
	=================================
	"""
	.stripIndent()



// Input for preprocessing (untar, seperate)
ch_input = params.benchmarking ? Channel.empty() : Channel.fromPath(params.input, checkIfExists: true)

// Benchmarking inputs
benchmarking_in_ch = params.benchmarking ? Channel
        .fromPath(params.input, type: 'file', checkIfExists: true) : Channel.empty()

// Check all tools work well on the platform
process UntarBenchmarking {
	// echo true
	tag 'UntarBench'
	cache  'lenient'
	// teminate all later processes if this process is not passed
	errorStrategy 'terminate'

	when:
	params.benchmarking

	input:
	file x from benchmarking_in_ch

	output:
	file "untar_dir/BenchmarkingData/MB*" into benchmarking_out_ch

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
    ! params.benchmarking

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

fast5_tar_ch = Channel.fromPath(params.input_fast5_tar)

process UntarOnlineData{
	tag 'UntarOnlineData'
	cache  'lenient'

	when:
	params.online

	input:
	file x from fast5_tar_ch

	output:
	file "M${x}/*" into online_out_ch

	"""
	time tar -xf ${x} -C M${x}
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

    output:
    file 'basecall_dir/M*' into basecall_out_ch
    file "basecall_dir/${x}/*-sequencing_summary.txt" into qc_ch

    when:
    ! params.debug

    """
    mkdir -p basecall_dir/${x}
    ${params.GuppyDir}/bin/guppy_basecaller --input_path $x \
        --save_path "basecall_dir/${x}" \
        --config ${params.GUPPY_BASECALL_MODEL} \
        --num_callers 3 \
        --fast5_out \
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



// resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
	tag "${x}"
	cache  'lenient'

	//clusterOptions '-q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem=50G'

	input:
    file x from resquiggle_in_ch.flatten()

    output:
    file 'resquiggle_dir/M*' into resquiggle_out_ch

    """
    set -x
	mkdir -p resquiggle_dir/${x}/workspace
	### only copy workspace files, due to nanpolish modify base folder at x
	cp -rf ${x}/workspace/* resquiggle_dir/${x}/workspace

    tombo resquiggle --dna --processes ${params.processors} \
        --corrected-group ${params.correctedGroup} \
		--basecall-group Basecall_1D_000 --overwrite \
		resquiggle_dir/${x}/workspace ${params.refGenome}
    """
}


// duplicate resquiggle results for DeepSignal and Tombo
resquiggle_out_ch.into { deepsignal_in_ch; tombo_in_ch }


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${x}"

	cache  'lenient'

	input:
    file x from deepsignal_in_ch.flatten()

    output:
    file '*.tsv' into deepsignal_out_ch

    """
	deepsignal call_mods --input_path ${x}/workspace \
	    --model_path ${params.deepsignalModel}/bn_17.sn_360.epoch_9.ckpt \
		--result_file "${params.dsname}-N${params.ntask}-DeepSignal.batch_${x}.CpG.deepsignal.call_mods.tsv" \
		--reference_path ${params.refGenome} \
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


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${x}"
	cache  'lenient'

	input:
    file x from megalodon_in_ch.flatten()

    output:
    file 'megalodon_results/*.per_read_modified_base_calls.txt' into megalodon_out_ch

    when:
    ! params.debug

    """
    mkdir -p indir/${x}
    cp -rf ${x}/* indir/${x}/

    megalodon \
	    indir/${x} \
	    --overwrite \
	    --outputs basecalls mod_basecalls mappings \
	    per_read_mods mods mod_mappings \
	    per_read_refs \
	    --guppy-server-path ${params.GuppyDir}/bin/guppy_basecall_server \
	    --guppy-config ${params.MEGALODON_MODEL_FOR_GUPPY_CONFIG} \
	    --guppy-params "--num_callers 5 --ipc_threads 80" \
	    --reads-per-guppy-batch ${params.READS_PER_GUPPY_BATCH} \
	    --guppy-timeout ${params.GUPPY_TIMEOUT} \
	    --samtools-executable ${params.SAMTOOLS_PATH} \
	    --sort-mappings \
	    --mappings-format bam \
	    --reference ${params.refGenome} \
	    --mod-motif m CG 0 \
	    --mod-output-formats bedmethyl wiggle \
	    --write-mods-text \
	    --write-mod-log-probs \
	    --devices 0 \
	    --processes ${params.processors}

	mv megalodon_results/per_read_modified_base_calls.txt megalodon_results/${params.analysisPrefix}.batch_${x}.per_read_modified_base_calls.txt
    """
}


// DeepMod runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepMod {
	tag "${x}"

	cache  'lenient'

	errorStrategy 'ignore'

	input:
    file x from deepmod_in_ch.flatten()

    output:
    file 'mod_output/batch_*_num' into deepmod_out_ch

    when:
    ! params.top3 && params.runDeepMod

    """
    DeepMod.py detect \
			--wrkBase ${x}/workspace --Ref ${params.refGenome} \
			--Base C --modfile ${workflow.projectDir}/model_params/${params.deepModModel} \
			--FileID batch_${x}_num \
			--threads ${params.processors} --move

	# rm -rf mod_output/batch_*.done
    """
}


// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Nanopolish {
	tag "${x}"
	cache  'lenient'

	input:
    file x from nanopolish_in_ch.flatten()

    output:
    file '*.tsv' into nanopolish_out_ch

    """
    set -x

    ### We put all fq and bam files into working dir, DO NOT affect the basecall dir
	##fastqFile=${x}/reads.fq
	fastqFile=${x}.reads.fq
	fastqNoDupFile="\${fastqFile}.noDups.fq"
	bamFileName="${params.analysisPrefix}.batch_${x}.sorted.bam"


	echo \${fastqFile}
	echo \${fastqNoDupFile}

	##rm -rf \${fastqFile} \${fastqNoDupFile}
	##rm -rf ${x}/\${bamFileName} \${fastqNoDupFile}

	## Do alignment firstly
	touch \${fastqFile}
	for f in \$(ls -1 ${x}/*.fastq)
	do
		cat \$f >> \$fastqFile
		# echo "cat \$f >> \$fastqFile - COMPLETED"
	done

	python ${workflow.projectDir}/model_params/nanopore_nanopolish_preindexing_checkDups.py \${fastqFile} \${fastqNoDupFile}

	# Index the raw read with fastq
	nanopolish index -d ${x}/workspace \${fastqNoDupFile}

	minimap2 -t ${params.processors} -a -x map-ont ${params.refGenome} \${fastqNoDupFile} | samtools sort -T tmp -o \${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ threads \${bamFileName}
	echo "### samtools finished"

	echo "### Alignment step DONE"

	nanopolish call-methylation -t ${params.processors} -r \${fastqNoDupFile} -b \${bamFileName} -g ${params.refGenome} > ${params.analysisPrefix}.batch_${x}.nanopolish.methylation_calls.tsv
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

    output:
    file '*.combine.tsv' into deepmod_combine_out_ch

    when:
    x.size() >= 1

    """
    set -x

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
		${params.genomeMotifC} \
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
