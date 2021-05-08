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


// We collect all folders of fast5 files, and send into Channels for pipelines
if (params.input.endsWith(".filelist.txt")) { // filelist
	Channel.fromPath( params.input )
	    .splitCsv(header: false)
	    .map { file(it[0]) }
	    .toList()
	    .set{ fast5_tar_ch }

	// emit one time for each fast5.tar file
	fast5_tar_ch.flatten().into{fast5_tar_ch1; fast5_tar_ch2}
} else { // single file
	Channel.fromPath( params.input ).into {fast5_tar_ch1; fast5_tar_ch2}
}


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


// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
	tag "${fast5_tar.simpleName}"

	label 'with_gpus'

	input: // input here is not passed
    file fast5_tar from fast5_tar_ch1

    output:
    file "${fast5_tar.simpleName}_basecalled" into basecall_out_ch  // try to fix the christina proposed problems
    file "${fast5_tar.simpleName}_basecalled/${fast5_tar.simpleName}-sequencing_summary.txt" into qc_ch

    when:
    params.runBasecall

    """
    set -x

    mkdir -p untarDir

    infn=${fast5_tar}

	if [ "\${infn##*.}" == "tar" ]; then ### deal with tar
		tar -xf ${fast5_tar} -C untarDir
	elif [ "\${infn##*.}" == "gz" ]; then ### deal with tar.gz
		tar -xzf ${fast5_tar} -C untarDir
	else ### deal with ready folder
		mv  ${fast5_tar} untarDir
	fi

	# move fast5 files in tree folders into a single folder
	mkdir -p untarDir1
    find untarDir -name "*.fast5" -type f -exec mv {} untarDir1/ \\;

    mkdir -p ${fast5_tar.simpleName}_basecalled
    guppy_basecaller --input_path untarDir1 \
        --save_path "${fast5_tar.simpleName}_basecalled" \
        --config ${params.GUPPY_BASECALL_MODEL} \
        --num_callers ${params.GuppyNumCallers} \
        --fast5_out \
        --recursive \
        --verbose_logs ${params.GuppyGPUOptions}

    ## After basecall, rename and publishe summary file names
	mv ${fast5_tar.simpleName}_basecalled/sequencing_summary.txt ${fast5_tar.simpleName}_basecalled/${fast5_tar.simpleName}-sequencing_summary.txt

	## Clean
	## rm -rf untarDir untarDir1
    """
}


// Collect and output QC results for basecall
process QCExport {
	publishDir "${params.outputDir}" , mode: "copy"

	input:
    file flist from qc_ch.collect()

    output:
    file "${params.dsname}-qc-report/*-sequencing_summary.txt" into qc_out_ch

    """
	mkdir -p ${params.dsname}-qc-report
    cp -L -f *-sequencing_summary.txt ${params.dsname}-qc-report/
    """
}


// Duplicates basecall outputs
basecall_out_ch
	.into { resquiggle_in_ch; nanopolish_in_ch; deepmod_in_ch }


// resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
	tag "${basecallIndir.simpleName}"

	input:
	file basecallIndir from resquiggle_in_ch
	each file(reference_genome_tar) from Channel.fromPath(params.reference_genome_tar)

    output:
    file "${basecallIndir.simpleName}_resquiggle_dir" into resquiggle_out_ch

    when:
    params.runMethcall

    """
    set -x

    ## untar reference_genome
    tar -xzf ${reference_genome_tar}
	refGenome=${params.referenceGenome}

	### copy basecall workspace files, due to tombo resquiggle modify base folder
	mkdir -p ${basecallIndir.simpleName}_resquiggle_dir
	cp -rf ${basecallIndir}/* ${basecallIndir.simpleName}_resquiggle_dir/

    tombo resquiggle --dna --processes ${params.processors} \
        --corrected-group ${params.correctedGroup} \
		--basecall-group Basecall_1D_000 --overwrite \
		${basecallIndir.simpleName}_resquiggle_dir/workspace \${refGenome}

	echo "### Tombo methylation calling DONE"
    """
}


// duplicate resquiggle results for DeepSignal and Tombo
resquiggle_out_ch.into { deepsignal_in_ch; tombo_in_ch }


deepsignal_in2_ch =
	deepsignal_in_ch.flatten().combine(Channel.fromPath(params.reference_genome_tar))


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${ttt[0].simpleName}"

	label 'with_gpus'

	input:
	file ttt from deepsignal_in2_ch //[basecallDir, reference_genome]
	each file(deepsignal_model_tar) from Channel.fromPath(params.deepsignel_model_tar)

    output:
    file "DeepSignal.batch_${ttt[0].simpleName}.CpG.deepsignal.call_mods.tsv" into deepsignal_out_ch

    when:
    params.runMethcall

    """
    set -x

    tar -xzf ${ttt[1]}
    refGenome=${params.referenceGenome}

    tar -xzf ${deepsignal_model_tar}

	deepsignal call_mods --input_path ${ttt[0].simpleName}/workspace \
	    --model_path ./model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt \
		--result_file "DeepSignal.batch_${ttt[0].simpleName}.CpG.deepsignal.call_mods.tsv" \
		--reference_path \${refGenome} \
		--corrected_group ${params.correctedGroup} \
		--nproc ${params.processors} \
		--is_gpu ${params.DeepSignal_isgpu}
    """
}


// Tombo runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Tombo {
	tag "${resquiggleDir.simpleName}"

	input:
	each file(reference_genome_tar) from Channel.fromPath(params.reference_genome_tar)
	file resquiggleDir from tombo_in_ch

    output:
    file "batch_${resquiggleDir.simpleName}.CpG.tombo.per_read_stats.bed" into tombo_out_ch

    when:
    params.runMethcall

    """
    set -x
    tar -xzf ${reference_genome_tar}

	tombo detect_modifications alternative_model \
		--fast5-basedirs ${resquiggleDir}/workspace \
		--dna --standard-log-likelihood-ratio \
		--statistics-file-basename \
		batch_${resquiggleDir.simpleName} \
		--per-read-statistics-basename batch_${resquiggleDir.simpleName} \
		--alternate-bases CpG \
		--processes ${params.processors} \
		--corrected-group ${params.correctedGroup}

	python ${workflow.projectDir}/utils/tombo_extract_per_read_stats.py \
		${params.chromSizesFile} \
		"batch_${resquiggleDir.simpleName}.CpG.tombo.per_read_stats" \
		"batch_${resquiggleDir.simpleName}.CpG.tombo.per_read_stats.bed"
    """
}


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${fast5_tar.simpleName}"

	label 'with_gpus'

	input:
	file fast5_tar from fast5_tar_ch2
	each file (reference_genome_tar) from Channel.fromPath(params.reference_genome_tar)
	each file (megalodonModelTar) from Channel.fromPath(params.megalodon_model_tar)

    output:
    file "batch_${fast5_tar.simpleName}.per_read_modified_base_calls.txt" into megalodon_out_ch

    when:
    params.runMethcall && (params.runon == "gpu")

    """
    set -x

	## Get raw fast5 files
    mkdir -p untarDir
    infn=${fast5_tar}
	if [ "\${infn##*.}" == "tar" ]; then ### deal with tar
		tar -xf ${fast5_tar} -C untarDir
	elif [ "\${infn##*.}" == "gz" ]; then ### deal with tar.gz
		tar -xzf ${fast5_tar} -C untarDir
	else ### deal with ready folder
		mv  ${fast5_tar} untarDir
	fi

	## Get reference genome dir
	tar -xzf ${reference_genome_tar}
	refGenome=${params.referenceGenome}

	## Get megalodon model dir
	tar -xzf ${megalodonModelTar}

    megalodon \
	    untarDir \
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

	mv megalodon_results/per_read_modified_base_calls.txt batch_${fast5_tar.simpleName}.per_read_modified_base_calls.txt
    """
}


// DeepMod runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepMod {
	tag "${basecallDir.simpleName}"

	errorStrategy 'ignore'
	label 'with_gpus'

	input:
	each file (reference_genome_tar) from Channel.fromPath(params.reference_genome_tar)
	file basecallDir from deepmod_in_ch

    output:
    file "mod_output/batch_${basecallDir.simpleName}_num" into deepmod_out_ch

    when:
    params.runDeepMod

    """
    set -x

    wget ${params.DeepModGithub} --no-verbose
    tar -xzf v0.1.3.tar.gz
    DeepModProjectDir="DeepMod-0.1.3"

    tar -xzf ${reference_genome_tar}
	refGenome=${params.referenceGenome}

    DeepMod.py detect \
			--wrkBase ${basecallDir}/workspace --Ref \${refGenome} \
			--Base C --modfile \${DeepModProjectDir}/train_deepmod/${params.deepModModel} \
			--FileID batch_${basecallDir.simpleName}_num \
			--threads ${params.processors} ${params.DeepModMoveOptions}  #--move
    """
}


// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Nanopolish {
	tag "${basecallDir.simpleName}"

	input:
	file basecallDir from nanopolish_in_ch
	each file (reference_genome_tar) from Channel.fromPath(params.reference_genome_tar)

    output:
    file "batch_${basecallDir.simpleName}.nanopolish.methylation_calls.tsv" into nanopolish_out_ch

    when:
    params.runMethcall

    """
    set -x

	tar -xzf ${reference_genome_tar}
	refGenome=${params.referenceGenome}

	echo ${basecallDir}

    ### We put all fq and bam files into working dir, DO NOT affect the basecall dir
	fastqFile=reads.fq
	fastqNoDupFile="\${fastqFile}.noDups.fq"
	bamFileName="${params.dsname}.batch_${basecallDir.simpleName}.sorted.bam"

	echo \${fastqFile}
	echo \${fastqNoDupFile}

	## Do alignment firstly
	touch \${fastqFile}
	for f in \$(ls -1 ${basecallDir}/*.fastq)
	do
		cat \$f >> \$fastqFile
		# echo "cat \$f >> \$fastqFile - COMPLETED"
	done

	python ${workflow.projectDir}/utils/nanopore_nanopolish_preindexing_checkDups.py \${fastqFile} \${fastqNoDupFile}

	# Index the raw read with fastq
	nanopolish index -d ${basecallDir}/workspace \${fastqNoDupFile}

	minimap2 -t ${params.processors} -a -x map-ont \${refGenome} \${fastqNoDupFile} | samtools sort -T tmp -o \${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ threads \${bamFileName}
	echo "### samtools finished"
	echo "### Alignment step DONE"

	nanopolish call-methylation -t ${params.processors} -r \
		\${fastqNoDupFile} -b \${bamFileName} -g \${refGenome} \
		 > batch_${basecallDir.simpleName}.nanopolish.methylation_calls.tsv

	echo "### Nanopolish methylation calling DONE"
    """
}


// prepare combining results
nanopolish_combine_in_ch = nanopolish_out_ch.toList()
deepmod_combine_in_ch=deepmod_out_ch.toList()
megalodon_combine_in_ch = megalodon_out_ch.toList()
tombo_combine_in_ch = tombo_out_ch.toList()
deepsignal_combine_in_ch = deepsignal_out_ch.toList()


// Combine DeepSignal runs' all results together
process DpSigComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
    file x from deepsignal_combine_in_ch

    output:
    file "${params.dsname}.DeepSignal.combine.tsv.gz" into deepsignal_combine_out_ch

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
process TomboComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
    file x from tombo_combine_in_ch // list of tombo bed files

    output:
    file "${params.dsname}.Tombo.combine.bed.gz" into tombo_combine_out_ch

    when:
    x.size() >= 1

    """
	touch ${params.dsname}.Tombo.combine.bed
	cat ${x} > ${params.dsname}.Tombo.combine.bed

	gzip ${params.dsname}.Tombo.combine.bed
    """
}


// Combine Megalodon runs' all results together
process MgldnComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
    file x from megalodon_combine_in_ch

    output:
    file "${params.dsname}.Megalodon.combine.bed.gz" into megalodon_combine_out_ch

    when:
    x.size() >= 1

    """
	> ${params.dsname}.Megalodon.combine.tsv

	for fn in $x
	do
		break
	done
	#sed -n '1p' \${fn} > ${params.dsname}.Megalodon.combine.bed

    for fn in $x
	do
		sed '1d' \${fn} >> ${params.dsname}.Megalodon.combine.bed
	done

	gzip ${params.dsname}.Megalodon.combine.bed
    """
}


// Combine Nanopolish runs' all results together
process NplshComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"

	input:
    file x from nanopolish_combine_in_ch

    output:
    file "${params.dsname}.Nanopolish.combine.tsv.gz" into nanopolish_combine_out_ch

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
process DpmodComb {
	publishDir "${params.outputDir}/${params.dsname}-methylation-callings" , mode: "copy"
	label 'with_gpus'  // cluster model need gpu

	input:
    file x from deepmod_combine_in_ch
    file deepmod_c_tar_file from Channel.fromPath(params.deepmod_ctar)

    output:
    file '*.combine.bed.gz' into deepmod_combine_out_ch

    when:
    x.size() >= 1

    """
    set -x

    wget ${params.DeepModGithub} --no-verbose
    tar -xzf v0.1.3.tar.gz
    DeepModProjectDir="DeepMod-0.1.3"

    tar -xzf ${deepmod_c_tar_file}

	mkdir -p indir
    for dx in $x
    do
        mkdir -p indir/\$dx
        cp -rf \$dx/* indir/\$dx/
    done

    python ${workflow.projectDir}/utils/sum_chr_mod.py \
        indir/ C ${params.dsname}.deepmod ${params.DeepModSumChrSet}

	## TODO: deepmod_project_dir is not correct now
	## Only apply to human genome
	python ${workflow.projectDir}/utils/hm_cluster_predict.py \
		indir/${params.dsname}.deepmod \
		./C \
		\${DeepModProjectDir}/train_deepmod/${params.clusterDeepModModel}  || true

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
    megalodonFile=\$(find . -maxdepth 1 -name '*.Megalodon.combine.bed.gz')

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
    megalodonFile=\$(find . -maxdepth 1 -name '*.Megalodon.combine.bed.gz')

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
