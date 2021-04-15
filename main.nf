#!/usr/bin/env nextflow

Channel
    .fromPath(params.indata)
    .ifEmpty { exit 1, "Cannot find input file"}
    .set {ch_input}

aa = Channel
    .from( "Megalodon", "Tombo", "DeepMod", "DeepSignal", "Nanopolish" )
    .toList()

print(aa)


process FirstCheck {
	input:
	val aa1 from aa

	tag 'checkEnv'

	"""
	set -x

	echo $aa1
	echo ${aa1[0]}
	echo ${aa1[1]}
	echo ${aa1[2]}
	echo ${aa1[3]}
	echo ${aa1[4]}

    echo hello
	pwd
    ls -la
    conda env list

	echo ${params.genomeMotifC}
    ls -la ${params.genomeMotifC}

    tombo -v
    nanopolish --version
    megalodon -v
    deepsignal
    DeepMod.py

    ${params.GuppyDir}/bin/guppy_basecaller -v
    """
}


// untar file, seperate into N folders named 'M1', ..., 'M10', etc.
process Preprocess {
    input:
    file fast5_tar from ch_input

	output:
    file 'sept_dir/M*' into fast5Inputs

    when:
    true

    """
    set -x

    echo ${workflow.projectDir}/${params.chromSizesFile}


    echo ${fast5_tar}
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
    python ${workflow.projectDir}/src/nanocompare/methcall/FilesSeparatorNew.py untar_dir ${params.ntask} sept_dir
    echo "### Seperation fast5 files done. ###"
    """
}


fast5Inputs.into { fast5Inputs1ForBasecall; fast5Inputs2ForMegalodon }


// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
	tag "${x}"

	input:
    file x from fast5Inputs1ForBasecall.flatten()

    output:
    file 'basecall_dir/M*' into basecallOutputs

    """
    mkdir -p basecall_dir/${x}
    ${params.GuppyDir}/bin/guppy_basecaller --input_path $x \
        --save_path "basecall_dir/${x}" \
        --config ${params.GUPPY_BASECALL_MODEL} \
        --gpu_runners_per_device ${params.processors} \
        --num_callers 3 \
        --fast5_out \
        --verbose_logs \
        --device auto
    """
}


basecallOutputs.into { basecallOutputs1ForResquiggle; basecallOutputs2ForNanopolish; basecallOutputs3ForDeepMod }


// resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process Resquiggle {
	tag "${x}"

	input:
    file x from basecallOutputs1ForResquiggle.flatten()

    output:
    file 'resquiggle_dir/M*' into resquiggleOutputs

    """
    set -x
	mkdir -p resquiggle_dir/${x}
	cp -rf ${x}/* resquiggle_dir/${x}

    tombo resquiggle --dna --processes ${params.processors} \
        --corrected-group ${params.correctedGroup} \
		--basecall-group Basecall_1D_000 --overwrite \
		resquiggle_dir/${x}/workspace ${params.refGenome}
    """
}


// duplicate resquiggle results for DeepSignal and Tombo
resquiggleOutputs.into { resquiggleOutputs1ForDeepSignal; resquiggleOutputs2ForTombo }


// DeepSignal runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process DeepSignal {
	tag "${x}"

	input:
    file x from resquiggleOutputs1ForDeepSignal.flatten()

    output:
    file '*.tsv' into deepsignalOutput

    """
	deepsignal call_mods --input_path ${x}/workspace \
	    --model_path ${params.deepsignalModel}/bn_17.sn_360.epoch_9.ckpt \
		--result_file "${params.dsname}-N${params.ntask}-DeepSignal.batch_${x}.CpG.deepsignal.call_mods.tsv" \
		--reference_path ${params.refGenome} \
		--corrected_group ${params.correctedGroup} \
		--nproc ${params.processors} \
		--is_gpu ${params.isGPU}
    """
}


// Tombo runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Tombo {
	tag "${x}"

	input:
    file x from resquiggleOutputs2ForTombo.flatten()

    output:
    file '*.per_read_stats.bed' into tomboOutput

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

	python ${workflow.projectDir}/src/nanocompare/methcall/Tombo_extract_per_read_stats.py ${workflow.projectDir}/model_params/${params.chromSizesFile} \
				"${params.analysisPrefix}.batch_${x}.CpG.tombo.per_read_stats" \
				"${params.analysisPrefix}.batch_${x}.CpG.tombo.per_read_stats.bed"
    """
}


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	tag "${x}"

	input:
    file x from fast5Inputs2ForMegalodon.flatten()

    output:
    file 'megalodon_results/*.per_read_modified_base_calls.txt' into megalodonOutput

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
	    --guppy-config ${params.GUPPY_MOD_CONFIG} \
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

	input:
    file x from basecallOutputs3ForDeepMod.flatten()

    output:
    file 'mod_output/batch_*_num' into deepmodOutput

    when:
    params.RunDeepMod == "true"

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

	input:
    file x from basecallOutputs2ForNanopolish.flatten()

    output:
    file '*.tsv' into nanopolishOutput

    """
    set -x
	fastqFile=${x}/reads.fq
	fastqNoDupFile="\${fastqFile}.noDups.fq"
	bamFileName="${params.analysisPrefix}.batch_${x}.sorted.bam"


	echo \${fastqFile}
	echo \${fastqNoDupFile}

	rm -rf \${fastqFile} \${fastqNoDupFile}
	rm -rf ${x}/\${bamFileName} \${fastqNoDupFile}

	## Do alignment firstly
	touch \${fastqFile}
	for f in \$(ls -1 ${x}/*.fastq)
	do
		cat \$f >> \$fastqFile
		# echo "cat \$f >> \$fastqFile - COMPLETED"
	done

	python ${workflow.projectDir}/src/nanocompare/methcall/nanopore_nanopolish.NA19240_pipeline.step_02.preindexing_checkDups.py \${fastqFile} \${fastqNoDupFile}

	nanopolish index -d ${x}/workspace \${fastqNoDupFile}

	minimap2 -t ${params.processors} -a -x map-ont ${params.refGenome} \${fastqNoDupFile} | samtools sort -T tmp -o ${x}/\${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ threads ${x}/\${bamFileName}
	echo "### samtools finished"

	echo "### Alignment step DONE"

	nanopolish call-methylation -t ${params.processors} -r \${fastqNoDupFile} -b ${x}/\${bamFileName} -g ${params.refGenome} > ${params.analysisPrefix}.batch_${x}.nanopolish.methylation_calls.tsv
    """
}


// prepare combining results
nanopolishResults = nanopolishOutput.toList()
deepmodResults=deepmodOutput.toList()
megalodonResults = megalodonOutput.toList()
tomboResults = tomboOutput.toList()
deepsignalResults = deepsignalOutput.toList()


// Combine DeepSignal runs' all results together
process DpSigCombine {
	input:
    file x from deepsignalResults

    output:
    file '*.combine.tsv' into deepsignalCombineResult

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
	input:
    file x from tomboResults

    output:
    file '*.combine.tsv' into tomboCombineResult

    when:
    x.size() >= 1

    """
	touch ${params.dsname}.Tombo.combine.tsv
	cat ${x} > ${params.dsname}.Tombo.combine.tsv
    """
}


// Combine Megalodon runs' all results together
process MgldnCombine {
	input:
    file x from megalodonResults

    output:
    file '*.combine.tsv' into megalodonCombineResult

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
	input:
    file x from nanopolishResults

    output:
    file '*.combine.tsv' into nanopolishCombineResult

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
	input:
    file x from deepmodResults

    output:
    file '*.combine.tsv' into deepmodCombineResult

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

deepsignalCombineResult.concat(tomboCombineResult,megalodonCombineResult, \
	nanopolishCombineResult,deepmodCombineResult.flatten())
	.toSortedList()
	.set{allCombinedResultsList}

allCombinedResultsList.into { allCombinedResultsList_ch1; allCombinedResultsList_ch2 }


// Evaluation on combined results
process ReadLevelPerf {
	input: // TODO: I can not sort fileList by name, seems sorted by date????
    file fileList from allCombinedResultsList_ch1

    output:
    file "MethPerf-*" into ReadLevelPerfOut

    when:
    params.eval == "true"

    """
    set -x

	# Sort file by my self
    flist=(\$(ls *.combine.tsv))

    echo \${flist[@]}

    for fn in \${flist[@]}
    do
        echo File: \$fn
        head -n 3 \$fn
    done

    export PYTHONPATH=${workflow.projectDir}:\${PYTHONPATH}

    ## python ${workflow.projectDir}/src/nanocompare/read_level_eval.py --help

	## Read level evaluations
	python ${workflow.projectDir}/src/nanocompare/read_level_eval.py \
		--calls DeepSignal:\${flist[2]} \
				Tombo:\${flist[5]} \
				Nanopolish:\${flist[4]} \
				DeepMod.C:\${flist[1]} \
				Megalodon:\${flist[3]} \
		--bgtruth 'bismark:/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDA.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDF.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz' \
		--runid MethPerf-${params.runid} \
		--dsname ${params.dsname} \
		--min-bgtruth-cov ${params.bgtruthCov} --report-joined --mpi --enable-cache --using-cache -o .

	echo "### Read level analysis DONE"
    """
}

// Site level correlation analysis
process SiteLevelCorr {
	input:
    file perfDir from ReadLevelPerfOut
    file fileList from allCombinedResultsList_ch2

    output:
    file "MethCorr-*" into SiteLevelCorrOut

    when:
    params.eval == "true"

    """
    set -x

	# Sort file by my self
    flist=(\$(ls *.combine.tsv))

    echo \${flist[@]}

    for fn in \${flist[@]}
    do
        echo File: \$fn
        head -n 3 \$fn
    done

    export PYTHONPATH=${workflow.projectDir}:\${PYTHONPATH}

	## Site level evaluations
	python ${workflow.projectDir}/src/nanocompare/site_level_eval.py \
		--calls DeepSignal:\${flist[2]} \
				Tombo:\${flist[5]} \
				Nanopolish:\${flist[4]} \
				DeepMod.Cluster:\${flist[0]} \
				Megalodon:\${flist[3]} \
		--bgtruth 'bismark:/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDA.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/HL60/HL60_RRBS_ENCFF000MDF.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz' \
		--runid MethPerf-${params.runid} \
		--dsname ${params.dsname} \
		--min-bgtruth-cov ${params.bgtruthCov} --toolcov-cutoff ${params.toolCov} \
		--beddir ${perfDir} \
		--enable-cache --using-cache -o .

	echo "### Site level analysis DONE"
    """
}


