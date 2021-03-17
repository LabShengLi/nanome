#!/usr/bin/env nextflow

Channel
    .fromPath(params.indata)
    .ifEmpty { exit 1, "Cannot find input file"}
    .set {ch_input}

// untar file, seperate into N folders named 'M1', ..., 'M10', etc.
process Preprocess {
	// TODO: how to change conda env name here, or in whole nextflow pipeline?
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'

    input:
    set file(fast5_tar) from ch_input

	output:
    file 'sept_dir/M*' into fast5Inputs

    //when:
    //params.RunDeepMod == "true"

    """
    set -x
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


// fast5Inputs1ForBasecall.close()
// fast5Inputs2ForMegalodon.close()

// basecall of subfolders named 'M1', ..., 'M10', etc.
process Basecall {
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'

	input:
    file x from fast5Inputs1ForBasecall.flatten()

    output:
    file 'basecall_dir/M*' into basecallOutputs

    """
    mkdir -p basecall_dir/${x}
    ${params.guppy.basecall} --input_path $x \
        --save_path "basecall_dir/${x}" \
        --config dna_r9.4.1_450bps_hac.cfg \
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
	conda ${params.condaEnvName}
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'
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
	conda '/home/liuya/anaconda3/envs/nanoai'
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'
	input:
    file x from resquiggleOutputs1ForDeepSignal.flatten()

    output:
    file '*.tsv' into deepsignalOutput

    """
	deepsignal call_mods --input_path ${x}/workspace \
	    --model_path ${params.deepsignalModel} \
		--result_file "${params.dsname}-N${params.ntask}-DeepSignal.batch_${x}.CpG.deepsignal.call_mods.tsv" \
		--reference_path ${params.refGenome} \
		--corrected_group ${params.correctedGroup} \
		--nproc ${params.processors} \
		--is_gpu ${params.isGPU}
    """
}


// Tombo runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Tombo {
	conda '/home/liuya/anaconda3/envs/nanoai'
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'
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

	python ${workflow.projectDir}/src/nanocompare/methcall/Tombo_extract_per_read_stats.py ${params.chromSizesFile} \
				"${params.analysisPrefix}.batch_${x}.CpG.tombo.per_read_stats" \
				"${params.analysisPrefix}.batch_${x}.CpG.tombo.per_read_stats.bed"
    """
}


// Megalodon runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Megalodon {
	conda '/home/liuya/anaconda3/envs/nanoai'
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'
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
	conda '/home/liuya/anaconda3/envs/nanoai'
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'
	input:
    file x from basecallOutputs3ForDeepMod.flatten()

    output:
    file 'mod_output/batch_*_num' into deepmodOutput

    when:
    params.RunDeepMod == "true"

    """
    python ${params.DeepModDir}/bin/DeepMod.py detect \
			--wrkBase ${x}/workspace --Ref ${params.refGenome} \
			--Base C --modfile ${params.deepModModel} \
			--FileID batch_${x}_num \
			--threads ${params.processors} --move

	# rm -rf mod_output/batch_*.done
    """
}




// Nanopolish runs on resquiggled subfolders named 'M1', ..., 'M10', etc.
process Nanopolish {
	conda '/home/liuya/anaconda3/envs/nanoai'
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'

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

	${params.NanopolishDir}/nanopolish index -d ${x}/workspace \${fastqNoDupFile}

	minimap2 -t ${params.processors} -a -x map-ont ${params.refGenome} \${fastqNoDupFile} | samtools sort -T tmp -o ${x}/\${bamFileName}
	echo "### minimap2 finished"

	samtools index -@ threads ${x}/\${bamFileName}
	echo "### samtools finished"

	echo "### Alignment step DONE"

	${params.NanopolishDir}/nanopolish call-methylation -t ${params.processors} -r \${fastqNoDupFile} -b ${x}/\${bamFileName} -g ${params.refGenome} > ${params.analysisPrefix}.batch_${x}.nanopolish.methylation_calls.tsv

	# touch batch_${x}.nanopolish.tsv
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
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'

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
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'

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
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'

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
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'

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
	conda '/home/liuya/anaconda3/envs/nanoai'
	executor 'slurm'
	clusterOptions '-p gpu -q inference -n 8 --gres=gpu:1 --time=06:00:00 --mem-per-cpu=170G'

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

    python ${params.DeepModDir}/DeepMod_tools/sum_chr_mod.py \
        indir/ C ${params.dsname}.deepmod

	python ${params.DeepModDir}/DeepMod_tools/hm_cluster_predict.py \
		indir/${params.dsname}.deepmod \
		${params.genomeMotifC} \
		${params.clusterDeepModModel}  || true

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


//allCombineResults=deepsignalCombineResult.concat(tomboCombineResult, megalodonCombineResult, nanopolishCombineResult, deepmodCombineResult)
//allCombineResults.view()

