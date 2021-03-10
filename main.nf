#!/usr/bin/env nextflow

Channel
    .fromPath(params.indir)
    .ifEmpty { exit 1, "Cannot find input file"}
    .set {ch_input}

Channel
    .from( 'a', 'b', 'aa', 'bc', 3, 4.5 )
    .filter( ~/^a.*/ )
    .view()


process Preprocess {
    input:
    set file(fast5_tar) from ch_input

	output:
    file 'sept_dir/M*' into fast5Inputs

    """
    mkdir -p untar_dir
    tar xzf ${fast5_tar} -C untar_dir

    mkdir -p sept_dir
    python ${workflow.projectDir}/src/nanocompare/methcall/FilesSeparatorNew.py untar_dir ${params.ntask} sept_dir
    """
}

//fast5Inputs.into { fast5Inputs; fast5Inputs2 }
//fast5Inputs2.view()
//fast5Inputs.subscribe { println it }

process Basecall {
	input:
    file x from fast5Inputs.flatten()

    output:
    file 'basecall_dir/M*' into basecallOutputs

	/*
	TODO: how to use GPU and slurm task on this single process
	*/
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

//basecallOutputs.subscribe { println it }

process Resquiggle {
	input:
    file x from basecallOutputs.flatten()

    output:
    file 'resquiggle_dir/M*' into resquiggleOutputs

    """
	mkdir -p resquiggle_dir/${x}
	cp -rf ${x}/* resquiggle_dir/${x}

    tombo resquiggle --dna --processes ${params.processors} \
        --corrected-group ${params.correctedGroup} \
		--basecall-group Basecall_1D_000 --overwrite \
		resquiggle_dir/${x} ${params.refGenome}
    """
}

//resquiggleOutputs.subscribe { println it }

process DeepSignal {
	input:
    file x from resquiggleOutputs.flatten()

    output:
    file '*.tsv' into deepsignalOutput

    """
	deepsignal call_mods --input_path $x \
	    --model_path ${params.deepsignalModel} \
		--result_file "${params.dsname}-N${params.ntask}-DeepSignal.batch_${x}.CpG.deepsignal.call_mods.tsv" \
		--reference_path ${params.refGenome} \
		--corrected_group ${params.correctedGroup} \
		--nproc ${params.processors} \
		--is_gpu ${params.isGPU}
    """
}

deepsignalOutput.view()

/*
process DeepSignalPostprocess {
	input:
    file x from deepsignalOutput

    output:
    file 'combine.tsv' into deepsignalPostOutput

    """
	echo ${x}
	touch combine.tsv
    """
}
*/


//deepsignalOutput.subscribe { println it }


