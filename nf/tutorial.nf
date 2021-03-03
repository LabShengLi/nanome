#!/usr/bin/env nextflow

params.str = 'Hello world!'

params.indir = '/projects/li-lab/yang/workspace/nano-compare/data/raw-fast5/demo.fast5.reads.tar.gz'
params.ntask = 10

process preProcess {
	output:
    file 'sept_dir/*' into fast5Inputs

    """
    mkdir -p untar_dir
    mkdir -p sept_dir
    tar xzvf ${params.indir} -C untar_dir
    python /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/FilesSeparatorNew.py untar_dir ${params.ntask} sept_dir
    """
}

fast5Inputs.subscribe { println it }

process baseCall {
	input:
    file x from fast5Inputs.flatten()

    output:
    file 'basecall_dir/*' into basecallOutputs

    """
    echo hello
    mkdir -p basecall_dir
    echo ${x}

    echo $it
    /projects/li-lab/software/ont-guppy-gpu_4.2.2/bin/guppy_basecaller --input_path $x \
        --save_path basecall_dir --config dna_r9.4.1_450bps_hac.cfg \
        --gpu_runners_per_device 8 --num_callers 3 --fast5_out --verbose_logs --device auto

    """
}


