#!/usr/bin/env nextflow

params.dsname = 'DEMODATA'
params.indir = '/projects/li-lab/yang/workspace/nano-compare/data/raw-fast5/demo.fast5.reads.tar.gz'
params.ntask = 10

params.analysisPrefix=${params.dsname}-N${params.ntask}

params.isGPU=no
params.processors = 8

params.correctedGroup="RawGenomeCorrected_000"
params.refGenome="/projects/li-lab/reference/hg38/hg38.fasta"

params.deepsignalModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt"


params.guppy.basecall='/projects/li-lab/software/ont-guppy-gpu_4.2.2/bin/guppy_basecaller'



process Preprocess {
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

/*
[/pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/1,
/pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/10,
/pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/2, /pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/3, /pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/4, /pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/5, /pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/6, /pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/7, /pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/8, /pod/2/li-lab/yang/workspace/nano-compare/nf/work/10/40a474166d37b7e9128c8767eea9bd/sept_dir/9]
*/

process Basecall {
	input: /* TODO: how to get the exact folder name and number of task id of previous process here? */
    file x from fast5Inputs.flatten()

    output:
    file 'basecall_dir/*' into basecallOutputs

	/*
	TODO: how to use GPU and slurm task on it
	*/
    """
    echo hello
    mkdir -p basecall_dir
    echo ${x}

    echo $it
    ${params.guppy.basecall} --input_path $x \
        --save_path basecall_dir --config dna_r9.4.1_450bps_hac.cfg \
        --gpu_runners_per_device ${params.processors} --num_callers 3 --fast5_out --verbose_logs --device auto

    """
}


process Resquiggle {
	input: /* TODO: how to get the exact folder name of previous process here? */
    file x from basecallOutputs.flatten()

    output:
    file 'resquiggle_dir/*' into resquiggleOutputs

    """
	# TODO: how to get $x folder and num of task
	mkdir -p resquiggle_dir/numofx
    cp $x resquiggle_dir/numofx

    tombo resquiggle --dna --processes ${params.processors} --corrected-group ${params.correctedGroup} \
		--basecall-group Basecall_1D_000 --overwrite $x ${params.refGenome}
    """
}

process DeepSignal {
	input: /* TODO: how to get the exact folder name of previous process here? */
    file x from resquiggleOutputs.flatten()

    output:
    file '*.tsv' into deepsignalOutput

    """
	# TODO: how to get $x folder and num of task
	deepsignal call_mods --input_path $x --model_path ${params.deepsignalModel} \
		--result_file ${params.analysisPrefix}-DeepSignal.batch_${numofx}.CpG.deepsignal.call_mods.tsv \
		--reference_path ${params.refGenome} --corrected_group ${params.correctedGroup} --nproc ${params.processors} --is_gpu ${params.isGPU}

    """
}




