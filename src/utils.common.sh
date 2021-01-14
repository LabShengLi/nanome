#!/bin/bash

# Specify your home bash, this required by conda setup
HOME_BASH_FN=/home/liuya/.bash_profile


set +u
set +e
source ${HOME_BASH_FN}

#set -x

## Export slurm path
export PATH=/cm/shared/apps/slurm/18.08.8/bin:${PATH}

# Export DepMod and Nanopolish dir
export NanopolishDir=/projects/li-lab/yang/tools/latest-version/nanopolish
export DeepModDir=/projects/li-lab/yang/tools/latest-version/DeepMod
export GuppyDir=/projects/li-lab/software/ont-guppy-gpu_4.2.2


# Specify global vars for nano-compare project
### Reference file path configuration, used by each base or meth calling
correctedGroup="RawGenomeCorrected_000"
refGenome="/projects/li-lab/reference/hg38/hg38.fasta"
chromSizesFile="${NanoCompareDir}/data/genome-annotation/hg38.chrom.sizes"

#deepsignalModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/model.CpG.R9.4_1D.human_hx1.bn17.sn360/bn_17.sn_360.epoch_7.ckpt"
deepsignalModel="${NanoCompareDir}/data/dl-model/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt"

#deepModModel="/projects/li-lab/yang/workspace/nano-compare/data/dl-model/rnn_conmodC_P100wd21_f7ne1u0_4/mod_train_conmodC_P100wd21_f3ne1u0"
deepModModel="${DeepModDir}/train_deepmod/rnn_conmodC_P100wd21_f7ne1u0_4train_mod/rnn_conmodC_P100wd21_f7ne1u0_4/mod_train_conmodC_P100wd21_f3ne1u0"

clusterDeepModModel="${DeepModDir}/train_deepmod/na12878_cluster_train_mod-keep_prob0.7-nb25-chr1/Cg.cov5.nb25"


get_arrayjob_ids(){
	## Generate dependency pattern format from sbatch job submission return string
	## Input is <ret-str> <num>
	## Ouput is like: ":12345_1:12345_2:12345_3:12345_4:12345_5:12345_6:12345_7:12345_8:12345_9:12345_10"
	## Example:
	## ret='Submitted batch job 5812809'
	## num=10
	## task_ids=$(get_arrayjob_ids "${ret}" "${num}")
	## Later can be used by "--dependency=afterok${task_ids}"
	usage="get_arrayjob_ids '<ret-str>' '<num>'"
	local ret=${1:?"undefined '<ret-str>': $usage"}
	local num=${2:?"undefined '<num>': $usage"}

	local arrayjob_id=${ret##* }
	local taskids=

	set +x
	for ((i=1;i<=${num};i++));
	do
		taskids=${taskids}:${arrayjob_id}_$i
	done
	echo ${taskids}
	return 0
}
