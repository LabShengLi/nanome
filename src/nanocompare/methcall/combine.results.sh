#!/bin/bash
#SBATCH --job-name=combine.tombo
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 1-23:00:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

################################################################################
# Combine nanopore methylation call results to a single file
# Need to populate the parameters into this script
################################################################################

set -x

if [ "${Tool}" = "Tombo" ] ; then
	ls ${methCallsDir}/*perReadsStats.bed | wc -l

	cat ${methCallsDir}/*perReadsStats.bed > ${methCallsDir}/${dsname}.tombo.perReadsStats.combine.bed

	wc -l ${methCallsDir}/${dsname}.tombo.perReadsStats.combine.bed
	echo "### Tombo combine results DONE. ###"
fi

if [ "${Tool}" = "DeepSignal" ] ; then
	ls ${methCallsDir}/*.deepsignal.call_mods.tsv | wc -l

	cat ${methCallsDir}/*.deepsignal.call_mods.tsv > ${methCallsDir}/${dsname}.deepsignal.call_mods.combine.tsv

	wc -l ${methCallsDir}/${dsname}.deepsignal.call_mods.combine.tsv
	echo "### DeepSignal combine results DONE. ###"
fi

if [ "${Tool}" = "DeepMod" ] ; then
	## Step:  join results from different batches
	cd ${methCallsDir}
	python /projects/li-lab/yang/tools/DeepMod/tools/sum_chr_mod.py ${methCallsDir}/ C ${dsname}.C

	## Step: CpG index in a human genome: (must be done only once per genome)
	## TODO: check why only once, generate common file for all dataset used? CHeck results firstly running. DeepMod N70
	## Must firstly generate these files to a folder like:/projects/li-lab/yang/workspace/nano-compare/data/genome_motif/C
	## No need any modifications later
	python /projects/li-lab/yang/tools/DeepMod/tools/generate_motif_pos.py ${refGenome} /projects/li-lab/yang/workspace/nano-compare/data/genome_motif/C C CG 0

	set +x
	source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh
	conda activate nanoai
	set -x
	## Step: to consider cluster effect
	## Need nanoai env, using tf 1.8.0, or will be some compilation error
	python /projects/li-lab/yang/tools/DeepMod/tools/hm_cluster_predict.py ${methCallsDir}/${dsname}.C /projects/li-lab/yang/workspace/nano-compare/data/genome_motif/C ${clusterDeepModModel}

	set +x
	conda deactivate
	set -x

	rm -rf ${dsname}.deepmod.C.combined.bed
	touch ${dsname}.deepmod.C.combined.bed

	for f in $(ls -1 ${dsname}.C.chr*)
	do
	  cat $f >> ${dsname}.deepmod.C.combined.bed
	done

	wc -l ${dsname}.deepmod.C.combined.bed

	rm -rf ${dsname}.deepmod.C_clusterCpG.combined.bed
	touch ${dsname}.deepmod.C_clusterCpG.combined.bed

	for f in $(ls -1 ${dsname}.C_clusterCpG.chr*)
	do
	  cat $f >> ${dsname}.deepmod.C_clusterCpG.combined.bed
	done

	wc -l ${dsname}.deepmod.C_clusterCpG.combined.bed

	echo "### DeepMod combine results DONE. ###"
fi

if [ "${Tool}" = "Nanopolish" ] ; then
	### TODO: how to cut first line: sed '1d' test.txt > result.txt
	ls ${methCallsDir}/*.nanopolish.methylation_calls.tsv | wc -l

	rm -rf ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv

	sed -n '1p' ${methCallsDir}/*.batch_0.nanopolish.methylation_calls.tsv > ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv

	for fn in ${methCallsDir}/*.nanopolish.methylation_calls.tsv; do
		echo "Append ${fn}"
		sed '1d' ${fn} >> ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv
	done

	wc -l ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv
	echo "### Nanopolish combine results DONE. ###"
fi