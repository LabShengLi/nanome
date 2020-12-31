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

	cat ${methCallsDir}/*perReadsStats.bed > ${methCallsDir}/${analysisPrefix}.tombo.perReadsStats.combine.bed

	wc -l ${methCallsDir}/${analysisPrefix}.tombo.perReadsStats.combine.bed
	echo "### Tombo combine results DONE. ###"
fi

if [ "${Tool}" = "DeepSignal" ] ; then
	ls ${methCallsDir}/*.deepsignal.call_mods.tsv | wc -l

	cat ${methCallsDir}/*.deepsignal.call_mods.tsv > ${methCallsDir}/${analysisPrefix}.deepsignal.call_mods.combine.tsv

	wc -l ${methCallsDir}/${analysisPrefix}.deepsignal.call_mods.combine.tsv
	echo "### DeepSignal combine results DONE. ###"
fi

if [ "${Tool}" = "DeepMod" ] ; then
	## Step:  join results from different batches
	cd ${methCallsDir}
	python /projects/li-lab/yang/tools/DeepMod/tools/sum_chr_mod.py ${methCallsDir}/ C ${dsname}.C

	## Step: to consider cluster effect
	python /projects/li-lab/yang/tools/DeepMod/tools/hm_cluster_predict.py ${methCallsDir}/${dsname}.C genome_motif/C ${clusterDeepModModel}

	for f in $(ls -1 ${dsname}.C.chr*)
	do
	  cat $f >> ${dsname}.deepmod.C.combined.bed
	done

	for f in $(ls -1 ${dsname}.C_clusterCpG.chr*)
	do
	  cat $f >> ${dsname}.deepmod.C_clusterCpG.combined.bed
	done

	echo "### DeepMod combine results DONE. ###"
fi

if [ "${Tool}" = "Nanopolish" ] ; then
	### TODO: how to cut first line: sed '1d' test.txt > result.txt
	ls ${methCallsDir}/*.nanopolish.methylation_calls.tsv | wc -l

	rm -rf ${methCallsDir}/${analysisPrefix}.nanopolish.methylation_calls.combine.tsv

	sed -n '1p' ${methCallsDir}/*.batch_0.nanopolish.methylation_calls.tsv > ${methCallsDir}/${analysisPrefix}.nanopolish.methylation_calls.combine.tsv

	for fn in ${methCallsDir}/*.nanopolish.methylation_calls.tsv; do
		echo "Append ${fn}"
		sed '1d' ${fn} >> ${methCallsDir}/${analysisPrefix}.nanopolish.methylation_calls.combine.tsv
	done
#	cat ${methCallsDir}/*.nanopolish.methylation_calls.tsv > ${methCallsDir}/${analysisPrefix}.nanopolish.methylation_calls.combine.tsv

	wc -l ${methCallsDir}/${analysisPrefix}.nanopolish.methylation_calls.combine.tsv
	echo "### Nanopolish combine results DONE. ###"
fi