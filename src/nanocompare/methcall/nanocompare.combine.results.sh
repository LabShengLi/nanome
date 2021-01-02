#!/bin/bash
#SBATCH --job-name=combine.tombo
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=200g # memory pool for all cores
#SBATCH --time=13:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

################################################################################
# Combine nanopore methylation call results and postpone processing of the final single file
# Need to populate the parameters into this script
# Output file name is "<dsname>.<tool>.*.combine.tsv"
################################################################################
set -e
set -x

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "methCallsDir: ${methCallsDir}"
echo "clusterDeepModModel: ${clusterDeepModModel}"
echo "##################"
set +u

if [ "${Tool}" = "Tombo" ] ; then
	ls ${methCallsDir}/*perReadsStats.bed | wc -l

	> ${methCallsDir}/${dsname}.tombo.perReadsStats.combine.tsv

	cat ${methCallsDir}/*perReadsStats.bed > ${methCallsDir}/${dsname}.tombo.perReadsStats.combine.tsv

	wc -l ${methCallsDir}/${dsname}.tombo.perReadsStats.combine.tsv

	echo "### Tombo combine results DONE. ###"

	# we need do filter out non-CG patterns also
	bash /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/meth_stats_tool_array_job_pipe.sh tombo-add-seq 300 ${methCallsDir}/${dsname}.tombo.perReadsStats.combine.tsv ${methCallsDir}

	echo "### Tombo filter out NON-CG patterns post-process jobs submitted. ###"
fi

if [ "${Tool}" = "DeepSignal" ] ; then
	ls ${methCallsDir}/*.deepsignal.call_mods.tsv | wc -l

	> ${methCallsDir}/${dsname}.deepsignal.call_mods.combine.tsv

	cat ${methCallsDir}/*.deepsignal.call_mods.tsv > ${methCallsDir}/${dsname}.deepsignal.call_mods.combine.tsv

	wc -l ${methCallsDir}/${dsname}.deepsignal.call_mods.combine.tsv
	echo "### DeepSignal combine results DONE. ###"
fi

if [ "${Tool}" = "DeepMod" ] ; then
	## Step:  join results from different batches, based on ref: https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md
	## We need firstly use DeepMod script to merge different runs of modification detection
	cd ${methCallsDir}

	# Step: Detect modifications from FAST5 files, ref https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md#1-how-to-detect-modifications-from-fast5-files
	# Usage: python {} pred_folder-of-DeepMod Base-of-interest unique-fileid-in-sum-file [chr-list]
	python /projects/li-lab/yang/tools/DeepMod/tools/sum_chr_mod.py ${methCallsDir}/ C ${dsname}.C

	## Step: Output C in CpG motifs in a genome, i.e. CpG index in a human genome: (must be done only once per genome)
	## TODO: check why only once, generate common file for all dataset used? CHeck results firstly running. DeepMod N70
	## Must firstly generate these files to a folder like:/projects/li-lab/yang/workspace/nano-compare/data/genome_motif/C
	## No need any modifications later, I failed to generate with correct log, so use WR results instead.
	#	python /projects/li-lab/yang/tools/DeepMod/tools/generate_motif_pos.py ${refGenome} /projects/li-lab/yang/workspace/nano-compare/data/genome_motif/C C CG 0

	set +x
	source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh
	conda activate nanoai
	set -x
	## Step: Generated clustered results to consider cluster effect, ref: https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md#generated-clustered-results
	## Need nanoai env, using tf 1.8.0, or will be some compilation error
	## $1-sys.argv[1]+'.%s.C.bed', (save to sys.argv[1]+'_clusterCpG.%s.C.bed'),
	## $2-gmotfolder ('%s/motif_%s_C.bed')   $3-not used
	python /projects/li-lab/yang/tools/DeepMod/tools/hm_cluster_predict.py ${methCallsDir}/${dsname}.C /projects/li-lab/yang/workspace/nano-compare/data/genome_motif/C1 ${clusterDeepModModel}

	set +x
	conda deactivate
	set -x

	> ${dsname}.deepmod.C.combine.tsv

	for f in $(ls -1 ${dsname}.C.chr*.C.bed)
	do
	  cat $f >> ${dsname}.deepmod.C.combine.tsv
	done

	wc -l ${dsname}.deepmod.C.combine.tsv

	> ${dsname}.deepmod.C_clusterCpG.combine.tsv

	for f in $(ls -1 ${dsname}.C_clusterCpG.chr*.C.bed)
	do
	  cat $f >> ${dsname}.deepmod.C_clusterCpG.combine.tsv
	done

	wc -l ${dsname}.deepmod.C_clusterCpG.combine.tsv

	echo "### DeepMod combine results DONE. ###"

	bash /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/meth_stats_tool_array_job_pipe.sh deepmod-add-seq 100 ${methCallsDir}/${dsname}.deepmod.C.combine.tsv ${methCallsDir}

	echo "### DeepMod filter out NON-CG patterns post-process jobs submitted. ###"

fi

if [ "${Tool}" = "Nanopolish" ] ; then
	### TODO: how to cut first line: sed '1d' test.txt > result.txt
	ls ${methCallsDir}/*.nanopolish.methylation_calls.tsv | wc -l

	sed -n '1p' ${methCallsDir}/*.batch_0.nanopolish.methylation_calls.tsv > ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv

	for fn in ${methCallsDir}/*.nanopolish.methylation_calls.tsv; do
		echo "Append ${fn}"
		sed '1d' ${fn} >> ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv
	done

	wc -l ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv
	echo "### Nanopolish combine results DONE. ###"
fi