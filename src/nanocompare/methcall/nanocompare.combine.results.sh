#!/bin/bash
#SBATCH --job-name=combine.results
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=200g # memory pool for all cores
#SBATCH --time=13:00:00 # time (D-HH:MM:SS)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR

################################################################################
# Combine nanopore methylation call results and postpone processing of the final single file
# Need to populate the parameters into this script
# Output file name is "<dsname>.<tool>.*.combine.tsv"
################################################################################
cd "$(dirname "$0")"

set -e

set +x
#source /home/liuya/.bash_profile
#export PATH=/cm/shared/apps/slurm/18.08.8/bin:${PATH}
source utils.common.sh

set -x

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "analysisPrefix: ${analysisPrefix}"
echo "methCallsDir: ${methCallsDir}"
echo "clusterDeepModModel: ${clusterDeepModModel}"
echo "chromSizesFile: ${chromSizesFile}"
echo "run_clean: ${run_clean}"
echo "##################"
set +u

filter_taskid=

if [ "${Tool}" = "Tombo" ] ; then
	set +x
	conda activate nanoai
	set -x
	for fn in $(ls $methCallsDir/$analysisPrefix.batch_*.*.tombo.per_read_stats); do
		## Postprocess per read methylation calls for each run complete
		time python Tombo_extract_per_read_stats.py \
				${chromSizesFile} ${fn} ${fn}.bed
		wc -l ${fn}.bed
	done
	set +x
	conda deactivate
	set -x
	echo "###   Tombo extract per-read-stats bed files DONE"

	ls ${methCallsDir}/${analysisPrefix}.batch_*.*.tombo.per_read_stats.bed | wc -l

	> ${methCallsDir}/${dsname}.tombo.perReadsStats.combined.tsv

	cat ${methCallsDir}/${analysisPrefix}.batch_*.*.tombo.per_read_stats.bed > ${methCallsDir}/${dsname}.tombo.perReadsStats.combined.tsv

	wc -l ${methCallsDir}/${dsname}.tombo.perReadsStats.combined.tsv

	echo "### Tombo combining batches results finished. ###"

	# we need do filter out non-CG patterns also, now using MPI version
	filter_ret=$(sbatch --nodes=1 --ntasks=60 --job-name=flt.noCG.${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/meth_stats_tool_mpi.sh tombo-add-seq -i ${methCallsDir}/${dsname}.tombo.perReadsStats.combined.tsv --mpi --processors 60 -o ${methCallsDir})

	filter_taskid=:$(echo ${filter_ret} |grep -Eo '[0-9]+$')

	# wait the filter task finished, then start clean, so block here
#	srun --dependency=afterok:${filter_taskid} echo "Wait Filter-out NON-CG task finished."

	echo "### Tombo filter out NON-CG patterns post-process jobs submitted. ###"
#	echo "### DeepMod combine results DONE. ###"

elif [ "${Tool}" = "DeepSignal" ] ; then
	ls ${methCallsDir}/*.deepsignal.call_mods.tsv | wc -l

	> ${methCallsDir}/${dsname}.deepsignal.call_mods.combine.tsv

	cat ${methCallsDir}/*.deepsignal.call_mods.tsv > ${methCallsDir}/${dsname}.deepsignal.call_mods.combine.tsv

	wc -l ${methCallsDir}/${dsname}.deepsignal.call_mods.combine.tsv
	echo "### DeepSignal combine results DONE. ###"

elif [ "${Tool}" = "DeepMod" ] ; then

	## Extract read-level DeepMod results

	filter_ret=$(sbatch --nodes=1 --ntasks=50 --job-name=extr.rldepmd.${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err ../meth_stats/meth_stats_tool_mpi.sh deepmod-read-level --processors 50 --basecallDir ${basecallOutputDir} --methcallDir ${methCallsDir} -o ${methCallsDir}/${dsname}.deepmod_read_level.read_level.combine.tsv --o2 ${methCallsDir}/${dsname}.deepmod_read_level.base_level.combined.C.bed)
	filter_taskid=":$(echo ${filter_ret} |grep -Eo '[0-9]+$')"

	## Step:  join results from different batches, based on ref: https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md
	## We need firstly use DeepMod script to merge different runs of modification detection
	cd ${methCallsDir}

	# Step: Detect modifications from FAST5 files, ref https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md#1-how-to-detect-modifications-from-fast5-files
	# Usage: python {} pred_folder-of-DeepMod Base-of-interest unique-fileid-in-sum-file [chr-list]
#	time python /projects/li-lab/yang/tools/DeepMod/tools/sum_chr_mod.py ${methCallsDir}/ C ${dsname}.C
	time python deepmod_sum_chr_mod.py ${methCallsDir}/ C ${dsname}.deepmod


	## Step: Output C in CpG motifs in a genome, i.e. CpG index in a human genome: (must be done only once per genome)
	## TODO: check why only once, generate common file for all dataset used? CHeck results firstly running. DeepMod N70
	## Must firstly generate these files to a folder like:/projects/li-lab/yang/workspace/nano-compare/data/genome_motif/C
	## No need any modifications later, I failed to generate with correct log, so use WR results instead.
	#	python /projects/li-lab/yang/tools/DeepMod/tools/generate_motif_pos.py ${refGenome} /projects/li-lab/yang/workspace/nano-compare/data/genome_motif/C C CG 0

	set +x
#	source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh
	conda activate nanoai
	set -x
	## Step: Generated clustered results to consider cluster effect, ref: https://github.com/WGLab/DeepMod/blob/master/docs/Usage.md#generated-clustered-results
	## Need nanoai env, using tf 1.8.0, or will be some compilation error
	## $1-sys.argv[1]+'.%s.C.bed', (save to sys.argv[1]+'_clusterCpG.%s.C.bed'),
	## $2-gmotfolder ('%s/motif_%s_C.bed')   $3-not used
	time python /projects/li-lab/yang/tools/DeepMod/tools/hm_cluster_predict.py ${methCallsDir}/${dsname}.deepmod /projects/li-lab/yang/workspace/nano-compare/data/genome_motif/C ${clusterDeepModModel}

	set +x
	conda deactivate
	set -x

	> ${dsname}.deepmod.C.combined.tsv

	for f in $(ls -1 ${dsname}.deepmod.chr*.C.bed)
	do
	  cat $f >> ${dsname}.deepmod.C.combined.tsv
	done

	wc -l ${dsname}.deepmod.C.combined.tsv

	> ${dsname}.deepmod.C_clusterCpG.combine.tsv

	for f in $(ls -1 ${dsname}.deepmod_clusterCpG.chr*.C.bed)
	do
	  cat $f >> ${dsname}.deepmod.C_clusterCpG.combine.tsv
	done

	wc -l ${dsname}.deepmod.C_clusterCpG.combine.tsv

	echo "### DeepMod combine all batches results. ###"

#	bash /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/meth_stats_tool_array_job_pipe.sh deepmod-add-seq 100 ${methCallsDir}/${dsname}.deepmod.C.combined.tsv ${methCallsDir}
	# we need do filter out non-CG patterns also, now using MPI version
	filter_ret=$(sbatch --nodes=1 --ntasks=60 --job-name=flt.noCG.${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err ../meth_stats/meth_stats_tool_mpi.sh deepmod-add-seq -i ${methCallsDir}/${dsname}.deepmod.C.combined.tsv --mpi --processors 60 -o ${methCallsDir})

	filter_taskid=${filter_taskid}:$(echo ${filter_ret} |grep -Eo '[0-9]+$')
	# wait the filter task finished, then start clean, so block here
#	srun --dependency=afterok:${filter_taskid} date
	echo "### DeepMod filter out NON-CG patterns post-process jobs submitted. ###"
#	echo "### DeepMod combine results DONE. ###"

elif [ "${Tool}" = "Nanopolish" ] ; then
	### TODO: how to cut first line: sed '1d' test.txt > result.txt
	ls ${methCallsDir}/*.nanopolish.methylation_calls.tsv | wc -l

	sed -n '1p' ${methCallsDir}/*.batch_1.nanopolish.methylation_calls.tsv > ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv

	for fn in ${methCallsDir}/*.nanopolish.methylation_calls.tsv; do
		echo "Append ${fn}"
		sed '1d' ${fn} >> ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv
	done

	wc -l ${methCallsDir}/${dsname}.nanopolish.methylation_calls.combine.tsv
	echo "### Nanopolish combine results DONE. ###"
else
	echo "Tool=${Tool} is not support now"
	exit -1
fi

# After combine task finished, or just filer-out task running, we will start task of clean, it will depend on filter-out task also if needed


################################################################################
# Step 5: Clean intermediate dirs
################################################################################
if [ "${run_clean}" = true ] ; then
	echo "Step5: clean dirs"

	# If previous step need to depend on
	depend_param=""
	if [ -n "${filter_taskid}" ] ; then
		depend_param="afterok${filter_taskid}"
	fi
	sbatch --job-name=clen.${analysisPrefix} --output=${methCallsDir}/log/%x.%j.out --error=${methCallsDir}/log/%x.%j.err --dependency=${depend_param} --export=dsname=${dsname},Tool=${Tool},targetNum=${targetNum},analysisPrefix=${analysisPrefix},untaredInputDir=${untaredInputDir},septInputDir=${septInputDir},basecallOutputDir=${basecallOutputDir},methCallsDir=${methCallsDir},outbasedir=${outbasedir},clean_preprocessing=${clean_preprocessing},clean_basecall=${clean_basecall},tar_basecall=${tar_basecall},tar_methcall=${tar_methcall},run_clean=${run_clean} nanocompare.clean.intermediate.sh
	echo "Submitted clean dirs task for ${analysisPrefix}."
fi
