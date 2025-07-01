/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : RESQUIGGLE.nf
 @Software : NANOME project
 @Organization : Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// Resquiggle on basecalled subfolders named 'M1', ..., 'M10', etc.
process RESQUIGGLE {
	tag "${basecallDir.baseName}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Features-${params.dsname}",
		mode: "copy",
		pattern: "${basecallDir.baseName}.deepsignal1_batch_features.tsv.gz",
		enabled: params.feature_extract

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Features-${params.dsname}",
		mode: "symlink",
		pattern: "${basecallDir.baseName}.resquiggle",
		enabled: params.publishResquiggle

	publishDir "${params.outdir}/${params.dsname}-run-log",
		mode: "copy", pattern: "*.Resquiggle.run.log"

	input:
	tuple 	val(id), path (untarDir), path (basecallDir)
	each 	path(reference_genome)

	output:
	path "${basecallDir.baseName}.resquiggle", 	emit: resquiggle
	path "${basecallDir.baseName}.deepsignal1_batch_features.tsv.gz", 	emit: feature_extract, optional: true
	path "*.Resquiggle.run.log", optional:true,	emit: runlog


	when:
	(params.runMethcall && ((params.runDeepSignal1 && ! params.stopDeepSignal) || params.runTombo || params.runDeepSignal)) || params.runResquiggle

	shell:
	cores = task.cpus * params.highProcTimes
	samtools_cores = task.cpus * params.mediumProcTimes
	resquiggle_cores = (task.cpus*params.reduceProcTimes).intValue()
	'''
	### copy basecall workspace files, due to tombo resquiggle modify base folder
	rm -rf !{basecallDir.baseName}.resquiggle
	mkdir -p !{basecallDir.baseName}.resquiggle/workspace

	### original basecalled results will be parrallelly used by other processes
	cp -f !{basecallDir}/batch_basecall_combine_fq_*.fq.gz  \
		!{basecallDir.baseName}.resquiggle/

	## cp -rf !{untarDir}/*.fast5  !{basecallDir.baseName}.resquiggle/
	find !{untarDir}/ -name '*.fast5' -type f| \
		parallel -j!{task.cpus * params.highProcTimes}  \
		'cp {}   !{basecallDir.baseName}.resquiggle/workspace/'
	echo "### Duplicate from untar DONE"

	### Prerocessing, using combined fq.gz
	### ref: https://github.com/bioinfomaticsCSU/deepsignal#quick-start
	gunzip !{basecallDir.baseName}.resquiggle/batch_basecall_combine_fq_*.fq.gz
	tombo preprocess annotate_raw_with_fastqs\
		--fast5-basedir !{basecallDir.baseName}.resquiggle/workspace\
		--fastq-filenames !{basecallDir.baseName}.resquiggle/batch_basecall_combine_fq_*.fq\
		--sequencing-summary-filenames !{basecallDir}/!{untarDir.baseName}-sequencing_summary.txt \
		--basecall-group !{params.BasecallGroupName}\
		--basecall-subgroup !{params.BasecallSubGroupName}\
		--overwrite --processes  !{samtools_cores} \
		&>> !{params.dsname}.!{basecallDir.baseName}.Resquiggle.run.log
	echo "### tombo preprocess DONE"

	### Need to check Tombo resquiggle bugs, lots of users report long runtime and hang at nearly completion for large data
	### ref: https://github.com/nanoporetech/tombo/issues/139, https://github.com/nanoporetech/tombo/issues/111
	### ref: https://github.com/nanoporetech/tombo/issues/365, https://github.com/nanoporetech/tombo/issues/167
	### ref: https://nanoporetech.github.io/tombo/resquiggle.html?highlight=processes
	### Out of memory solution for large data, e.g. NA12878: --tomboResquiggleOptions '--signal-length-range 0 500000  --sequence-length-range 0 50000'
	tombo resquiggle\
		--processes !{resquiggle_cores} \
		--threads-per-process !{params.tomboThreadsPerProcess} \
		--corrected-group !{params.ResquiggleCorrectedGroup} \
		--basecall-group !{params.BasecallGroupName} \
		--basecall-subgroup !{params.BasecallSubGroupName}\
		--ignore-read-locks !{params.tomboResquiggleOptions ? params.tomboResquiggleOptions : ''}\
		--overwrite \
		!{basecallDir.baseName}.resquiggle/workspace \
		!{params.referenceGenome} &>> !{params.dsname}.!{basecallDir.baseName}.Resquiggle.run.log

	echo "### tombo resquiggle DONE"

	## Start to extract features
	if [[ !{params.feature_extract} == true ]] ; then
		deepsignal extract \
			--fast5_dir !{basecallDir.baseName}.resquiggle/workspace/ \
			--reference_path !{params.referenceGenome} \
			--write_path !{basecallDir.baseName}.deepsignal1_batch_features.tsv \
			--corrected_group !{params.ResquiggleCorrectedGroup} \
			--nproc !{cores}
		gzip -f !{basecallDir.baseName}.deepsignal1_batch_features.tsv
	fi
	'''
}