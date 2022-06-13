/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : BASECALL.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// basecall of subfolders named 'M1', ..., 'M10', etc.
process BASECALL {
	tag "${fast5Untar.baseName}"

	input:
	path fast5Untar

	output:
	path "${fast5Untar.baseName}.basecall", optional:true,	emit: basecall
	tuple val(fast5Untar.baseName), path ("${fast5Untar.baseName}.basecall"),	optional:true,  emit: basecall_tuple  // must use Name value, not file var, or will failed for B

	when:
	params.runBasecall

	shell:
	'''
	date; hostname; pwd

	echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"
	if [[ ${CUDA_VISIBLE_DEVICES:-} == "" ]] ; then
		echo "Detect no GPU, using CPU commandType"
		commandType='cpu'
		gpuOptions=" "
	else
		echo "Detect GPU, using GPU commandType"
		commandType='gpu'
		gpuOptions="-x auto"
	fi

	which guppy_basecaller
	guppy_basecaller -v
	mkdir -p !{fast5Untar.baseName}.basecall

	if [[ !{params.skipBasecall} == false ]] ; then
		## CPU/GPU version command
		guppy_basecaller --input_path !{fast5Untar} \
			--save_path  !{fast5Untar.baseName}.basecall \
			--config !{params.GUPPY_BASECALL_MODEL} \
			--num_callers !{task.cpus} \
			--fast5_out --compress_fastq\
			--verbose_logs  ${gpuOptions} &>> !{params.dsname}.!{fast5Untar.baseName}.Basecall.run.log
	else
		## Just use user's basecalled input
		cp -rf !{fast5Untar}/*   !{fast5Untar.baseName}.basecall/
	fi

	## Combine fastq
	touch !{fast5Untar.baseName}.basecall/batch_basecall_combine_fq_!{fast5Untar.baseName}.fq.gz

	## Below is compatable with both Guppy v4.2.2 (old) and newest directory structures
	## !{fast5Untar.baseName}.basecall/fail

	if [[ params.filter_fail_fq == true ]] ; then
		failDir=
	else
		failDir=!{fast5Untar.baseName}.basecall/fail
	fi

	find !{fast5Untar.baseName}.basecall/ !{fast5Untar.baseName}.basecall/pass \
	 	${failDir}  -maxdepth 1 -name '*.fastq.gz' -type f\
	 	-print0 2>/dev/null | \
	 	while IFS= read -r -d '' file;  do
	 		cat $file >> \
	 			!{fast5Untar.baseName}.basecall/batch_basecall_combine_fq_!{fast5Untar.baseName}.fq.gz
	 	done
	echo "### Combine fastq.gz DONE"

	## Remove fastq.gz
	find !{fast5Untar.baseName}.basecall/   !{fast5Untar.baseName}.basecall/pass/ \
	 	!{fast5Untar.baseName}.basecall/fail/ -maxdepth 1 -name '*.fastq.gz' -type f 2>/dev/null |\
	 	parallel -j!{task.cpus * params.highProcTimes} 'rm -f {}'

	## After basecall, rename and publish summary filenames, summary may also be used by resquiggle
	mv !{fast5Untar.baseName}.basecall/sequencing_summary.txt \
		!{fast5Untar.baseName}.basecall/!{fast5Untar.baseName}-sequencing_summary.txt

    ## Clean
    if [[ !{params.cleanStep} == "true" ]]; then
    	echo "### No need to clean"
    fi
	echo "### Basecalled by Guppy DONE"
	'''
}
