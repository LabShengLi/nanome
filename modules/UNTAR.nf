/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : UNTAR.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// Untar of subfolders named 'M1', ..., 'M10', etc.
process UNTAR {
	tag "${fast5Input.baseName}"

	input:
	path fast5Input // can be a folder/tar/tar.gz file

	output:
	// untar dir or basecalled dir structure
	path "${fast5Input.baseName}.untar", emit: untar,  optional: true
	tuple val(fast5Input.baseName), path ("${fast5Input.baseName}.untar"),	optional:true,  emit: untar_tuple

	shell:
	cores = task.cpus * params.highProcTimes
	if (!params.skipBasecall) { // normal case: need basecall later
		'''
		date; hostname; pwd
		echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"

		## Extract input files tar/tar.gz/folder
		mkdir -p untarTempDir
		if [[ !{fast5Input} == *.tar && -f !{fast5Input} ]] ; then
			### deal with tar
			tar -xf !{fast5Input} -C untarTempDir
		elif [[ !{fast5Input} == *.tar.gz && -f !{fast5Input} ]] ; then
			### deal with tar.gz
			tar -xzf !{fast5Input} -C untarTempDir
		elif [[ -d !{fast5Input} ]]; then
			## Copy files, do not change original files such as old analyses data
			find !{fast5Input}/ -name '*.fast5' | \
				parallel -j!{cores}  cp {} untarTempDir/
		else
			echo "### Untar error for input=!{fast5Input}"
		fi

		if [[ !{params.multi_to_single_fast5} == true ]] ; then
			echo "### Do multi_to_single_fast5"
			untarTempDir=untarTempDir2
			mkdir -p ${untarTempDir}
			multi_to_single_fast5 -i untarTempDir -s untarTempDir2 -t !{cores}  --recursive
		else
			untarTempDir=untarTempDir
		fi

		## Move fast5 raw/basecalled files into XXX.untar folder
		mkdir -p !{fast5Input.baseName}.untar

		find ${untarTempDir} -name "*.fast5" -type f | \
			parallel -j!{cores}  mv {}  !{fast5Input.baseName}.untar/

		## Clean temp files
		rm -rf untarTempDir  untarTempDir2

		## Clean old basecalled analyses in input fast5 files
		if [[ "!{params.cleanAnalyses}" == true ]] ; then
			echo "### Start cleaning old analysis"
			## python -c 'import h5py; print(h5py.version.info)'
			clean_old_basecall_in_fast5.py \
				-i !{fast5Input.baseName}.untar --is-indir --verbose\
				--processor !{cores}
		fi

		totalFiles=$( find !{fast5Input.baseName}.untar -name "*.fast5" -type f | wc -l )
		echo "### Total fast5 input files:${totalFiles}"
		if (( totalFiles==0 )); then
			echo "### no fast5 files at !{fast5Input.baseName}.untar, skip this job"
			rm -rf !{fast5Input.baseName}.untar
		fi
		echo "### Untar DONE"
		'''
	} else { // deal with basecalled input data
		'''
		date; hostname; pwd

		## Extract input files tar/tar.gz/folder
		mkdir -p untarTempDir
		if [[ !{fast5Input} == *.tar && -f !{fast5Input} ]] ; then
			### deal with tar
			tar -xf !{fast5Input} -C untarTempDir
		elif [[ !{fast5Input} == *.tar.gz && -f !{fast5Input} ]] ; then
			### deal with tar.gz
			tar -xzf !{fast5Input} -C untarTempDir
		elif [[ -d !{fast5Input} ]] ; then
			## user provide basecalled input dir, just cp them
			mkdir -p untarTempDir/test
			cp -rf !{fast5Input}/*   untarTempDir/test/
		else
			echo "### Untar error for input=!{fast5Input}"
		fi

		## Move fast5 raw/basecalled files into XXX.untar folder
		mkdir -p !{fast5Input.baseName}.untar
		## Keep the directory structure for basecalled input
		mv untarTempDir/*/*   !{fast5Input.baseName}.untar/

		## Clean temp files
		rm -rf untarTempDir

		totalFiles=$( find !{fast5Input.baseName}.untar -name "*.fast5" -type f | wc -l )
		echo "### Total fast5 input files:${totalFiles}"
		if (( totalFiles==0 )); then
			echo "### no fast5 files at !{fast5Input.baseName}.untar, skip this job"
			rm -rf !{fast5Input.baseName}.untar
		fi
		echo "### Untar DONE"
		'''
	}
}
