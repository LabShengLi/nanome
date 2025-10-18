/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : DORADO_UNTAR.nf, etc.
 @Software : NANOME project
 @Organization : Sheng Li Lab
----------------------------------------------------------------------------------------
*/
process DORADO_UNTAR {
	tag "${params.dsname}"

	input:
	path fast5InputList // can be a folder/tar/tar.gz file

	output:
	// untar dir
	path "${params.dsname}.untar", emit: untar,  optional: true

	when:
	params.dorado

	shell:
	cores = task.cpus * params.highProcTimes
	'''
		date; hostname; pwd

		echo fast5InputList=!{fast5InputList}

		## Extract input files tar/tar.gz/folder
		mkdir -p untarTempDir

		for fast5Input in !{fast5InputList}; do
			echo "Unpacking $fast5Input..."
			## fast5Input=!{fast5InputList}
			if [[ $fast5Input == *.tar && -f ${fast5Input} ]] ; then
				### deal with tar
				tar -xf ${fast5Input} -C untarTempDir
			elif [[ ${fast5Input} == *.tar.gz && -f ${fast5Input} ]] ; then
				### deal with tar.gz
				tar -xzf ${fast5Input} -C untarTempDir
			elif [[ -d ${fast5Input} ]]; then
				## For dir, should copy files, we do not want to change original files such as old analyses data in fast5
				find ${fast5Input}/ -name "*.!{params.file_format}"   | \
					parallel -j!{cores}  cp -L -f {} untarTempDir/
			else
				echo "### Untar error for input=${fast5Input}"
			fi
		done
		## Move fast5 raw/basecalled files into XXX.untar folder
		mkdir -p !{params.dsname}.untar

		find untarTempDir -name "*.!{params.file_format}" -type f | \
			parallel -j!{cores}  mv {}  !{params.dsname}.untar/

		## Clean temp files
		rm -rf untarTempDir

		echo "### Untar DONE"
	'''
}


process DORADO_CALL {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		pattern: "${params.dsname}.dorado_call",
		mode: "copy"

	input:
	path untar_dir // can be a folder/tar/tar.gz file
	path reference_genome

	output:
	// untar dir
	path "${params.dsname}.dorado_call", emit: dorado_call,  optional: true

	when:
	params.dorado

	shell:
	cores = task.cpus * params.highProcTimes
	kitOpt = params.demux ? "--kit-name ${params.kit_name}" : ""
	'''
		date; hostname; pwd

		echo untar_dir=!{untar_dir}
		echo reference_genome=!{reference_genome}
		echo !{kitOpt}
		dorado -vv

		ls /models

		mkdir -p !{params.dsname}.dorado_call

		dorado basecaller \
			!{params.dorado_model_dir}/!{params.dorado_basecall_model} \
			!{untar_dir}/ \
			!{kitOpt} \
			--models-directory !{params.dorado_model_dir} \
			--modified-bases-models !{params.dorado_model_dir}/!{params.dorado_methcall_model} \
			-x auto --verbose -r \
			--reference  !{params.referenceGenome} \
			--output-dir !{params.dsname}.dorado_call \
			--batchsize 96 2> >(tee -a !{params.dsname}.dorado_calls.run.log)

		mv !{params.dsname}.dorado_call/*.bam !{params.dsname}.dorado_call/!{params.dsname}.dorado_call.bam
		mv !{params.dsname}.dorado_call/*.bam.bai !{params.dsname}.dorado_call/!{params.dsname}.dorado_call.bam.bai
		echo "### Dorado call DONE"
	'''
}


process DORADO_DEMUX {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		pattern: "${params.dsname}.dorado_call_demux*",
		mode: "copy"

	input:
	path dorado_call // can be a folder/tar/tar.gz file

	output:
	// untar dir
	path "${params.dsname}.dorado_call_demux*", emit: dorado_demux,  optional: true

	when:
	params.dorado

	shell:
	cores = task.cpus * params.highProcTimes
	kitOpt = params.kit_name ? "--kit-name ${params.kit_name}" : ""
	'''
		date; hostname; pwd

		bam_fn=!{dorado_call}/*.bam

		mkdir -p !{params.dsname}.dorado_call_demux_!{params.kit_name}
		dorado demux \
			--output-dir !{params.dsname}.dorado_call_demux_!{params.kit_name} \
			--no-classify ${bam_fn} -v \
			2> >(tee -a !{params.dsname}_dorado_demux.run.log)
		echo "### Dorado demux DONE"
	'''
}


process DORADO_QC {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-DoradoQC",
		mode: "copy", overwrite: true

	input:
	path dorado_call
	path reference_genome

	output:
	path "${params.dsname}.dorado_qc",				emit: dorado_qc

	when:
	params.dorado

	shell:
	cores = task.cpus * params.highProcTimes
	samtools_cores = task.cpus * params.mediumProcTimes
	'''

	mkdir -p !{params.dsname}.dorado_qc

	bam_file=!{dorado_call}/*.bam
	NanoComp --bam ${bam_file} \
			--names !{params.dsname} --outdir !{params.dsname}.dorado_qc -t !{cores} \
			--raw  -f pdf -p !{params.dsname}_   &>> !{params.dsname}.DORADO_QC.run.log

	## extract QC metrics to tsv
	infer_input=(!{params.dsname}.dorado_qc/*_NanoStats.txt)

	input="${infer_input[0]}"
	output=!{params.dsname}.dorado_qc/!{params.dsname}_NanoStatsQC.tsv

	# Extract only the top 13 lines (the summary section)
	# Parse top 13 lines into one row (header + values)
	head -n 13 -- "$input" | \
	awk -F: -v output="$output" '
	{
	  key=$1
	  gsub(/^[ \t]+|[ \t]+$/, "", key)   # trim
	  gsub(/ /, "_", key)                # spaces -> _
	  val=$2
	  gsub(/,/, "", val)                 # remove commas in numbers
	  gsub(/^[ \t]+|[ \t]+$/, "", val)   # trim
	  keys   = keys   ? keys   "\t" key : key
	  values = values ? values "\t" val : val
	}
	END {
	  print keys  >  output   # header (overwrite)
	  print values >> output  # one row of values (append)
	}'
	echo "TSV file created: $output"

    echo "### Dorado QC all DONE"
	'''
}


// extract BAM as per-read results
process DORADO_CALL_EXTRACT {
	tag "${call_tagname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		pattern: "${params.dsname}.${call_tagname}.dorado_call_extract.tsv.gz",
		mode: "copy"

	input:
	val call_tagname // call's tagname
	// val bam_fn // BAM file location
	path dorado_call // can be a folder contains files of BAM and BAM.BAI
	path reference_genome
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}.${call_tagname}.dorado_call_extract.tsv.gz", emit: dorado_call_extract,  optional: true

	when:
	params.dorado

	shell:
	cores = task.cpus * params.highProcTimes
	'''
		date; hostname; pwd

		echo !{call_tagname}
		echo !{dorado_call}

		bam_fn=!{dorado_call}/*.bam

		python utils/modbam2bed_extract_read_cpg.py   \
				-r !{params.referenceGenome} \
				-i ${bam_fn} \
				-o !{params.dsname}.!{call_tagname}.dorado_call_extract.tsv \
				-a !{params.guppy_canon_threshold} -b !{params.guppy_mod_threshold}

		gzip -f  !{params.dsname}.!{call_tagname}.dorado_call_extract.tsv

		echo "### Dorado call extract DONE"
	'''
}


// convert per-read results to per-site/MethylKit results
process UNIFY {
	tag "${call_tagname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy"

	input:
	val tool_name	// e.g., Dorado
	val encode_name	// e.g., NANOME
	val call_tagname // e.g., all, H1, H2, chr1_22,etc.
	path per_read_fn // per-read file
	path reference_genome
	path ch_src
	path ch_utils

	output:
	path "Read_Level-${params.dsname}_${call_tagname}/*-perRead-score*.gz",	emit: read_unify, optional: true
	path "Site_Level-${params.dsname}_${call_tagname}/*-perSite-cov*.gz",	emit: site_unify, optional: true
	path "Site_Level-${params.dsname}_${call_tagname}/*.methylKit.CpG.txt.gz",	emit: methylkit_unify, optional: true

	when:
	params.dorado

	shell:
	cores = task.cpus * params.highProcTimes
	'''
		date; hostname; pwd

		echo !{per_read_fn}

		per_read_fn=!{per_read_fn}
		if [[ !{params.deduplicate} == true ]] ; then
			echo "### Deduplicate for read-level outputs"
			## sort order: Chr, Start, (End), ID, Strand
			zcat !{per_read_fn} |\
				sort -V -u -k2,2 -k3,3n -k1,1 -k4,4 |\
				gzip -f > !{per_read_fn}.sort.tsv.gz

			per_read_fn=!{per_read_fn}.sort.tsv.gz
		fi

		## Unify format output
		echo "### generate read/site level results"
		bash utils/unify_format_for_calls.sh \
			!{params.dsname}_!{call_tagname}  !{tool_name}  !{encode_name} \
			!{per_read_fn} \
			.  !{task.cpus}  12  !{params.sort ? true : false}  "!{params.chrSet1.replaceAll(',', ' ')}"

		## MethylKit format
		site_fn="Site_Level-!{params.dsname}_!{call_tagname}/!{params.dsname}_!{call_tagname}_!{tool_name}-perSite-cov1.sort.bed.gz"
		methylkit_fn="Site_Level-!{params.dsname}_!{call_tagname}/!{params.dsname}_!{call_tagname}_!{tool_name}.methylKit.CpG.txt.gz"

		python utils/nanome2methylkit.py \
			-i $site_fn \
			-o $methylkit_fn \
			--chrSet !{params.chrSet1.replaceAll(',', ' ')}

		echo "### Unify DONE"
	'''
}


process CLAIR3_dorado {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-phasing",
		mode: "copy", pattern: "${params.dsname}_clair3_out"

	publishDir "${params.outdir}/${params.dsname}-phasing",
		mode: "copy", pattern: "${params.dsname}_phased_intermediate"

	publishDir "${params.outdir}/${params.dsname}-phasing",
		mode: "copy", pattern: "${params.dsname}_phased_bam"

	publishDir "${params.outdir}/${params.dsname}-run-log",
		mode: "copy", pattern: "*.run.log"

	input:
	path dorado_call
	path reference_genome

	output:
	path "${params.dsname}_clair3_out",	emit:	clair3_out_ch, optional: true
	path "${params.dsname}_phased_intermediate",	emit:	phased_intermediate_out_ch, optional: true
	path "${params.dsname}_phased_bam",	emit:	phased_bam_out_ch, optional: true
	path "*.Clair3.run.log", optional:true,	emit: runlog
	path "${params.dsname}_phased_bam/${params.dsname}_HP1",	emit:	phased_bam_hp1_out_ch, optional: true
	path "${params.dsname}_phased_bam/${params.dsname}_HP2",	emit:	phased_bam_hp2_out_ch, optional: true

	"""
	run_clair3.sh --version

	bam_fn=${dorado_call}/*.bam
	MODEL_NAME="${params.CLAIR3_MODEL_NAME}"  ##"r941_prom_sup_g5014"

	mkdir -p ${params.dsname}_clair3_out
	run_clair3.sh \
	  --bam_fn=\${bam_fn} \
	  --ref_fn=${params.referenceGenome} \
	  --threads=${task.cpus} \
	  --platform="ont" \
	  --model_path="/opt/models/\${MODEL_NAME}" \
	  --enable_phasing --var_pct_phasing=${params.CLAIR3_var_pct_phasing} \
	  --output=${params.dsname}_clair3_out  ${params.ctg_name ? "--ctg_name=${params.ctg_name}": " "} \
	  &> ${params.dsname}.Clair3.run.log

	echo "### Clair3 for variant calling DONE"

	# haplotag
	# run whatshap haplotag
	mkdir -p ${params.dsname}_phased_intermediate
	phased_vcf_fn=${params.dsname}_clair3_out/phased_merge_output.vcf.gz
    tsvFile="${params.dsname}_phased_intermediate/${params.dsname}_whatshap_haplotag_read_list.tsv"
    haplotagBamFile="${params.dsname}_phased_intermediate/${params.dsname}_whatshap_haplotag_bam.bam"

    ## Phasing tag extraction for each chromosome
    ## older version lacks: --skip-missing-contigs  --output-threads ${task.cpus}
    whatshap  haplotag \
        --ignore-read-groups\
        --reference ${params.referenceGenome} \
        --output-haplotag-list \${tsvFile} \
        -o \${haplotagBamFile} \
        \${phased_vcf_fn}  \${bam_fn} \
        2>&1 | tee -a ${params.dsname}_phased_intermediate/${params.dsname}_whatshap_haplotag.run.log && \
    touch ${params.dsname}_phased_intermediate/${params.dsname}_whatshap_haplotag.done && \
    echo "### DONE for haplotag"

	# split bam
	mkdir -p ${params.dsname}_phased_bam
	outdir=${params.dsname}_phased_bam
	whatshap split \
		--output-h1 \${outdir}/${params.dsname}_split_HP1.bam \
		--output-h2 \${outdir}/${params.dsname}_split_HP2.bam \
		--output-untagged \${outdir}/${params.dsname}_split_untagged.bam  \
		\${bam_fn} \
		\${tsvFile} \
		2>&1 | tee \${outdir}/${params.dsname}_whatshap_split.run.log

	# Index bam files
	samtools index -@ ${task.cpus} \${outdir}/${params.dsname}_split_HP1.bam
	samtools index -@ ${task.cpus} \${outdir}/${params.dsname}_split_HP2.bam
	samtools index -@ ${task.cpus} \${outdir}/${params.dsname}_split_untagged.bam
	echo "### whatshap split DONE"

	mkdir -p \${outdir}/${params.dsname}_HP1 &&
		mv \${outdir}/${params.dsname}_split_HP1.bam* \${outdir}/${params.dsname}_HP1

	mkdir -p \${outdir}/${params.dsname}_HP2 &&
		mv \${outdir}/${params.dsname}_split_HP2.bam* \${outdir}/${params.dsname}_HP2
	"""
}
