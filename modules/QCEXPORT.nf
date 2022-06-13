// Collect and output QC results for basecall, and report ONT coverage
process QCEXPORT {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-basecallings",
		mode: "copy", enabled: params.outputQC, overwrite: true

	input:
	path basecall_list
	path alignment_list
	path reference_genome

	output:
	path "${params.dsname}_basecall_report.html",	optional: true, emit: qc_html
	path "${params.dsname}_QCReport",				emit: qc_report
	path "${params.dsname}_bam_data",				optional: true,	 emit: bam_data

	shell:
	cores = task.cpus * params.highProcTimes
	samtools_cores = task.cpus * params.mediumProcTimes
	'''
	## Combine all sequencing summary files
	touch !{params.dsname}_combine_sequencing_summary.txt.gz
	firstFile=true
	find *.basecall/ -name '*-sequencing_summary.txt' -type f -print0 |\
		while read -d $'\0' file ; do
			if $firstFile ; then
				awk 'NR>=1' $file | \
					gzip -f >> !{params.dsname}_combine_sequencing_summary.txt.gz
				firstFile=false
			else
				awk 'NR>1' $file | \
					gzip -f >> !{params.dsname}_combine_sequencing_summary.txt.gz
			fi
		done

	mkdir -p !{params.dsname}_QCReport
	if [[ !{params.skipQC} == false ]]; then
		## Perform QC report by NanoComp
		NanoComp --summary !{params.dsname}_combine_sequencing_summary.txt.gz  \
			--names !{params.dsname} --outdir !{params.dsname}_QCReport -t !{cores} \
			--raw  -f pdf -p !{params.dsname}_   &>> !{params.dsname}.QCReport.run.log
	fi

	if [[ !{params.outputBam} == true  || !{params.outputONTCoverage} == true || !{params.phasing} == true ]]; then
		## Combine all bam files
		samtools merge -@ !{samtools_cores}  !{params.dsname}_merge_all_bam.bam  *.alignment/*_bam.bam  &&\
			samtools index -@ !{samtools_cores}   !{params.dsname}_merge_all_bam.bam
		echo "### Samtools merge done"
	fi

	if [[ !{params.outputONTCoverage} == true ]]; then
		## calculates the sequence coverage at each position
		## reporting genome coverage for all positions in BEDGRAPH format.
		bedtools genomecov -ibam !{params.dsname}_merge_all_bam.bam -bg -strand + |
			awk '$4 = $4 FS "+"' |
			gzip -f > !{params.dsname}.coverage.positivestrand.bed.gz

		bedtools genomecov -ibam !{params.dsname}_merge_all_bam.bam -bg -strand - |
			awk '$4 = $4 FS "-"' |
			gzip -f > !{params.dsname}.coverage.negativestrand.bed.gz

		cat !{params.dsname}.coverage.positivestrand.bed.gz > \
			!{params.dsname}_ONT_coverage_combine.bed.gz
		cat !{params.dsname}.coverage.negativestrand.bed.gz >> \
			!{params.dsname}_ONT_coverage_combine.bed.gz

		mv !{params.dsname}_ONT_coverage_combine.bed.gz !{params.dsname}_QCReport/
	fi

	[ -f !{params.dsname}_combine_sequencing_summary.txt.gz ] && \
		mv -f !{params.dsname}_combine_sequencing_summary.txt.gz !{params.dsname}_QCReport/
	[ -f !{params.dsname}_QCReport/!{params.dsname}_NanoComp-report.html ] && \
		mv -f !{params.dsname}_QCReport/!{params.dsname}_NanoComp-report.html \
		 	!{params.dsname}_basecall_report.html

	## Clean
	if [[ !{params.cleanStep} == "true" ]]; then
		rm -f !{params.dsname}.coverage.positivestrand.bed.gz \
		 	!{params.dsname}.coverage.negativestrand.bed.gz
		rm -f merge_all_fq.fq.gz
		if [[ !{params.outputBam} == false && !{params.phasing} == false ]]; then
			rm -f !{params.dsname}_merge_all_bam.bam*
		else
			mkdir -p !{params.dsname}_bam_data
			mv  !{params.dsname}_merge_all_bam.bam*  !{params.dsname}_bam_data/
		fi
	fi
    echo "### ONT coverage done!"
    echo "### QCReport all DONE"
	'''
}