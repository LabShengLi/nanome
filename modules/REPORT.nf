/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : REPORT.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// Not cache due to the script contains run information, each time of resume run will need updated
process REPORT {
	tag "${params.dsname}"

	publishDir "${params.outdir}",
		mode: "copy", pattern: "README_${params.dsname}.txt"

	publishDir "${params.outdir}",
		mode: "copy", pattern: "${params.dsname}_nanome_report.html"

	publishDir "${params.outdir}/MultiQC",
		mode: "copy", pattern: "multiqc_report.html"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy", pattern: "GenomeBrowser-${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings",
		mode: "copy",
		pattern: "Site_Level-${params.dsname}/*-perSite-cov*.gz"

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Raw_Results-${params.dsname}",
		mode: "copy",
		pattern: "${params.dsname}_nanome_*_per_read_combine.*.gz",
		enabled: params.outputRaw

	input:
	path site_fileList
	path read_fileList
	path tools_version_tsv
	path qc_report
	path reference_genome
	path ch_src
	path ch_utils

	output:
	path "${params.dsname}_nanome_report.html",	emit:	report_out_ch
	path "README_${params.dsname}.txt",	emit: 	readme_out_ch
	path "multiqc_report.html",	emit: 	lbt_report_ch
	path "GenomeBrowser-${params.dsname}", emit:  genome_browser_ch, optional: true
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify, optional: true
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify, optional: true
	path "${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.*.gz", emit: nanome_combine_out, optional: true

	when:
	params.runNANOME

	"""
	if [[ ${params.NANOME_MODEL} == "NANOME3T_W" ]] ; then
		echo "### NANOME3T_W"

		cat *.deepsignal1_batch_features.tsv.gz > \
			${params.dsname}_deepsignal1_feature_combine.tsv.gz

		MegalodonReadReport=\$(find . -maxdepth 1 -name '*Megalodon-perRead-score.tsv.gz')
		if [[ -z \$MegalodonReadReport ]] ; then
			echo "### Not found Megalodon read-level outputs"
			MegalodonOptions=" "
		else
			MegalodonOptions="--megalodon_file \$MegalodonReadReport"
		fi

		NanopolishReadReport=\$(find . -maxdepth 1 -name '*Nanopolish-perRead-score.tsv.gz')
		if [[ -z \$NanopolishReadReport ]] ; then
			echo "### Not found Nanopolish read-level outputs"
			NanopolishOptions=" "
		else
			NanopolishOptions="--nanopolish_file \$NanopolishReadReport"
		fi

		DeepSignalReadReport=\$(find . -maxdepth 1 -name '*DeepSignal-perRead-score.tsv.gz')
		if [[ -z \$DeepSignalReadReport ]] ; then
			echo "### Not found DeepSignal read-level outputs"
			DeepSignalOptions=" "
		else
			DeepSignalOptions="--deepsignal_file \$DeepSignalReadReport"
		fi

		FeatureFile=\$(find . -maxdepth 1 -name '*_deepsignal1_feature_combine.tsv.gz')
		if [[ -z \$FeatureFile ]] ; then
			echo "### Not found Feature file"
			FeatureOptions=" "
		else
			FeatureOptions="--feature_file \$FeatureFile"
		fi

		importmodelfile=${params.NANOME_MODEL_FILE}

		if [[ ${params.consensus_by_chr} == true ]] ; then
			mkdir -p consensus_by_chr
			for chr in ${params.chrSet1.replaceAll(',', ' ')} ; do
				echo "### consensus for chr=\${chr}"
				PYTHONPATH=src python src/nanome/xgboost/weight/model_predict.py \
					--dsname ${params.dsname} \
					--importmodelfile \${importmodelfile} \
					--base_model_type ${params.NANOME_MODEL_BASE} \
					--specific_model_type ${params.NANOME_SPECIFIC_MODEL_TYPE} \
					--one_hot_sequence true  --chrs \${chr} \
					\${MegalodonOptions}  \${NanopolishOptions} \${DeepSignalOptions} \
					\${FeatureOptions} --outdir consensus_by_chr\
					&>> ${params.dsname}.Report.run.log
			done
			## combine all chrs
			touch ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz
			zcat \$(ls consensus_by_chr/*.gz| head -n 1) | head -n 1 | \
				awk -v FS='\t' -v OFS='\t' '{print \$4,\$1,\$2,\$3,\$5,\$6,\$7,\$8,\$9,\$10}' |\
				gzip -f >> \
				${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz

			find consensus_by_chr -name '*.tsv.gz' -print0 |
				sort -V -z |
				while IFS= read -r -d '' infn; do
					zcat \${infn}  | \
					awk -v FS='\t' -v OFS='\t' 'NR>1 {print \$4,\$1,\$2,\$3,\$5,\$6,\$7,\$8,\$9,\$10}' | \
					gzip -f >> \
						${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz
				done
			echo "### consensus combine chr DONE"
		else
			echo PYTHONPATH=src python src/nanome/xgboost/weight/model_predict.py \
				--dsname ${params.dsname} \
				--importmodelfile \${importmodelfile} \
				--base_model_type ${params.NANOME_MODEL_BASE} \
				--specific_model_type ${params.NANOME_SPECIFIC_MODEL_TYPE} \
				--one_hot_sequence true  \
				\${MegalodonOptions}  \${NanopolishOptions} \${DeepSignalOptions} \
				\${FeatureOptions} --outdir .
			echo "### Not support yet"
			exit 1
		fi
		echo "### TODO DONE"
	elif [[ ${params.NANOME_MODEL} == "nanome_cs" ]] ; then
		echo "### nanome_cs"

		if test -n "\$(find . -maxdepth 1 -name '*.deepsignal1_batch_features.tsv.gz' -print -quit)"
		then
			echo "### found deepsignal1_batch_features"
			cat *.deepsignal1_batch_features.tsv.gz > \
				${params.dsname}_deepsignal1_feature_combine.tsv.gz
		else
			echo "### not found deepsignal1_batch_features"
		fi

		MegalodonReadReport=\$(find . -maxdepth 1 -name '*Megalodon-perRead-score.tsv.gz')
		if [[ -z \$MegalodonReadReport ]] ; then
			echo "### Not found Megalodon read-level outputs"
			MegalodonOptions=" "
		else
			MegalodonOptions="--megalodon \$MegalodonReadReport"
		fi

		NanopolishReadReport=\$(find . -maxdepth 1 -name '*Nanopolish-perRead-score.tsv.gz')
		if [[ -z \$NanopolishReadReport ]] ; then
			echo "### Not found Nanopolish read-level outputs"
			NanopolishOptions=" "
		else
			NanopolishOptions="--nanopolish \$NanopolishReadReport"
		fi

		DeepSignalReadReport=\$(find . -maxdepth 1 -name '*DeepSignal-perRead-score.tsv.gz')
		if [[ -z \$DeepSignalReadReport ]] ; then
			echo "### Not found DeepSignal read-level outputs"
			DeepSignalOptions=" "
		else
			DeepSignalOptions="--deepsignal \$DeepSignalReadReport"
		fi

		FeatureFile=\$(find . -maxdepth 1 -name '*_deepsignal1_feature_combine.tsv.gz')
		if [[ -z \$FeatureFile ]] ; then
			echo "### Not found Feature file"
			FeatureOptions=" "
		else
			FeatureOptions="--feature \$FeatureFile"
		fi

		if [[ ${params.consensus_by_chr} == true ]] ; then
			mkdir -p consensus_by_chr
			for chr in ${params.chrSet1.replaceAll(',', ' ')} ; do
				echo "### consensus for chr=\${chr}"
				PYTHONPATH=src python src/nanome/xgboost/cs_predict.py \
					--dsname ${params.dsname} \
					-m ${params.CS_MODEL_FILE} --model_specific ${params.CS_MODEL_SPEC}  \
					--chrs \${chr} \
					\${MegalodonOptions}  \${NanopolishOptions} \${DeepSignalOptions} \
					\${FeatureOptions} \
					-o  consensus_by_chr/${params.dsname}_nanome_${params.CS_MODEL_SPEC}_per_read_\${chr}.tsv.gz   \
					&>> ${params.dsname}.Report.run.log
			done
			## combine all chrs
			touch ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz
			zcat \$(ls consensus_by_chr/*.gz| head -n 1) | head -n 1 | \
				gzip -f >> \
				${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz

			find consensus_by_chr -name '*.tsv.gz' -print0 |
				sort -V -z |
				while IFS= read -r -d '' infn; do
					zcat \${infn}  | awk 'NR>1' |\
					gzip -f >> \
						${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz
				done
			echo "### consensus combine chr DONE"
		else
			PYTHONPATH=src python src/nanome/xgboost/cs_predict.py \
				--dsname ${params.dsname} \
				-m ${params.CS_MODEL_FILE} --model_specific ${params.CS_MODEL_SPEC}  \
				--chrs \${chr} \
				\${MegalodonOptions}  \${NanopolishOptions} \${DeepSignalOptions} \
				\${FeatureOptions} \
				-o  ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz   \
				&>> ${params.dsname}.Report.run.log
		fi
	elif [[ ${params.NANOME_MODEL} == "NANOME3T" ]] ; then
		## NANOME XGBoost method
		modelContentTSVFileName=${params.dsname}_nanome_${params.NANOME_MODEL}_model_content.tsv
		> \$modelContentTSVFileName
		passModelTsv=false
		if [[ "${params.NANOME_CONSENSUS_TOOLS}" == *"Nanopolish"* ]]; then
			NanopolishReadReport=\$(find . -maxdepth 1 -name '*Nanopolish-perRead-score.tsv.gz')
			if [[ -z \$NanopolishReadReport ]] ; then
				echo "### Not found Nanopolish read-level outputs"
				NanopolishReadReport="None"
			else
				passModelTsv=true
			fi
			printf '%s\t%s\n' nanopolish \${NanopolishReadReport} >> \$modelContentTSVFileName
		fi

		if [[ "${params.NANOME_CONSENSUS_TOOLS}" == *"Megalodon"* ]]; then
			MegalodonReadReport=\$(find . -maxdepth 1 -name '*Megalodon-perRead-score.tsv.gz')
			if [[ -z \$MegalodonReadReport ]] ; then
				echo "### Not found Megalodon read-level outputs"
				MegalodonReadReport="None"
			else
				passModelTsv=true
			fi
			printf '%s\t%s\n' megalodon \${MegalodonReadReport} >> \$modelContentTSVFileName
		fi

		if [[ "${params.NANOME_CONSENSUS_TOOLS}" == *"DeepSignal"* ]]; then
			DeepSignalReadReport=\$(find . -maxdepth 1 -name '*DeepSignal-perRead-score.tsv.gz')
			if [[ -z \$DeepSignalReadReport ]] ; then
				echo "### Not found DeepSignal read-level outputs"
				DeepSignalReadReport="None"
			else
				passModelTsv=true
			fi
			printf '%s\t%s\n' deepsignal \${DeepSignalReadReport} >> \$modelContentTSVFileName
		fi

		if [[ "\$passModelTsv" == true ]] ; then
			## NANOME XGBoost model results, if there are model results exists
			echo "### NANOME XGBoost predictions"

			## 0.23.2 version work both for NANOME>=0.23.2 and METEORE<=0.23.2
			## pip install -U scikit-learn==0.23.2

			pip show scikit-learn

			if [[ ${params.consensus_by_chr} == true ]] ; then
				mkdir -p consensus_by_chr
				for chr in ${params.chrSet1.replaceAll(',', ' ')} ; do
					echo "### consensus for chr=\${chr}"
					PYTHONPATH=src python src/nanome/xgboost/xgboost_predict.py \
						--tsv-input\
						--dsname ${params.dsname} -i \${modelContentTSVFileName}\
						-m ${params.NANOME_MODEL}  \
						-o consensus_by_chr/${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_\${chr}.tsv.gz \
						--chrs \${chr} \
						&>> ${params.dsname}.Report.run.log
				done
				## combine all chrs
				touch ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz
				zcat \$(ls consensus_by_chr/*.gz| head -n 1) | head -n 1 | gzip -f >> \
					${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz

				find consensus_by_chr -name '*.tsv.gz' -print0 |
					sort -V -z |
   					while IFS= read -r -d '' infn; do
						zcat \${infn}  | awk 'NR>1' | gzip -f >> \
							${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz
					done
				echo "### consensus combine chr DONE"
			else
				PYTHONPATH=src python src/nanome/xgboost/xgboost_predict.py \
					--tsv-input\
					--dsname ${params.dsname} -i \${modelContentTSVFileName}\
					-m ${params.NANOME_MODEL}  \
					-o ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz \
					&>> ${params.dsname}.Report.run.log
			fi
		fi
	fi

	if [[ -f ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz ]] ; then
		if [[ ${params.deduplicate} == true ]] ; then
			echo "### Deduplicate for read-level outputs"
			## sort order: Chr, Start, (End), ID, Strand
			zcat ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz |\
				sort -V -u -k2,2 -k3,3n -k1,1 -k4,4 |\
				gzip -f > ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.sort.tsv.gz
			rm ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz &&\
				mv ${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.sort.tsv.gz\
					${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz
		fi

		## Unify format output
		echo "### NANOME read/site level results"
		bash utils/unify_format_for_calls.sh \
			${params.dsname}  NANOME NANOME\
			${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.tsv.gz \
			.  $task.cpus  12  ${params.sort ? true : false}  "${params.chrSet1.replaceAll(',', ' ')}"
		ln -s Site_Level-${params.dsname}/${params.dsname}_NANOME-perSite-cov1.sort.bed.gz\
			${params.dsname}_NANOME-perSite-cov1.sort.bed.gz
	fi

	## Generate NF pipeline running information tsv
	> running_information.tsv
	printf '%s\t%s\n' Title Information >> running_information.tsv
	printf '%s\t%s\n' dsname ${params.dsname} >> running_information.tsv
	printf '%s\t%s\n' projectDir ${workflow.projectDir} >> running_information.tsv
	printf '%s\t%s\n' workDir ${workflow.workDir} >> running_information.tsv
	printf '%s\t%s\n' commandLine "${workflow.commandLine}" >> running_information.tsv
	printf '%s\t%s\n' runName ${workflow.runName} >> running_information.tsv
	printf '%s\t%s\n' start ${workflow.start} >> running_information.tsv
	printf '%s\t%s\n' input "${params.input}" >> running_information.tsv
	printf '%s\t%s\n' outputs ${params.outdir} >> running_information.tsv

	## Note that the reason of report process can not be cached, is due to
	## Above script codes will be changed each time, so report can not apply old cached script

	## Get basecalling results from NanoComp
	basecallOutputFile=\$(find ${params.dsname}_QCReport/ -name "*NanoStats.txt" -type f)

	if [[ -z "\${basecallOutputFile}" ]] ; then
		basecallOutputFile=None
	fi

	## Generate report dir and html utilities
	if [ -d /opt/nanome ]; then
		nanome_dir=/opt/nanome
	else
		nanome_dir="."
	fi
	mkdir -p ${params.dsname}_NANOME_report
	cp \${nanome_dir}/src/nanome/nanocompare/report/style.css ${params.dsname}_NANOME_report/
	cp -rf \${nanome_dir}/src/nanome/nanocompare/report/icons ${params.dsname}_NANOME_report/
	cp -rf \${nanome_dir}/src/nanome/nanocompare/report/js ${params.dsname}_NANOME_report/

	## Generate html NANOME report
	PYTHONPATH=src python src/nanome/nanocompare/report/gen_html_report.py\
		${params.dsname} \
		running_information.tsv \
		\${basecallOutputFile} \
		. \
		${params.dsname}_NANOME_report \
		./src/nanome/nanocompare/report\
		${tools_version_tsv}  &>> ${params.dsname}.Report.run.log

	## Combine a single html report
	## No tty usage, ref: https://github.com/remy/inliner/issues/151
	script -qec "inliner ${params.dsname}_NANOME_report/nanome_report.html" /dev/null \
	  	> ${params.dsname}_nanome_report.html

	## Used for lifebit rendering feature
	cp ${params.dsname}_nanome_report.html   multiqc_report.html

	## Generate readme.txt
	PYTHONPATH=src PYTHONIOENCODING=UTF-8 python src/nanome/nanocompare/report/gen_txt_readme.py\
		src/nanome/nanocompare/report/readme.txt.template ${params.dsname} ${params.outdir}\
		${workflow.projectDir} ${workflow.workDir} "${workflow.commandLine}"\
		${workflow.runName} "${workflow.start}"\
		> README_${params.dsname}.txt   2>> ${params.dsname}.Report.run.log

	## Output BigWig format for IGV
	if [[ ${params.outputGenomeBrowser} == true ]] ; then
		if command -v bedGraphToBigWig  &> /dev/null ; then
			mkdir -p GenomeBrowser-${params.dsname}
			find . -maxdepth 1 -name '*-perSite-cov1.sort.bed.gz' -print0 | \
				while IFS= read -r -d '' infn ; do
					echo "### processing infn=\$infn"
					## methfreq bw generation
					zcat \${infn} | awk '{printf "%s\\t%d\\t%d\\t%2.5f\\n" , \$1,\$2,\$3,\$7}' > \
						GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.bedgraph}
					LC_COLLATE=C sort -u -k1,1 -k2,2n \
						GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.bedgraph} > \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.sorted.bedgraph}

					## Check if bedgraph is empty, issue ref: https://biostar.usegalaxy.org/p/6794/
					if [[ ! -s GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.sorted.bedgraph} ]] ; then
						continue
					fi

					bedGraphToBigWig GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.sorted.bedgraph} \
						reference_genome/chrom.sizes   GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.bw}
					rm -f GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.bedgraph}  \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_methfreq.sorted.bedgraph}

					## coverage bw generation
					zcat \${infn} | \
						awk '{printf "%s\\t%d\\t%d\\t%d\\n" , \$1,\$2,\$3,\$8}' > \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.bedgraph}
					LC_COLLATE=C sort -u -k1,1 -k2,2n \
						GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.bedgraph} > \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.sorted.bedgraph}

					## Check if bedgraph is empty
					if [[ ! -s GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.sorted.bedgraph} ]] ; then
						continue
					fi
					bedGraphToBigWig GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.sorted.bedgraph} \
						reference_genome/chrom.sizes   GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.bw}
					rm -f GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.bedgraph}  \
							GenomeBrowser-${params.dsname}/\${infn/-perSite-cov1.sort.bed.gz/_coverage.sorted.bedgraph}
				done
		else
			echo "### ERROR: No bedGraphToBigWig in PATH, please install it"
		fi
	fi
	echo "### report html DONE"
	"""
}
