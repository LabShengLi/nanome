/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : CONSENSUS.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
process CONSENSUS {
	tag "${params.dsname}"

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
	path read_fileList
	path ch_src
	path ch_utils

	output:
	path "Read_Level-${params.dsname}/${params.dsname}_*-perRead-score*.gz",	emit: read_unify, optional: true
	path "Site_Level-${params.dsname}/*-perSite-cov*.gz",	emit: site_unify, optional: true
	path "${params.dsname}_nanome_${params.NANOME_MODEL}_per_read_combine.*.gz", emit: nanome_combine_out, optional: true

	when:
	params.runNANOME && (params.runNanopolish || params.runDeepSignal2 || params.runMegalodon)

	"""
	if [[ ${params.NANOME_MODEL} == "nanome_cs" ]] ; then
		echo "### nanome_cs"
		## check if consensus method input exists
		canRun="false"

		if test -n "\$(find . -maxdepth 1 -name '*.deepsignal1_batch_features.tsv.gz' -print -quit)"
		then
			echo "### found deepsignal1_batch_features"
			cat *.deepsignal1_batch_features.tsv.gz > \
				${params.dsname}_deepsignal_feature_combine.tsv.gz
		else
			echo "### not found deepsignal_batch_features"
		fi

		MegalodonReadReport=\$(find . -maxdepth 1 -name '*Megalodon-perRead-score.tsv.gz')
		if [[ -z \$MegalodonReadReport ]] ; then
			echo "### Not found Megalodon read-level outputs"
			MegalodonOptions=" "
		else
			MegalodonOptions="--megalodon \$MegalodonReadReport"
			canRun="true"
		fi

		NanopolishReadReport=\$(find . -maxdepth 1 -name '*Nanopolish-perRead-score.tsv.gz')
		if [[ -z \$NanopolishReadReport ]] ; then
			echo "### Not found Nanopolish read-level outputs"
			NanopolishOptions=" "
		else
			NanopolishOptions="--nanopolish \$NanopolishReadReport"
			canRun="true"
		fi

		DeepSignalReadReport=\$(find . -maxdepth 1 -name '*DeepSignal*-perRead-score.tsv.gz' | head -n 1)
		if [[ -z \$DeepSignalReadReport ]] ; then
			echo "### Not found DeepSignal read-level outputs"
			DeepSignalOptions=" "
		else
			DeepSignalOptions="--deepsignal \$DeepSignalReadReport"
			canRun="true"
		fi

		FeatureFile=\$(find . -maxdepth 1 -name '*_deepsignal*_feature_combine.tsv.gz' | head -n 1)
		if [[ -z \$FeatureFile ]] ; then
			echo "### Not found Feature file"
			FeatureOptions=" "
		else
			FeatureOptions="--feature \$FeatureFile"
		fi

		if [[ \$canRun == false ]] ; then
			echo "No input for NANOME, exit"
			exit 0
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
					&>> ${params.dsname}.Consensus.run.log
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
				&>> ${params.dsname}.Consensus.run.log
		fi
	elif [[ ${params.NANOME_MODEL} == "NANOME3T" ]] ; then
		## NANOME XGBoost method, will be deprecated
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
						&>> ${params.dsname}.Consensus.run.log
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
					&>> ${params.dsname}.Consensus.run.log
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

	echo "### NANOME Consensus model DONE"
	"""
}
