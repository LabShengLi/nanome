/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : PHASING.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
process CLAIR3 {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-phasing",
		mode: "copy", pattern: "${params.dsname}_clair3_out"

	input:
	path merged_bam
	path reference_genome

	output:
	path "${params.dsname}_clair3_out",	emit:	clair3_out_ch, optional: true

	"""
	run_clair3.sh --version

	MODEL_NAME="r941_prom_sup_g5014"
	mkdir -p ${params.dsname}_clair3_out
	run_clair3.sh \
	  --bam_fn=${params.dsname}_bam_data/${params.dsname}_merge_all_bam.bam \
	  --ref_fn=${params.referenceGenome} \
	  --threads=${task.cpus} \
	  --platform="ont" \
	  --model_path="/opt/models/\${MODEL_NAME}" \
	  --output=${params.dsname}_clair3_out  ${params.ctg_name ? "--ctg_name=${params.ctg_name}": " "} \
	  &> ${params.dsname}.Clair3.run.log

	echo "### Clair3 for variant calling DONE"

	## haplotag
	whatshap --version

	## print header for file list tags
	head -n 1 \
		\$(find ${params.dsname}_clair3_out  -name '*_whatshap_haplotag_read_list_chr*.tsv' | head -n 1) \
    	> ${params.dsname}_clair3_out/${params.dsname}_haplotag_read_list_combine.tsv

	for chr in chr{1..22} chrX chrY; do
		if [[ ${params.ctg_name} != "null" &&  "${params.ctg_name}", != *"\$chr",* ]] ; then
			continue
		fi
		echo "### haplotag chr=\$chr"
		# run whatshap haplotag
		tsvFile="${params.dsname}_clair3_out/${params.dsname}_whatshap_haplotag_read_list_\$chr.tsv"
		haplotagBamFile="${params.dsname}_clair3_out/${params.dsname}_whatshap_haplotag_bam_\$chr.bam"
		phasingGZFile="${params.dsname}_clair3_out/tmp/phase_output/phase_vcf/phased_\$chr.vcf.gz"

		## Phasing tag extraction for each chromosome
		## older version lacks: --skip-missing-contigs  --output-threads ${task.cpus}
		whatshap  haplotag \
			--ignore-read-groups\
			--regions \${chr}\
			--reference ${params.referenceGenome}\
			--output-haplotag-list \${tsvFile} \
			-o \${haplotagBamFile} \
			\${phasingGZFile}  ${params.dsname}_bam_data/${params.dsname}_merge_all_bam.bam \
			&>> ${params.dsname}.Clair3.run.log

		if [[ ! -z "\${tsvFile}" ]]; then
			awk 'NR>1' \${tsvFile} \
				>> ${params.dsname}_clair3_out/${params.dsname}_haplotag_read_list_combine.tsv
		fi
		echo "### DONE for haplotag chr=\$chr"
	done

	# Extract h1 and h2 haplotype reads
	whatshap split \
		--output-h1 ${params.dsname}_clair3_out/${params.dsname}_split_h1.bam \
		--output-h2 ${params.dsname}_clair3_out/${params.dsname}_split_h2.bam \
		--output-untagged ${params.dsname}_clair3_out/${params.dsname}_split_untagged.bam  \
		${params.dsname}_bam_data/${params.dsname}_merge_all_bam.bam \
		${params.dsname}_clair3_out/${params.dsname}_haplotag_read_list_combine.tsv \
		&>> ${params.dsname}.Clair3.run.log
	echo "### Split by haplotag DONE"

	# Index bam files
	samtools index -@ ${task.cpus} ${params.dsname}_clair3_out/${params.dsname}_split_h1.bam
	samtools index -@ ${task.cpus} ${params.dsname}_clair3_out/${params.dsname}_split_h2.bam
	samtools index -@ ${task.cpus} ${params.dsname}_clair3_out/${params.dsname}_split_untagged.bam
	"""
}


process PHASING {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-phasing",
		mode: "copy"

	input:
	path mega_and_nanome_raw_list
	path clair3_out
	path ch_src
	path merged_bam
	path reference_genome

	output:
	path "hp_split_${params.dsname}*",	emit: hp_split_ch, 	optional: true
	path "${params.dsname}*mock_bam", 	emit: mock_bam_ch, 		optional: true

	"""
	echo "### hello phasing"

	## TODO: change hmc filename in Megalodon raw output
	toolList=(${params.hmc? "megalodon" : "megalodon"}  "nanome_${params.NANOME_MODEL}")
	encodeList=("megalodon" "nanome")
	numClassList=(${params.hmc? "3" : "2"}  2)

	for i in "\${!toolList[@]}"; do
		tool="\${toolList[i]}"
    	encode="\${encodeList[i]}"
    	numClass="\${numClassList[i]}"

		infn=\$(find . -name "${params.dsname}_\${tool}_per_read_combine*.gz")
		if [[ -z \${infn} ]] ; then
			continue
		fi

		echo "### tool=\${tool}, encode=\${encode}, infn=\${infn}"

		## Split methylation results by phasing tag for each chromosome
		for chr in chr{1..22} chrX chrY; do
			if [[ ${params.ctg_name} != "null" &&  "${params.ctg_name}", != *"\$chr",* ]] ; then
				continue
			fi

			## Step1: HP split meth data
			echo "### HP split for chr=\${chr}"
			PYTHONPATH=src python src/nanome/other/phasing/hp_split.py \
				--dsname ${params.dsname}\
				--tool \${tool}\
				--encode \${encode}\
				--num-class \${numClass}\
				-i \${infn}\
				--haplotype-list ${params.dsname}_clair3_out/${params.dsname}_whatshap_haplotag_read_list_\${chr}.tsv\
				--region \${chr}\
				-o .  --save-unified-read  &>> ${params.dsname}.Phasing.run.log

			## Start generate mocked BAM files
			## Step2: methcall2bed
			## hp_split_NA12878_CHR22_200_megalodon
			outdir=${params.dsname}_\${tool}_methcall2bed
			mkdir -p \${outdir}
			find hp_split_${params.dsname}_\${tool} -name "${params.dsname}*_perReadScore_\${chr}_H*.tsv.gz" -print0 |
				while IFS= read -r -d '' infn2; do
					basefn=\$(basename \$infn2)
					outfn=\${outdir}/\${basefn/.tsv.gz/_methcall2bed.bed.gz}
					PYTHONPATH=src  python src/nanome/other/phasing/methcall2bed.py \
						-i \${infn2} \
						-o \${outfn} \
						--verbose  &>> ${params.dsname}.Phasing.run.log

					zcat \${outfn} | sort -V -k1,1 -k2,2n -k3,3n |
						bgzip -f >\${outfn/.bed.gz/.sort.bed.gz} &&
						tabix -p bed \${outfn/.bed.gz/.sort.bed.gz}
					rm -f \${outfn}
					touch \${outfn/.bed.gz/.sort.bed.gz}.DONE
				done

			## Step3: bam2bis
			outdir2=${params.dsname}_\${tool}_mock_bam
			mkdir -p \${outdir2}
			bamFile=\$(find ${merged_bam}/ -name "*.bam")
			for hapType in H1 H2 H1_5hmc H2_5hmc; do
				methCallFile=\$(find \${outdir} -name "${params.dsname}_\${tool,,}_perReadScore_\${chr}_\${hapType}_methcall2bed.sort.bed.gz")
				if [ ! -e "\${methCallFile}" ] ; then
					continue
				fi
				PYTHONPATH=src  python  src/nanome/other/phasing/nanomethphase.py bam2bis \
					--bam \${bamFile} \
					--reference ${params.referenceGenome} \
					--methylcallfile \${methCallFile} \
					--output \${outdir2}/${params.dsname}_\${tool}_\${chr}_\${hapType} \
					-t ${task.cpus} --window \${chr} --overwrite  &>> ${params.dsname}.Phasing.run.log

				infn3=\$(find \${outdir2} -name "${params.dsname}_\${tool}_\${chr}_\${hapType}*.bam")
				if [ ! -e "\${infn}" ] ; then
					continue
				fi

				samtools sort -@ ${task.cpus} \$infn3 -o \${infn3/.bam/.sort.bam} &&
					samtools index -@ ${task.cpus} \${infn3/.bam/.sort.bam} &&
					rm -f \${infn3} &&
					touch \${infn3/.bam/.sort.bam}.DONE
			done
		done
	done
	"""
}
