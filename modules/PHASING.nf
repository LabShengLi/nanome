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
 @Organization : Sheng Li Lab
----------------------------------------------------------------------------------------
*/
process CLAIR3 {
	tag "${params.dsname}"

	publishDir "${params.outdir}/${params.dsname}-phasing",
		mode: "copy", pattern: "${params.dsname}_clair3_out"

	publishDir "${params.outdir}/${params.dsname}-run-log",
		mode: "copy", pattern: "*.Clair3.run.log"

	input:
	path merged_bam
	path reference_genome

	output:
	path "${params.dsname}_clair3_out",	emit:	clair3_out_ch, optional: true
	path "*.Clair3.run.log", optional:true,	emit: runlog

	"""
	run_clair3.sh --version

	MODEL_NAME="${params.CLAIR3_MODEL_NAME}"  ##"r941_prom_sup_g5014"
	mkdir -p ${params.dsname}_clair3_out
	run_clair3.sh \
	  --bam_fn=${params.dsname}_bam_data/${params.dsname}_merge_all_bam.bam \
	  --ref_fn=${params.referenceGenome} \
	  --threads=${task.cpus} \
	  --platform="ont" \
	  --model_path="/opt/models/\${MODEL_NAME}" \
	  --enable_phasing --var_pct_phasing=${params.CLAIR3_var_pct_phasing} \
	  --output=${params.dsname}_clair3_out  ${params.ctg_name ? "--ctg_name=${params.ctg_name}": " "} \
	  &> ${params.dsname}.Clair3.run.log

	echo "### Clair3 for variant calling DONE"

    ## combine phasing vcf files
    phase_dir=${params.dsname}_clair3_out/tmp/phase_output/phase_vcf
    > \${phase_dir}/${params.dsname}_phasing_vcf.vcf

    hfn=\$(find \${phase_dir} -name '*.vcf.gz' | head -n 1)
    zcat \${hfn} | awk '\$1 ~ /^#/' >> \
        \${phase_dir}/${params.dsname}_phasing_vcf.vcf

    find \${phase_dir} -name '*.vcf.gz'  -print0 |
        sort -V -z |
	    while IFS= read -r -d '' infn2; do
	        zcat \${infn2} | awk '\$1 !~ /^#/' >> \
	            \${phase_dir}/${params.dsname}_phasing_vcf.vcf
	    done

	> \${phase_dir}/${params.dsname}_phasing_vcf_QUAL_${params.CLAIR3_phasing_qual}.vcf
	awk '\$1 ~ /^#/' \${phase_dir}/${params.dsname}_phasing_vcf.vcf >> \
        \${phase_dir}/${params.dsname}_phasing_vcf_QUAL_${params.CLAIR3_phasing_qual}.vcf
    awk "\\\$6 >= ${params.CLAIR3_phasing_qual}" \${phase_dir}/${params.dsname}_phasing_vcf.vcf >> \
        \${phase_dir}/${params.dsname}_phasing_vcf_QUAL_${params.CLAIR3_phasing_qual}.vcf

    bgzip -c \${phase_dir}/${params.dsname}_phasing_vcf.vcf > \
        \${phase_dir}/${params.dsname}_phasing_vcf.vcf.gz
    tabix -p vcf \${phase_dir}/${params.dsname}_phasing_vcf.vcf.gz

    bgzip -c \${phase_dir}/${params.dsname}_phasing_vcf_QUAL_${params.CLAIR3_phasing_qual}.vcf > \
        \${phase_dir}/${params.dsname}_phasing_vcf_QUAL_${params.CLAIR3_phasing_qual}.vcf.gz
    tabix -p vcf \${phase_dir}/${params.dsname}_phasing_vcf_QUAL_${params.CLAIR3_phasing_qual}.vcf.gz

	echo "### combine phasing vcf DONE"

	## haplotag
	whatshap --version

	# TODO:
	## print header for file list tags,  there will be no this file at firsts
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

		# index phased vcf.gz
	 	tabix -p vcf \${phasingGZFile}

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

	publishDir "${params.outdir}/${params.dsname}-run-log",
		mode: "copy", pattern: "*.Phasing.run.log"

	input:
	path meth_for_phasing_inputs
	path clair3_out
	path ch_src
	path merged_bam
	path reference_genome

	output:
	path "hp_split_${params.dsname}*",	emit: hp_split_ch, 	optional: true
	path "${params.dsname}*mock_bam", 	emit: mock_bam_ch, 		optional: true
    path "${params.dsname}*methcall2bed", 	emit: methcall2bed_ch, 		optional: true
	path "${params.dsname}*meth_phasing", 	emit: meth_phasing_ch, 		optional: true
	path "*.Phasing.run.log", optional:true,	emit: runlog

	"""
	echo "### start phasing"

	# manner1
	if [[ ${params.phase_manner1} == true ]] ; then
    ## deal with meth2bed+nanomethphase_phase , phased meth looks good
    ## phaseToolList=("nanopolish" "megalodon" "nanome" "deepsignal" "guppy")
    phaseToolList=(${params.phasing_tools.replaceAll(',',' ')})
    for i in "\${!phaseToolList[@]}"; do
		tool="\${phaseToolList[i]}"

		if [[ \${tool} == "nanopolish" ]] ; then
            infn=\$(find . -maxdepth 1 \
                -name "${params.dsname}_\${tool}*_per_read_combine*.gz")
        else
            ## find file: NA19240_chr11_chr15_n49_Megalodon-perRead-score.tsv.gz
            infn=\$(find . -maxdepth 1 \
                -iname "${params.dsname}_\${tool}*perRead-score.*.gz")
        fi
		if [[ -z \${infn} ]] ; then
			continue
		fi

		echo "deal with \${infn}"

        outdir="${params.dsname}_\${tool}_meth_phasing"
        mkdir -p \${outdir}
		if [[ \${tool} == "nanopolish" ]] ; then
			## nanopolish way, keep same with NanomethPhase
            PYTHONPATH=src python src/nanome/other/phasing/nanomethphase.py \
                methyl_call_processor -mc \${infn} -t ${task.cpus} |
                sort -k1,1 -k2,2n -k3,3n |
                bgzip -f > \${outdir}/${params.dsname}_\${tool}_MethylationCall.bed.gz &&
                tabix -p bed \${outdir}/${params.dsname}_\${tool}_MethylationCall.bed.gz
		else
            PYTHONPATH=src  python src/nanome/other/phasing/methcall2bed.py \
                -i \${infn} \
                -o \${outdir}/${params.dsname}_\${tool}_MethylationCall1.bed.gz \
                --score-cutoff ${params.PHASE_meth_score_cutoff} \
                --verbose  &>> ${params.dsname}.Phasing.run.log

            zcat \${outdir}/${params.dsname}_\${tool}_MethylationCall1.bed.gz | \
                sort -V -k1,1 -k2,2n -k3,3n |
                bgzip -f >\${outdir}/${params.dsname}_\${tool}_MethylationCall.bed.gz &&
                tabix -p bed \${outdir}/${params.dsname}_\${tool}_MethylationCall.bed.gz  &&
                rm -f \${outdir}/${params.dsname}_\${tool}_MethylationCall1.bed.gz
		fi

		## phase meth, results looks good
		vcffn=${params.dsname}_clair3_out/tmp/phase_output/phase_vcf/${params.dsname}_phasing_vcf_QUAL_${params.CLAIR3_phasing_qual}.vcf.gz
		ref=${params.referenceGenome}
		bamfn=${params.dsname}_bam_data/${params.dsname}_merge_all_bam.bam

		PYTHONPATH=src python src/nanome/other/phasing/nanomethphase.py \
		    phase -mc \${outdir}/${params.dsname}_\${tool}_MethylationCall.bed.gz \
            -o \${outdir}/${params.dsname}_\${tool} \
            -of bam,methylcall,bam2bis \
            -b \${bamfn} -r \${ref} -v \${vcffn} -t ${task.cpus}  --overwrite

        ## gzip tsv files
        ls \${outdir}/*.tsv | \
        	parallel -j0 gzip {}
	done
	fi

	# manner2
	if [[ ${params.phase_manner2} == true ]] ; then
    ## deal with hpsplit+meth2bed+nanomethphase_vis_bam
	## TODO: change hmc filename in Megalodon raw output
	toolList=("megalodon" "nanome_${params.NANOME_MODEL}")
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
					 	--score-cutoff ${params.PHASE_meth_score_cutoff} \
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
			done # hapType
		done # chr
	done # i
	fi
	"""
}
