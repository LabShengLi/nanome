/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : EVAL.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// Check all tools work well
process EVAL {
	tag "${params.dsname}"

 	// debug true

	publishDir "${params.outdir}/${params.dsname}-methylation-callings/Evaluation-${params.dsname}",
		mode: "copy"

	input:
	path tools_read_unify
	path bg_list
	path ch_src
	path ch_utils
	path ch_ga

	output:
	path "MethPerf-*",	emit: read_eval, optional: true
	path "MethCorr-*",	emit: site_eval, optional: true

	when:
	(bg_list.size() >= 2) && (params.runEval)

	shell:
	// construct bgtruth params
	bgParams = bg_list[1..(bg_list.size()-1)].join(";")

	// construct genome annotation dir params
	if (params.genome_annotation_dir) {
		gaParams = "--genome-annotation  ${ch_ga}"
	} else {
		gaParams = ""
	}

	// construct call tools params
	callParams = ""
	for (int i = 1; i < tools_read_unify.size(); i++) {
		infn = tools_read_unify[i].toString()
		cutoff_str = null
		if (infn.toLowerCase().contains('nanopolish')) {
			tool = "Nanopolish"
			if (params.llr_cutoff_nanopolish) {
				cutoff_str = params.llr_cutoff_nanopolish
			}
		}
		if (infn.toLowerCase().contains('megalodon')) {
			tool = "Megalodon"
			if (params.llr_cutoff_megalodon) {
				cutoff_str = params.llr_cutoff_megalodon
			}
		}
		if (infn.toLowerCase().contains('deepsignal')) {
			tool = "DeepSignal"
			if (params.llr_cutoff_deepsignal) {
				cutoff_str = params.llr_cutoff_deepsignal
			}
		}
		if (infn.toLowerCase().contains('nanome')) {
			tool = "NANOME"
			if (params.llr_cutoff_nanome) {
				cutoff_str = params.llr_cutoff_nanome
			}
		}

		param1 = "${tool}:UNIREAD:${infn}"
		if (cutoff_str) {
			param1 = param1 + cutoff_str
		}
		// println(param1)
		callParams += " "+param1
	}
	// print(callParams)
	'''
	date; hostname; pwd

	PYTHONPATH=src  src/nanome/nanocompare/read_level_eval.py \
		--calls \
		!{callParams} \
		--bgtruth "!{params.bg_encode}:!{bgParams}" \
		--runid MethPerf-!{params.dsname}_!{params.runidSuffix} \
		--dsname !{params.dsname} \
		--min-bgtruth-cov !{params.min_bgtruth_cov} \
		--processors !{task.cpus} \
		!{gaParams} !{params.readEvalOptions ? params.readEvalOptions : " "}\
		-o .  &>> !{params.dsname}.Eval.run.log
	echo "### Read level eval DONE"

	PYTHONPATH=src  src/nanome/nanocompare/site_level_eval.py \
		--calls \
				!{callParams} \
		--bgtruth "!{params.bg_encode}:!{bgParams}" \
		--runid MethCorr-!{params.dsname}_!{params.runidSuffix} \
		--dsname !{params.dsname} \
		--min-bgtruth-cov !{params.min_bgtruth_cov} \
		--toolcov-cutoff !{params.toolcov_cutoff} \
		--processors !{task.cpus} \
		!{gaParams} !{params.siteEvalOptions ? params.siteEvalOptions : " "}\
		-o .    &>> !{params.dsname}.Eval.run.log
	echo "### Site level eval DONE"
	echo "### Eval all DONE"
	'''
}
