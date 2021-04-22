#!/bin/bash

prjBaseDir=/projects/li-lab/yang/workspace/nano-compare
pythonFile=${prjBaseDir}/src/nanocompare/tss_eval_submit.py

# input file
inputTsvFile=${1:-site_level_tss_analysis_in_paper.tsv}

#otherOptions=${2}
otherOptions=${2:-"--enable-cache --using-cache"}

PYTHONPATH=${prjBaseDir}/src python ${pythonFile} ${inputTsvFile} ${otherOptions}
