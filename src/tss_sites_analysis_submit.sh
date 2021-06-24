#!/bin/bash

pythonFile=nanocompare/tss_eval_submit.py

# input file
inputTsvFile=${1:-site_level_tss_analysis_in_paper.tsv}

#otherOptions=${2:-""}
otherOptions=${2:-"--enable-cache --using-cache"}

PYTHONPATH=. python ${pythonFile} ${inputTsvFile} ${otherOptions}
