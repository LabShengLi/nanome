#!/bin/bash

pythonFile=nanocompare/tss_eval_submit.py

# input file
inputTsvFile=${1:-site_level_meteore_analysis_in_paper.tsv}
scriptName=${3:-"meteore_output.sbatch"}
otherOptions=${3:-"--output-meteore"}
#otherOptions=${2:-"--enable-cache --using-cache"}

PYTHONPATH=. python ${pythonFile} ${inputTsvFile} ${scriptName} ${otherOptions}
