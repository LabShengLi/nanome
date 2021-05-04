#!/bin/bash

pythonFile=nanocompare/site_level_eval_submit.py

# input file
inputTsvFile=${1:-correlation_of_site_level_in_paper.tsv}
beddirOptions=${2:-"--beddir /projects/li-lab/Nanopore_compare/result/meth-exp/MethPerf-cut5"}
otherOptions=${3:-""}
#otherOptions=${3-:"--enable-cache --using-cache"}

# --beddir is dir for read-level results basedir, used for concordant and discordant sites BED files
PYTHONPATH=. python ${pythonFile} ${inputTsvFile} \
    ${beddirOptions}  ${otherOptions}\

