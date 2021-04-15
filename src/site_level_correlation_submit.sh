#!/bin/bash

prjBaseDir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${prjBaseDir}/src/nanocompare/site_level_eval_submit.py

# input file
inputTsvFile=correlation_of_site_level_in_paper.tsv

# --beddir is dir for read-level results basedir, used for concordant and discordant sites BED files
PYTHONPATH=${prjBaseDir}/src python ${pythonFile} ${inputTsvFile} --beddir /projects/li-lab/yang/results/2021-04-12 \
	--enable-cache --using-cache
