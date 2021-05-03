#!/bin/bash

pythonFile=nanocompare/site_level_eval_submit.py

# input file
inputTsvFile=correlation_of_site_level_in_paper.tsv

# --beddir is dir for read-level results basedir, used for concordant and discordant sites BED files
PYTHONPATH=. python ${pythonFile} ${inputTsvFile} \
		--beddir /projects/li-lab/yang/results/2021-04-16 \
		--enable-cache --using-cache
