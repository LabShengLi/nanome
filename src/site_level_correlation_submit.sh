#!/bin/bash

prjBaseDir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${prjBaseDir}/src/nanocompare/site_level_eval_submit.py

# input file
inputTsvFile=correlation_of_site_level_in_paper.tsv

PYTHONPATH=${prjBaseDir}/src python ${pythonFile} ${inputTsvFile}
