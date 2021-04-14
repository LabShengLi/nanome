#!/bin/bash

prjBaseDir=/projects/li-lab/yang/workspace/nano-compare
pythonFile=${prjBaseDir}/src/nanocompare/read_level_eval_submit.py

# input file
inputTsvFile=performance_of_read_level_in_paper.tsv

PYTHONPATH=${prjBaseDir}/src python ${pythonFile} ${inputTsvFile}