#!/bin/bash

prjBaseDir=/projects/li-lab/yang/workspace/nano-compare
pythonFile=${prjBaseDir}/src/nanocompare/meth_stats/Universal_meth_stats_evaluation_submit.py

# input file
inputTsvFile=NanoComarePerformance_paper.tsv

#Universal_meth_stats_evaluation_submit.py NanoComarePerformance_paper.tsv

PYTHONPATH=${prjBaseDir}/src python ${pythonFile} ${inputTsvFile}