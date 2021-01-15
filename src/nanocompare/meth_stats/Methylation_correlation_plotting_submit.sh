#!/bin/bash

prjBaseDir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${prjBaseDir}/src/nanocompare/meth_stats/Methylation_correlation_plotting_submit.py

# input file
inputTsvFile=${prjBaseDir}/src/nanocompare/meth_stats/NanoComareCorrelation_paper.tsv

PYTHONPATH=${prjBaseDir}/src python ${pythonFile} ${inputTsvFile}
