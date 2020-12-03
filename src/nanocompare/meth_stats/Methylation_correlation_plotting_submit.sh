#!/bin/bash

prj_dir=/projects/li-lab/yang/workspace/nano-compare

pythonFile=${prj_dir}/src/nanocompare/meth_stats/Methylation_correlation_plotting_submit.py
inputTsvFile=/projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/NanoComareCorrelation_paper.tsv

mkdir -p log

python ${pythonFile} ${inputTsvFile}
