#!/bin/bash

pythonFile=Methylation_correlation_plotting_submit.py
inputTsvFile=/projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComareCorrelation_configFile_NA19240.tsv

python ${pythonFile} ${inputTsvFile}
