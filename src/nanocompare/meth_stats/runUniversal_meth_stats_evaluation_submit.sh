#!/bin/bash

pythonFile=Universal_meth_stats_evaluation_submit.py
inputTsvFile=/projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComarePerformance_NA19240.tsv

python ${pythonFile} ${inputTsvFile}
