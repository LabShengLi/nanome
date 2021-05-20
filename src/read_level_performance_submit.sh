#!/bin/bash

pythonFile=nanocompare/read_level_eval_submit.py

# input file
inputTsvFile=${1:-performance_of_read_level_in_paper.tsv}
##otherOptions=${2:-"-mpi"}
otherOptions=${2:-"-mpi --enable-cache --using-cache"}

PYTHONPATH=. python ${pythonFile} ${inputTsvFile}  ${otherOptions}