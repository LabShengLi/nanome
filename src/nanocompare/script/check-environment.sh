#!/bin/bash

set -e

source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/conda_setup.sh

export PATH=/projects/li-lab/yang/tools/nanopolish:${PATH}

conda activate nanoai

# Check Albacore path
read_fast5_basecaller.py --version

# Check nanopore tool path
tombo --version

python /projects/li-lab/yang/tools/DeepMod/bin/DeepMod.py

deepsignal

pip show deepsignal

nanopolish --version


echo "### Enviroment check ok."