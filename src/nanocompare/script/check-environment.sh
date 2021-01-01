#!/bin/bash

set -e

source /home/liuya/.bash_profile

export PATH=/projects/li-lab/yang/tools/nanopolish:${PATH}

conda activate nanoai

echo "### Conda enviroment nanoai"
conda info --envs

echo "### Albacore"
# Check Albacore path
read_fast5_basecaller.py --version

echo "### Tombo"
# Check nanopore tool path
tombo --version

echo "### DeepMod"
echo python /projects/li-lab/yang/tools/DeepMod/bin/DeepMod.py
python /projects/li-lab/yang/tools/DeepMod/bin/DeepMod.py

echo "### DeepSignal"
deepsignal

pip show deepsignal

echo "### Nanopolish"
nanopolish --version


echo "### Enviroment check ok."

conda deactivate