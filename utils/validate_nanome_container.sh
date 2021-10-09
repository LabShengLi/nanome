#!/bin/bash
printf "### Nanome container check start\n"
printf "====================\n"
which conda
which python
python --version
which pip
echo "### PATHONPATH=$PATHONPATH"
echo "### PATH=$PATH"

printf "====================\n"
guppy_basecaller -v | head -1

printf "====================\n"
nanopolish --version | head -1

printf "====================\n"
megalodon -v
echo "### we need use tqdm >= 4.60 due to megalodon"
pip show tqdm | tail -n +1 | head -2 | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
pip show deepsignal | tail -n +1 | head -2 | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
fast5mod --version

# fast5mod and megalodon confilicts with it, but we are sure it is ok
# ont-fast5-api >= 3.0
echo "### we need use ont-fast5-api >= 3.0"
pip show ont-fast5-api | tail -n +1 | head -2 | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
tombo -v

# Tombo depend on < 3.0
echo "### we need use h5py < 3.0 due to Tombo"
pip show h5py | tail -n +1 | head -2 | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
pip show deepmod | tail -n +1 | head -2 | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
# METEORE depends on ==0.21.3
echo "### we need use pip install -U scikit-learn==0.21.3 due to METEORE"
pip show scikit-learn | tail -n +1 | head -2 | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
pip show nanome-jax | tail -n +1 | head -2 | awk -vRS="" -vOFS='  ' '$1=$1'

printf "### Nanome container check end\n"
