#!/bin/bash
# @Author   : Yang Liu
# @FileName : validate_nanome_container.sh
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

#versionFilename=${1:-"tools_version_table.tsv"}
versionFilename=${1:-}

printf "### Nanome container check start\n"
printf "====================\n"
which conda
which python
python --version
which pip
echo "### PATHONPATH=$PYTHONPATH"
echo "### PATH=$PATH"

printf "====================\n"
echo "### Check Guppy dir:"
which guppy_basecaller
guppy_basecaller -v | head -1
guppy_version=$(guppy_basecaller -v |  head -n 1 | getGuppyVersion.py)

printf "====================\n"
nanopolish --version | head -1
nanopolish_version=$(nanopolish --version | head -1 | awk '{print $NF}')

printf "====================\n"
megalodon -v
megalodon_version=$(megalodon -v | head -1 | awk '{print $NF}')

echo "### we need use tqdm >= 4.60 due to megalodon"
pip show tqdm | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
pip show deepsignal | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1'
deepsignal_version=$(pip show deepsignal | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1' | awk '{print $NF}')

printf "====================\n"
fast5mod --version
fast5mod_version=$(fast5mod --version | head -n 1 | awk '{print $NF}')

# fast5mod and megalodon confilicts with it, but we are sure it is ok
# ont-fast5-api >= 3.0
echo "### we need use ont-fast5-api >= 3.0"
pip show ont-fast5-api | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1'
ont_fast5_api_version=$(pip show ont-fast5-api | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1' | awk '{print $NF}')

printf "====================\n"
tombo -v
tombo_version=$(tombo -v | head -n 1 | awk '{print $NF}')

# Tombo depend on h5py < 3.0
echo "### we need use h5py < 3.0 due to Tombo"
pip show h5py | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
pip show deepmod | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1'
deepmod_version=$(pip show deepmod | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1' | awk '{print $NF}')

printf "====================\n"
# METEORE depends on ==0.21.3 <=0.23.2
echo "### we need use pip install -U scikit-learn==0.21.3 <=0.23.2 due to METEORE"
pip show scikit-learn | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
echo inliner version = `inliner -V`
pip show nanome-jax | awk 'NR<=2' | awk -vRS="" -vOFS='  ' '$1=$1'

printf "====================\n"
printf '### Tools version list\n'
if [[ -z "${versionFilename}" ]] ; then
    printf '%s\t%s\n' Tool Version
    printf '%s\t%s\n' Nanopolish ${nanopolish_version} 
    printf '%s\t%s\n' Megalodon ${megalodon_version} 
    printf '%s\t%s\n' DeepSignal ${deepsignal_version} 
    printf '%s\t%s\n' Guppy ${guppy_version} 
    printf '%s\t%s\n' Tombo ${tombo_version} 
    printf '%s\t%s\n' DeepMod ${deepmod_version} 
    printf '%s\t%s\n' METEORE 1.0.0 
else
    > ${versionFilename}
    printf '%s\t%s\n' Tool Version >> ${versionFilename}
    printf '%s\t%s\n' NANOME 1.0 >> ${versionFilename}
    printf '%s\t%s\n' Nanopolish ${nanopolish_version} >> ${versionFilename}
    printf '%s\t%s\n' Megalodon ${megalodon_version} >> ${versionFilename}
    printf '%s\t%s\n' DeepSignal ${deepsignal_version} >> ${versionFilename}
    printf '%s\t%s\n' Guppy ${guppy_version} >> ${versionFilename}
    printf '%s\t%s\n' Tombo ${tombo_version} >> ${versionFilename}
    printf '%s\t%s\n' METEORE 1.0.0 >> ${versionFilename}
    printf '%s\t%s\n' DeepMod ${deepmod_version} >> ${versionFilename}

    echo "### check tools version file:${versionFilename}"
    cat ${versionFilename}
fi
printf "### Nanome container check end\n"
