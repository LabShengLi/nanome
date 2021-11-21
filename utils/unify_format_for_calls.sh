#!/bin/bash
# @Author   : Yang Liu
# @FileName : unify_format_for_calls.sh
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

# Generate unified read-level and/or site-level format of calls
# Usage:
# [prog]  <dsname>   <toolname>  <call-encode> <call-fn> <outd-dir> <num-processors> <step12> <chr-options>
#             1       2               3       4               5            6          7    8
set -x
dsname=${1}
toolname=${2}
encode=${3}
callfn=${4}
outdir=${5}
processors=${6}
step=${7:-'12'}
chr_options=${8:-''}

if [[ "${chr_options}" != "" ]] ; then
	chr_options="--chrSet ${chr_options}"
fi

if [[ "$step" == *"1"* ]]; then
    ## Read level unify
    echo "### Read level unify"
    tss_eval.py \
        --calls \
            ${toolname}:${encode}:${callfn} \
        --runid Read_Level-${dsname} \
        --dsname ${dsname}\
        --read-level-format \
        --processors ${processors}	\
        -o ${outdir}   ${chr_options}
fi

if [[ "$step" == *"2"* ]]; then
    ## Site level unify
    echo "### Site level unify"
    tss_eval.py \
        --calls \
            ${toolname}:${encode}:${callfn} \
        --runid Site_Level-${dsname} \
        --dsname ${dsname} \
        --processors ${processors}	\
        -o ${outdir}   ${chr_options}

    ## Sort site level results
    fnlist=$(find Site_Level-${dsname} -name '*-perSite-cov1.bed.gz')
    for fn in $fnlist ; do
        outfn=${fn/cov1.bed.gz/cov1.sort.bed.gz}
        zcat ${fn} | sort -u -k1,1 -k2,2n | gzip -f > ${outfn}
        rm ${fn}
    done
fi

echo "### Unify calls DONE for ${encode}:${callfn}"
