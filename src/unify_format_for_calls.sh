#!/bin/bash
# @Author   : Yang Liu
# @FileName : unify_format_for_calls.sh
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

# Generate unified read-level and/or site-level format of calls
# Usage:
# [prog]  <dsname> <call-encode> <call-fn> <outd-dir> <num-processors> <step12> <chr-options>
#             1       2               3       4               5            6          7
dsname=${1}
encode=${2}
callfn=${3}
outdir=${4}
processors=${5}
step=${6:-'12'}
chr_options=${7:-''}

if [[ "${chr_options}" != "" ]] ; then
		chr_options="--chrSet ${chr_options}"
fi

if [[ "$step" == *"1"* ]]; then
    ## Read level unify
    echo "### Read level unify"
    tss_eval.py \
        --calls \
            ${encode}:${callfn} \
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
            ${encode}:${callfn} \
        --runid Site_Level-${dsname} \
        --dsname ${dsname} \
        --processors ${processors}	\
        -o ${outdir}   ${chr_options}

    ## Sort site level results
    fnlist=$(find Site_Level-${dsname} -name '*-perSite-cov1.bed.gz')
    for fn in $fnlist ; do
        outfn=${fn/cov1.bed.gz/cov1.sort.bed.gz}
        bedtools sort -i ${fn} | gzip -f > ${outfn}
        rm ${fn}
    done
fi

echo "### Unify calls DONE for ${encode}:${callfn}"
