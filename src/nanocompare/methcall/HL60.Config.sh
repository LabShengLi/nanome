#!/bin/bash

## Config file of datasets, all defined variables are feed into running script later

dsname=HL60

targetNum=50

# Input dir, note: do not use filename right now
# /projects/li-lab/NanoporeData/Leukemia_ONT/20180612_180601-18-li-004-GXB01102-002/HL60-Nanopore_GT18-07373.fast5.tar
inputDataDir=/projects/li-lab/NanoporeData/Leukemia_ONT/20180612_180601-18-li-004-GXB01102-002

# Output
untarDir=/fastscratch/liuya/nanocompare/${dsname}_untar
septDir=/fastscratch/liuya/nanocompare/${dsname}_sept
basecalledDir=/fastscratch/liuya/nanocompare/${dsname}_basecalled


## Sample usage:
##
## untar and seperate
## /projects/liuya/workspace/tcgajax/nanocompare/script/Preprocessing.sh NA19240.Config.sh
##
## Albacore basecall
## /projects/liuya/workspace/tcgajax/nanocompare/script/Basecall.Albacore.sh NA19240.Config.sh
##
## Guppy basecall
## /projects/liuya/workspace/tcgajax/nanocompare/script/Basecall.Guppy.sh NA19240.Config.sh