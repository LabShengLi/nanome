#!/bin/bash

## Config file of datasets, all defined variables are feed into running script later

dsname=APL

targetNum=100

# Input dir, note: do not use filename right now
# /projects/li-lab/AML-Nanopore/20180517_180508-18-li-001-GXB01186-001/APL-1750_GT18-06409.fast5.tar
inputDataDir=/projects/li-lab/AML-Nanopore/20180517_180508-18-li-001-GXB01186-001

# Output
untarDir=/fastscratch/liuya/nanocompare/APL_untar
septDir=/fastscratch/liuya/nanocompare/APL_sept
basecalledDir=/fastscratch/liuya/nanocompare/APL_basecalled


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