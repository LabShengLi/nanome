#!/bin/bash

## Config file of datasets, all defined variables are feed into running script later

dsname=NA19240

targetNum=6

# Input dir, note: do not use filename right now
inputDataDir=/fastscratch/liuya/nanocompare/NA19240

# Output
untarDir=/fastscratch/liuya/nanocompare/NA19240_untar
septDir=/fastscratch/liuya/nanocompare/NA19240_sept
basecalledDir=/fastscratch/liuya/nanocompare/NA19240_basecalled
methCalledDir=/fastscratch/liuya/nanocompare/NA19240_methcalled


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