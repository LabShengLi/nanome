#!/bin/bash

#####################################
# Config file of datasets,
# all defined variables are feed into running script later
#####################################

# data set name, for showing in job names
dsname=NA19240

# number of seperated folders
targetNum=100

# Input dir, note: do not use filename right now
inputDataDir=



# Output for untar, seperate, and basecall, etc.
untarDir=/fastscratch/liuya/nanocompare/NA19240_untar
septDir=/fastscratch/liuya/nanocompare/NA19240_sept
basecalledDir=/fastscratch/liuya/nanocompare/NA19240_basecalled


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