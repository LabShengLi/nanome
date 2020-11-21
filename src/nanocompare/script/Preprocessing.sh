#!/bin/bash

#####################################
#
# Preprocessing Tar Dir
#
# 1. untar all *.fast5.tar files;
# 2. seperate into # targetNum of folders (0-N-1). Note: job2 waiting for job1 finished and returned ok.
#
# Args:
#   $1  -   configuration file
#
# Usage:
#   /projects/liuya/workspace/tcgajax/nanocompare/script/Preprocessing.sh NA19240.Config.sh
#
#####################################

set -u

## Load config variables
source $1

#set -x

## Untar input
untarJobID=$(sbatch --job-name=UntarDir.${dsname} /projects/liuya/workspace/tcgajax/nanocompare/script/UntarDir.sbatch ${inputDataDir} ${untarDir})

echo submit untar job: ${untarJobID}

## Seperate input
sbatch --job-name=Sept.${dsname} --dependency=afterok:${untarJobID##* } /projects/liuya/workspace/tcgajax/nanocompare/script/SeperateFast5Files.sbatch ${untarDir} ${septDir} ${targetNum}

echo "### Preprocessing TarDir, 2 jobs submitted. ###"

