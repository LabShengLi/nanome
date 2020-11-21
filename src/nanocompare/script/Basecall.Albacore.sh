#!/bin/bash


#####################################
#
# Albacore sept data folder
#
# Args:
#   $1  -   configuration file
#
# Usage:
#   /projects/liuya/workspace/tcgajax/nanocompare/script/Basecall.Albacore.sh NA19240.Config.sh
#####################################

set -u

source $1

BasecallScript=/projects/liuya/workspace/tcgajax/nanocompare/script/Basecall.Albacore.Array.sbatch

sbatch  --job-name=Albacore.${dsname} --array=1-${targetNum} ${BasecallScript} ${septDir} ${basecalledDir}

echo "### Basecall Albacore, 1 job array (${targetNum} jobs) submitted. ###"
