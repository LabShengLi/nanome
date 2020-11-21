#!/bin/bash

#####################################
#
# Guppy basecall sept data folder
#
# Args:
#   $1  -   configuration file
#
# Usage:
#   /projects/liuya/workspace/tcgajax/nanocompare/script/Basecall.Guppy.sh NA19240.Config.sh
#####################################

set -u

source $1

BasecallScript=/projects/liuya/workspace/tcgajax/nanocompare/script/Basecall.Guppy.Array.sbatch

sbatch --job-name=Guppy.${dsname} --array=1-${targetNum} ${BasecallScript} ${septDir} ${basecalledDir}

echo "### Basecall Guppy, 1 job array (${targetNum} jobs) submitted. ###"
