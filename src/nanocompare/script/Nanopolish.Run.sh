#!/bin/bash
#####################################
#
# Nanopolish meth call on basecalled data folder
#
# Args:
#   $1  -   configuration file
#
# Usage:
#   /projects/liuya/workspace/tcgajax/nanocompare/script/Nanopolish.Run.sh NA19240.Config.sh
#####################################
set -u

source $1

targetNum=1
MethCallScript=/projects/liuya/workspace/tcgajax/nanocompare/script/Nanopolish.Array.sbatch

sbatch  --job-name=Nanopolish.${dsname} --array=1-${targetNum} ${MethCallScript} ${basecalledDir} ${methCalledDir} ${dsname}

echo "### Nanopolish Methylation call, 1 job array (${targetNum} jobs) submitted. ###"

