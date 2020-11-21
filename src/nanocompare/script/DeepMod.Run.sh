#!/bin/bash
#####################################
#
# DeepMod meth call on basecalled data folder
#
# Args:
#   $1  -   configuration file
#
# Usage:
#   /projects/liuya/workspace/tcgajax/nanocompare/script/DeepMod.Run.sh NA19240.Config.sh
#####################################
set -u

source $1

MethScript=/projects/liuya/workspace/tcgajax/nanocompare/script/DeepMod.Array.sbatch

sbatch  --job-name=DeepMod.${dsname} --array=1-${targetNum} ${MethScript} ${basecalledDir} ${methCalledDir}

echo "### Methylation call, 1 job array (${targetNum} jobs) submitted. ###"

