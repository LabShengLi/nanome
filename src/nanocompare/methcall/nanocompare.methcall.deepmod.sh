#!/bin/bash
#SBATCH --job-name=deepsignal.methcall
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem=250g # memory pool for all cores
#SBATCH --time=2-23:00:00
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
##SBATCH --array=1-11

################################################################################
# DeepMod methylation call workflow
# Need to populate the parameters into this script
################################################################################
#cd "$(dirname "$0")"

set -e
set +x
source ../../utils.common.sh

set -x

job_index=$((SLURM_ARRAY_TASK_ID))

## Modify directory for processed files after basecalling:
processedFast5DIR=${basecallOutputDir}/${job_index}/workspace

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "Tool: ${Tool}"
echo "targetNum: ${targetNum}"
echo "analysisPrefix: ${analysisPrefix}"
echo "basecall_name: ${basecall_name}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "processedFast5DIR: ${processedFast5DIR}"
echo "methCallsDir: ${methCallsDir}"
echo "refGenome: ${refGenome}"
echo "run_resquiggling: ${run_resquiggling}"
echo "deepModModel: ${deepModModel}"
echo "isGPU: ${isGPU}"
echo "processors: ${processors}"
echo "##################"
set +u

set +x
conda activate nanoai
set -x

rm -rf ${methCallsDir}/batch_${job_index} ${methCallsDir}/batch_${job_index}.done

## Call methylation from processed fast5 files:
## For guppy results, add --move for correct Events data location in fast5 files. ref:https://github.com/WGLab/DeepMod/issues/29#issuecomment-594121403

if [ "${basecall_name}" = "Albacore" ] ; then # Albacore using event data
	time python ${DeepModDir}/bin/DeepMod.py detect \
				--wrkBase ${processedFast5DIR} --Ref ${refGenome} --outFolder ${methCallsDir} \
				--Base C --modfile ${deepModModel} --FileID batch_${job_index} \
				--threads ${processors}
elif [ "${basecall_name}" = "Guppy" ] ; then # Guppy using move data
	time python ${DeepModDir}/bin/DeepMod.py detect \
			--wrkBase ${processedFast5DIR} --Ref ${refGenome} --outFolder ${methCallsDir} \
			--Base C --modfile ${deepModModel} --FileID batch_${job_index} \
			--threads ${processors} --move
fi
echo "###   DeepMod methylation calling DONE"

set +x
conda deactivate
set -x

echo "###   DeepMod Meth-call task is DONE    ###"
