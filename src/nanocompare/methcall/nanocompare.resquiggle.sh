#!/bin/bash
#SBATCH --job-name=rsquigl
#SBATCH --partition=compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
##SBATCH -p gpu
##SBATCH --gres=gpu:1
##SBATCH -q inference
#SBATCH --mem=150g
#SBATCH --time=1-16:00:00
#SBATCH -o log/%x.%j.figures
#SBATCH -e log/%x.%j.err
##SBATCH --array=1-11

################################################################################
# Tombo resquiggle
# Need to populate the parameters into this script
# Step: Tombo resquiglling
################################################################################

set +e
set +u
set +x
source ../../utils.common.sh
set -x

job_index=$((SLURM_ARRAY_TASK_ID))


## Original basecalled results dir, we do not want resquiggling modify it
jobkBasecallOutputDir=${basecallOutputDir}/${job_index}

## Resquiggle results, firstly copy from basecall, then do resquiggling
jobkResquiggleOutputDir=${resquiggleDir}/${job_index}

## Basecalled fast5 files with base info in it
processedFast5DIR=${jobkResquiggleOutputDir}/workspace

set -u
echo "##################"
echo "dsname: ${dsname}"
echo "targetNum: ${targetNum}"
echo "basecallOutputDir: ${basecallOutputDir}"
echo "jobkBasecallOutputDir: ${jobkBasecallOutputDir}"
echo "jobkResquiggleOutputDir: ${jobkResquiggleOutputDir}"
echo "processedFast5DIR: ${processedFast5DIR}"
echo "refGenome: ${refGenome}"
echo "correctedGroup: ${correctedGroup}"
echo "chromSizesFile: ${chromSizesFile}"
echo "processors: ${processors}"
echo "##################"
set +u

set +x
conda activate nanoai
set -x

## Re-squiggle the data:
## TODO: Only can be done once currently? Why --overwrite is not allow multiple runs? Ref https://github.com/nanoporetech/tombo/issues/5

# See ref of tombo resquiggle by tombo resquiggle -h

# Firstly clean it, then copy basecall results into resquiggle dir
rm -rf ${jobkResquiggleOutputDir}
mkdir -p ${jobkResquiggleOutputDir}

cp -rf ${jobkBasecallOutputDir}/* ${jobkResquiggleOutputDir}/

time tombo resquiggle --dna --processes ${processors} --corrected-group ${correctedGroup} \
	--basecall-group Basecall_1D_000 --overwrite ${processedFast5DIR} ${refGenome} #--quiet


set +x
conda deactivate
set -x

echo "###   Re-squiggling by Tombo DONE    ###"
