#!/bin/bash
#SBATCH --job-name=nanome.demo.hpc
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q training
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=20G # memory pool for all cores
#SBATCH --time=14-00:00:00 # time
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

date; hostname; pwd

baseDir=${1:-/fastscratch/li-lab/nanome}

workDir=${baseDir}/work
outputsDir=${baseDir}/outputs


########################################
########################################
# Ensure directories
export SINGULARITY_CACHEDIR="${baseDir}/singularity-cache"
mkdir -p  $SINGULARITY_CACHEDIR; chmod ugo+w $SINGULARITY_CACHEDIR

########################################
########################################
# Get nextflow and install it
if [ ! -f "nextflow" ]; then
    curl -s https://get.nextflow.io | bash
fi


########################################
# Clean old results
#rm -rf ${workDir} ${outputsDir}
mkdir -p ${workDir}; chmod ugo+w ${workDir}
mkdir -p ${outputsDir}; chmod ugo+w ${outputsDir}


########################################
########################################
# Running pipeline for demo human data
# More options: -with-report -with-timeline -with-trace -with-dag -resume
module load singularity
set -x
./nextflow run main.nf \
    -profile winter_singularity \
    -with-report -with-timeline -with-trace -with-dag -resume \
    -work-dir ${workDir} \
    --outputDir ${outputsDir} \
    --dsname TestData \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt

exit 0


########################################
########################################
# Running pipeline for demo human data
# More options: -with-report -with-timeline -with-trace -with-dag -resume
module load singularity
set -x
./nextflow run main.nf \
    -profile winter_singularity \
    --dsname TestData \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt

exit 0

./nextflow run https://github.com/liuyangzzu/nanome.git \
    -profile winter_singularity \
    --dsname TestData \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt

exit 0