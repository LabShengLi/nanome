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
rm -rf ${workDir} ${outputsDir}
mkdir -p ${workDir}; chmod ugo+w ${workDir}
mkdir -p ${outputsDir}; chmod ugo+w ${outputsDir}


########################################
########################################
# Running pipeline for demo human data
# More options: -with-report -with-timeline -with-trace -with-dag -resume

module load singularity
set -x
./nextflow run main.nf \
    -profile singularity,hpc \
    -work-dir ${workDir} \
    --outputDir ${outputsDir} \
    --dsname TestData \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt
exit 0



module load singularity
set -x
./nextflow run main.nf \
    -profile singularity,hpc \
    -config conf/jax_hpc.config \
    -work-dir ${workDir} \
    --outputDir ${outputsDir} \
    --dsname TestData \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt
exit 0


./nextflow run main.nf \
    -profile conda,hpc \
    -config conf/jax_hpc.config \
    -work-dir ${workDir} \
    --outputDir ${outputsDir} \
    --dsname TestData \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt \
    --guppyDir '/projects/li-lab/software/ont-guppy-gpu_4.2.2'

# Report
tree ${workDir} > ${baseDir}/work.tree.txt
tree ${outputsDir} > ${baseDir}/outputs.tree.txt
echo "### nanome pipeline demo DONE"
exit 0

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
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt \
    --bgTruth='/projects/li-lab/Nanopore_compare/data/NA19240/NA19240_RRBS_ENCFF000LZS_rep1.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz;/projects/li-lab/Nanopore_compare/data/NA19240/NA19240_RRBS_ENCFF000LZT_rep2.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz' \
    --eval true

# Report
tree ${workDir} > ${baseDir}/work.tree.txt
tree ${outputsDir} > ${baseDir}/outputs.tree.txt
echo "### nanome pipeline demo DONE"
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
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt \
    --processors 4

exit 0