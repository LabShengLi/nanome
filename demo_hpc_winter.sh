#!/bin/bash
#SBATCH --job-name=nanome.demo
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -q inference
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem=20G # memory pool for all cores
#SBATCH --time=06:00:00 # time
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

date; hostname; pwd

outDir=${1:-/fastscratch/li-lab/nanome}
mkdir -p ${outDir}; chmod ugo+w ${outDir}

########################################
########################################
# Get nextflow and install it
if [ ! -f "nextflow" ]; then
    curl -s https://get.nextflow.io | bash
fi

########################################
########################################
# Pull nanome singularity
module load singularity

export SINGULARITY_CACHEDIR="${outDir}/singularity-cache"
mkdir -p  $SINGULARITY_CACHEDIR; chmod ugo+w $SINGULARITY_CACHEDIR
sif_dir="${outDir}/sif"
mkdir -p $sif_dir; chmod ugo+w $sif_dir
nanome_singularity="${sif_dir}/nanome_v1.4.sif"

if [ ! -f ${nanome_singularity} ]; then
    singularity pull ${nanome_singularity} docker://quay.io/liuyangzzu/nanome:v1.4
fi

rm -rf ${outDir}/work
rm -rf ${outDir}/outputs

mkdir -p ${outDir}/work; chmod ugo+w ${outDir}/work
mkdir -p ${outDir}/outputs; chmod ugo+w ${outDir}/outputs

########################################
########################################
# Running pipeline
set -x
./nextflow run main.nf \
    -profile winter2 -resume \
    -with-report -with-timeline -with-trace -with-dag \
    -work-dir ${outDir}/work \
    -with-singularity ${nanome_singularity} \
    --input https://raw.githubusercontent.com/liuyangzzu/nanome/master/inputs/test.demo.filelist.txt \
    --outputDir ${outDir}/outputs \
    --eval true
