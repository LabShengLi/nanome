#!/bin/bash

## Pass all following vars to py script
set -x
set -u
set -e

echo ${NanoCompareDir}


Dataset="NA12878-CHR20"
chrName="chr20"
DeepMod_calls="DeepMod.C:/projects/li-lab/Nanopore_compare/data/NA12878/CHR20/NA12878CHR20.deepmod.C.combine.bed.gz DeepMod.Cluster:/projects/li-lab/Nanopore_compare/data/NA12878/CHR20/NA12878CHR20.deepmod.C_clusterCpG.combine.bed.gz"
bgTruth="/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz;/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz"
parser="encode"
otherOptions="--enable-cache --using-cache"

#sbatch sanity_check_deepmod.sbatch ${Dataset} ${chrName} "${DeepMod_calls}" ${bgTruth} ${parser} 0.5 "--is-report --enable-cache --using-cache"

sbatch sanity_check_deepmod.sbatch ${Dataset} ${chrName} "${DeepMod_calls}" ${bgTruth} ${parser} 0.5 "${otherOptions}"

sbatch sanity_check_deepmod.sbatch ${Dataset} ${chrName} "${DeepMod_calls}" ${bgTruth} ${parser} 0.1 "${otherOptions}"

sbatch sanity_check_deepmod.sbatch ${Dataset} ${chrName} "${DeepMod_calls}" ${bgTruth} ${parser} 0.25 "${otherOptions}"

Dataset="NA12878-CHR2"
chrName="chr2"
DeepMod_calls="DeepMod.C:/projects/li-lab/Nanopore_compare/data/NA12878/CHR2/NA12878-CHR2.deepmod.C.combine.bed.gz DeepMod.Cluster:/projects/li-lab/Nanopore_compare/data/NA12878/CHR2/NA12878-CHR2.deepmod.C_clusterCpG.combine.bed.gz"

sbatch sanity_check_deepmod.sbatch ${Dataset} ${chrName} "${DeepMod_calls}" ${bgTruth} ${parser} 0.5 "${otherOptions}"

Dataset="NA12878-CHR3"
chrName="chr3"
DeepMod_calls="DeepMod.C:/projects/li-lab/Nanopore_compare/data/NA12878/CHR3/NA12878-CHR3.deepmod.C.combine.bed.gz DeepMod.Cluster:/projects/li-lab/Nanopore_compare/data/NA12878/CHR3/NA12878-CHR3.deepmod.C_clusterCpG.combine.bed.gz"

sbatch sanity_check_deepmod.sbatch ${Dataset} ${chrName} "${DeepMod_calls}" ${bgTruth} ${parser} 0.5 "${otherOptions}"


