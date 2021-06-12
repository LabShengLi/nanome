#!/bin/bash

## Pass all following vars to py script
set -x
set -u
set -e

echo ${NanoCompareDir}

chrName=${1:-"chr20"}
predThreshold=${2:-"0.5"}
otherOptions="--enable-cache --using-cache"

chrTagname=${chrName^^}
Dataset="NA12878-${chrTagname}"

bgTruth="/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz;/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz"
parser="encode"
baseDir="/projects/li-lab/Nanopore_compare/data/NA12878"

DeepModFileName=$(find ${baseDir} -type f -name "*${chrTagname}.deepmod.C.combine*gz")
## DeepModFileName=
DeepModClusterFileName=$(find ${baseDir} -type f -name "*${chrTagname}.deepmod.C_clusterCpG.combine.*.gz")
NanopolishFileName=$(find ${baseDir} -type f -name "*${chrTagname}.nanopolish.*.combine.*.gz")
DeepSignalFileName=$(find ${baseDir} -type f -name "*${chrTagname}.deepsignal.*.combine.*.gz")
TomboFileName=$(find ${baseDir} -type f -name "*${chrTagname}.tombo.*.combine.*.gz")
MegalodonFileName=$(find ${baseDir} -type f -name "*${chrTagname}.megalodon.*.combine.*.gz")
GuppyFast5modFileName=$(find ${baseDir} -type f -name "*${chrTagname}.guppy.fast5mod*.combine.*.gz")
GuppyGcf52refFileName=$(find ${baseDir} -type f -name "*${chrTagname}.guppy.gcf52ref*.combine.*.gz")

callList="DeepMod.C:${DeepModFileName} DeepMod.Cluster:${DeepModClusterFileName} Nanopolish:${NanopolishFileName} DeepSignal:${DeepSignalFileName} Tombo:${TomboFileName} Megalodon:${MegalodonFileName} Guppy:${GuppyFast5modFileName} Guppy.gcf52ref:${GuppyGcf52refFileName}"

sbatch --job-name=sanity_check.${chrTagname} sanity_check_deepmod.sbatch ${Dataset} ${chrName} "${callList}" ${bgTruth} ${parser} ${predThreshold} "${otherOptions}"
exit 0

Dataset="NA12878-CHR4"
chrName="chr4"
DeepMod_calls="DeepMod.C:/projects/li-lab/Nanopore_compare/data/NA12878/CHR4/NA12878-CHR4.deepmod.C.combine.bed.gz DeepMod.Cluster:/projects/li-lab/Nanopore_compare/data/NA12878/CHR4/NA12878-CHR4.deepmod.C_clusterCpG.combine.bed.gz"

sbatch sanity_check_deepmod.sbatch ${Dataset} ${chrName} "${DeepMod_calls}" ${bgTruth} ${parser} 0.5 "${otherOptions}"

exit 0

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


