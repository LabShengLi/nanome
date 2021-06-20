#!/bin/bash
#SBATCH --job-name=meth.perf.NA12878_WGBS_2Reps
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem=280g
#SBATCH --time=72:00:00
#SBATCH -o log/%x.%j.out
#SBATCH -e log/%x.%j.err

set -x

bgTruth="/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz;/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz"
parser="encode"
baseDir="/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs"

DeepModFileName=$(find ${baseDir} -type f -name "*.deepmod.C.combine_allchrs.*gz")
DeepModClusterFileName=$(find ${baseDir} -type f -name "*.deepmod.C_clusterCpG.combine_allchrs.*.gz")
NanopolishFileName=$(find ${baseDir} -type f -name "*.nanopolish.*.combine_allchrs.*.gz")
DeepSignalFileName=$(find ${baseDir} -type f -name "*.deepsignal.*.combine_allchrs.*.gz")
TomboFileName=$(find ${baseDir} -type f -name "*.tombo.*.combine_allchrs.*.gz")
MegalodonFileName=$(find ${baseDir} -type f -name "*.megalodon.*.combine_allchrs.*.gz")
GuppyFast5modFileName=$(find ${baseDir} -type f -name "*.guppy.fast5mod*.combine_allchrs.*.gz")
GuppyGcf52refFileName=$(find ${baseDir} -type f -name "*.guppy.gcf52ref*.combine_allchrs.*.gz")

echo "########################################"
echo DeepModFileName=${DeepModFileName}
echo DeepModClusterFileName=${DeepModClusterFileName}
echo NanopolishFileName=${NanopolishFileName}
echo DeepSignalFileName=${DeepSignalFileName}
echo TomboFileName=${TomboFileName}
echo MegalodonFileName=${MegalodonFileName}
echo GuppyFast5modFileName=${GuppyFast5modFileName}
echo GuppyGcf52refFileName=${GuppyGcf52refFileName}
echo "########################################"

dsname=NA12878
RunPrefix="NA12878_WGBS_2Reps_Guppy_both"
minCov=5
otherOptions="--enable-cache --using-cache"

pythonFn=nanocompare/read_level_eval.py
#PYTHONPATH=. python ${pythonFn} --calls DeepSignal:${DeepSignalFileName} \
#	Tombo:${TomboFileName} Nanopolish:${NanopolishFileName} \
#	Megalodon:${MegalodonFileName} Guppy:${GuppyFast5modFileName} \
#	--bgtruth ${parser}:${bgTruth} \
#	--runid MethPerf-${RunPrefix} \
#	--dsname ${dsname} --min-bgtruth-cov ${minCov} \
#	--report-joined $otherOptions

python ${pythonFn} --calls DeepSignal:${DeepSignalFileName} \
	Tombo:${TomboFileName} Nanopolish:${NanopolishFileName} \
	Megalodon:${MegalodonFileName} Guppy.gcf52ref:${GuppyGcf52refFileName} \
	Guppy:${GuppyFast5modFileName} \
	--bgtruth ${parser}:${bgTruth} \
	--runid MethPerf-${RunPrefix} \
	--dsname ${dsname} --min-bgtruth-cov ${minCov} \
	--report-joined $otherOptions