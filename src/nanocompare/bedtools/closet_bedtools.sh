#!/bin/bash
#SBATCH --job-name=bedtools.closest
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 200g # memory pool for all cores
#SBATCH -t 2-23:00:00 # time (D-HH:MM:SS)
#SBATCH -o %x.%j.out # STDOUT
#SBATCH -e %x.%j.err # STDERR

#set -x

data_base_dir=/projects/li-lab/yang/results/11-25/K562_WGBS_Joined

tool_fn_list=(K562_WGBS_Joined-meth-cov-Tombo-baseCount0.bed K562_WGBS_Joined-meth-cov-DeepMod-baseCount0.bed K562_WGBS_Joined-meth-cov-DeepSignal-baseCount0.bed K562_WGBS_Joined-meth-cov-Nanopolish-baseCount0.bed)

bgtruth_fn=K562_WGBS_Joined-meth-cov-bgtruth-baseCount0.bed

bedtools sort -i ${data_base_dir}/${bgtruth_fn} > f2-sorted.bed

for tool_fn in "${tool_fn_list[@]}"
do
    outfn=${tool_fn/-baseCount0.bed/}-bgtruth-closest.bed
    bedtools sort -i ${data_base_dir}/${tool_fn} > f1-sorted.bed
	bedtools closest -a f1-sorted.bed -b f2-sorted.bed -d > ${outfn}
	rm -f f1-sorted.bed
	echo "Save results to ${outfn}"
done

rm -f f2-sorted.bed
echo "Closest task is finished"