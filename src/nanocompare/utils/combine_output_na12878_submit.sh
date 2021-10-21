#!/bin/bash

toolList=(Nanopolish Tombo DeepSignal Megalodon Guppy Guppy.gcf52ref DeepMod.C DeepMod.Cluster)


toolList=(DeepMod.C DeepMod.Cluster)

for tool in "${toolList[@]}"; do
	echo "Process combine ${tool}"
	sbatch --job-name=comb.allchrs.NA12878.${tool} combine_output_na12878.sbatch ${tool}
done


