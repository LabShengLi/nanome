#!/bin/bash

toolList=(Nanopolish Tombo DeepSignal Megalodon Guppy Guppy.gcf52ref DeepMod.C DeepMod.Cluster)


toolList=(DeepMod.C DeepMod.Cluster)

for tool in "${toolList[@]}"; do
	echo "Process combine ${tool}"
	sbatch --job-name=comb.allchrs.NA12878.${tool} combine_output_na12878.sbatch ${tool}
done

exit 0

#bash combine_output_na12878.sh Nanopolish
#bash combine_output_na12878.sh Tombo
#bash combine_output_na12878.sh DeepSignal
#bash combine_output_na12878.sh DeepMod.C
#bash combine_output_na12878.sh DeepMod.Cluster
#bash combine_output_na12878.sh Guppy
#bash combine_output_na12878.sh Guppy.gcf52ref
#bash combine_output_na12878.sh Megalodon
