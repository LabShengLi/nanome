#!/bin/bash

for i in {1..22} X Y; do
	echo "sanity check chr${i}"
	sbatch --job-name=reproduce_deepmod_chr${i}  reproduce_deepmod_na12878.sbatch chr${i}
done
