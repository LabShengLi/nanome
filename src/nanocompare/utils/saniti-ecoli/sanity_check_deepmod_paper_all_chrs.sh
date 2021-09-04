#!/bin/bash

for i in {1..22} X Y; do
	echo "sanity check chr${i}"
	sbatch --job-name=sanity_check_deepmod_chr${i}  sanity_check_deepmod_paper.sbatch chr${i}
done
