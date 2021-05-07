#!/bin/bash

#module load singularity

#singularity pull --name singularity-r.simg shub://nickjer/singularity-r

#singularity run --app Rscript singularity-r.simg test.R

module load singularity

baseDir=/projects/li-lab/yang/
singularity run -w ${baseDir}/singularity-r.simg

