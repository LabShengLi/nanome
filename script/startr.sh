#!/bin/bash

#module load singularity

#singularity pull --name singularity-r.simg shub://nickjer/singularity-r

#singularity run --app Rscript singularity-r.simg test.R

singularity run -w singularity-r.simg

