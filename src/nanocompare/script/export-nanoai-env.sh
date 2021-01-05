#!/bin/bash
set -x

source /home/liuya/.bash_profile

conda activate nanoai

conda env export > nanoai-environment.yml

conda deactivate