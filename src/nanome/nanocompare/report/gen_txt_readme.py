#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : gen_readme.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

import sys

from string import Template

with open(sys.argv[1], 'r', encoding="utf-8") as file:
    instr = file.read()

readme_template = Template(instr)

replace_dict = {'dsname': sys.argv[2], 'outdir': sys.argv[3],
                'projectDir': sys.argv[4], 'workDir': sys.argv[5],
                'commandLine': sys.argv[6], 'runName': sys.argv[7],
                'start': sys.argv[8]
                }
readme_out = readme_template.substitute(replace_dict)
print(readme_out)
