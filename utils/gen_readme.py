import sys

from string import Template

with open(sys.argv[1], 'r', encoding="utf-8") as file:
    instr=file.read()

readme_template = Template(instr)

replace_dict = {'dsname':sys.argv[2],'outdir':sys.argv[3]}
readme_out = readme_template.substitute(replace_dict)
print(readme_out)

