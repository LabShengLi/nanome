#!/usr/bin/env python3

import csv
import os
import subprocess
from sys import argv

from nanocompare.global_config import src_base_dir


if __name__ == '__main__':

    infile = open(argv[1], 'r')
    scriptName = argv[2] # "tss_eval.sbatch"
    scriptFileName = os.path.join(src_base_dir, "nanocompare", scriptName)

    others = ' '.join(argv[3:])
    print(f'Other options={others}')

    csvfile = csv.DictReader(infile, delimiter='\t')
    for row in csvfile:
        if row['status'] == "submit":
            outlogdir = os.path.join('.', 'log')
            os.makedirs(outlogdir, exist_ok=True)

            command = f"""
set -x; 

sbatch --job-name=tss-eval-{row['RunPrefix']} --output=log/%x.%j.out --error=log/%x.%j.err \
--export=ALL,dsname="{row['Dataset']}",DeepSignal_calls="{row['DeepSignal_calls']}",\
Tombo_calls="{row['Tombo_calls']}",Nanopolish_calls="{row['Nanopolish_calls']}",\
DeepMod_calls="{row['DeepMod_calls']}",Megalodon_calls="{row['Megalodon_calls']}",Guppy_calls="{row['Guppy_calls']}",\
bgTruth="{row['bgTruth']}",RunPrefix="{row['RunPrefix']}",parser="{row['parser']}",\
otherOptions="{others}" {scriptFileName}

echo DONE
"""
            print(f"RunPrefix={row['RunPrefix']}")

            # output sbatch submit a job's results to STDOUT
            print(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8"))
