#!/home/liuya/anaconda3/envs/nmf/bin/python

"""
Submit a list of jobs for running input of methylation_plotting script, for each line of tsv input, the output will be a tsv file with 4 tools and bgtruth methylation percentage and coverage info.

python Methylation_correlation_plotting_submit.py /projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComareCorrelation_paper.tsv

Methylation_correlation_plotting_submit.py NanoComareCorrelation_paper.tsv
"""
# example run command: python Methylation_correlation_plotting_submit.py <config file>
# python /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/Methylation_correlation_plotting_submit.py NanoComareCorrelation_deprecated.tsv
# python /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/Methylation_correlation_plotting_submit.py NanoComareCorrelation_paper.tsv

import csv
import os
import subprocess
from sys import argv

from nanocompare.global_config import prjsrc_base_dir

pyFile = os.path.join(prjsrc_base_dir, "nanocompare/meth_stats", "Methylation_correlation_plotting.sbatch")

if __name__ == '__main__':

    infile = open(argv[1], 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')
    for row in csvfile:
        if row['status'] == "submit":
            outlogdir = os.path.join('.', 'log')
            os.makedirs(outlogdir, exist_ok=True)

            command = f"""
set -x; sbatch --output=log/%x.%j.out --error=log/%x.%j.err --export=ALL,DeepSignal_calls="{row['DeepSignal_calls']}",Tombo_calls="{row['Tombo_calls']}",Nanopolish_calls="{row['Nanopolish_calls']}",DeepMod_calls="{row['DeepMod_calls']}",Megalodon_calls="{row['Megalodon_calls']}",bgTruth="{row['bgTruth']}",RunPrefix="{row['RunPrefix']}",parser="{row['parser']}" --job-name=meth-corr-{row['RunPrefix']} {pyFile}"""

            # print(command)
            print(f"RunPrefix={row['RunPrefix']}")

            # output sbatch submit a job's results to STDOUT
            print(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8"))

            # print(row['RunPrefix'], subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read())
