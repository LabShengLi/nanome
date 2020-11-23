#!/home/liuya/anaconda3/envs/nmf/bin/python

"""
Submit a list of jobs for running input of methylation_plotting script, for each line of tsv input, the output will be a tsv file with 4 tools and bgtruth methylation percentage and coverage info.

python Methylation_correlation_plotting_submit.py /projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComareCorrelation_paper.tsv

Methylation_correlation_plotting_submit.py NanoComareCorrelation_paper.tsv
"""

"""
modify tombo results in tsv as follows:

/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/NA19240_perChr/bk/NA19240_allChr.bed.CpGs.bed
"""
# example run command: python Methylation_correlation_plotting_submit.py <config file>
# python /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/Methylation_correlation_plotting_submit.py NanoComareCorrelation_deprecated.tsv
# python /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/Methylation_correlation_plotting_submit.py NanoComareCorrelation_paper.tsv


import csv
import os
import subprocess
from sys import argv

if __name__ == '__main__':
    nanocompare_prj = "/projects/li-lab/yang/workspace/nano-compare/src"

    baseDir = os.path.join(nanocompare_prj, "nanocompare/meth_stats")
    scriptFn = "Methylation_correlation_plotting.sbatch"

    infile = open(argv[1], 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')
    for row in csvfile:
        if row['status'] == "submit":
            # print(f"row={row}\n")

            command = """set -x; sbatch --export=DeepSignal_calls="{}",Tombo_calls="{}",Nanopolish_calls="{}",DeepMod_calls="{}",DeepMod_cluster_calls="{}",bgTruth="{}",RunPrefix="{}",parser="{}" --job-name=meth-corr-{} {}/{}""" \
                .format(row['DeepSignal_calls'], row['Tombo_calls'], row['Nanopolish_calls'],
                        row['DeepMod_calls'], row['DeepMod_cluster_calls'], row['bgTruth'],
                        row['RunPrefix'], row['parser'], row['RunPrefix'], baseDir, scriptFn)

            # print(f"command=[{command}]\n")
            print(f"RunPrefix={row['RunPrefix']}")

            # output sbatch submit a job's results to STDOUT
            print(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8"))

            # print(row['RunPrefix'], subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read())
