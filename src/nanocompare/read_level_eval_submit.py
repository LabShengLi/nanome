#!/home/liuya/anaconda3/envs/nmf/bin/python

"""
Submit evaluation jobs based on tsv of all tools results file name

python Universal_meth_stats_evaluation_submit.py /projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComarePerformance_NA19240.tsv


modify tombo results as

/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/NA19240_perChr/bk/NA19240_allChr.bed.CpGs.bed
"""

# example run command: python UniversalMethStatsEvaluation.standalone_01.submit.py <config file>

# python /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/UniversalMethStatsEvaluation.standalone_02.submit.py NanoComarePerformance_configFile_4programs.tsv
# python /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/UniversalMethStatsEvaluation.standalone_03.submit.py /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/NanoComarePerformance_03.joined_files.tsv
# python /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/UniversalMethStatsEvaluation.standalone_03.submit.py /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/NanoComarePerformance_03.joined_files.03.tsv

import csv
import os
import subprocess
from sys import argv

from nanocompare.global_config import src_base_dir

pyFile = os.path.join(src_base_dir, 'nanocompare', "read_level_eval.sbatch")

if __name__ == '__main__':

    infile = open(argv[1], 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')
    for row in csvfile:
        if row['status'] == "submit":
            logdir = os.path.join('.', 'log')
            os.makedirs(logdir, exist_ok=True)

            command = \
                f"""
set -x; sbatch --job-name=meth_perf_{row['RunPrefix']} --output=log/%x.%j.out --error=log/%x.%j.err \
 --export=ALL,DeepSignal_calls="{row['DeepSignal_calls']}",Tombo_calls="{row['Tombo_calls']}",Nanopolish_calls="{row['Nanopolish_calls']}",DeepMod_calls="{row['DeepMod_calls']}",Megalodon_calls="{row['Megalodon_calls']}",bgTruth="{row['bgTruth']}",\
RunPrefix="{row['RunPrefix']}",parser="{row['parser']}",minCov="{row['minCov']}",dsname="{row['Dataset']}" {pyFile}
"""

            print(row['RunPrefix'], subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8"))
