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
import subprocess
from sys import argv

if __name__ == '__main__':
    baseDir = "/projects/liuya/workspace/tcgajax/nanocompare/meth_stats"
    scriptFn = "Universal_meth_stats_evaluation.sbatch"

    infile = open(argv[1], 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')  # specify the delimiter here! I guess that by default it is space
    for row in csvfile:
        if row['status'] == "submit":
            print(f"read row={row}")
            for tool in ["DeepSignal", "DeepMod", "Nanopolish", "Tombo"]:
                command = """sbatch --export=DeepSignal_calls="{}",Tombo_calls="{}",Nanopolish_calls="{}",DeepMod_calls="{}",bgTruth="{}",RunPrefix="{}/{}",parser="{}",minCov="{}",tool="{}" {}/{}""" \
                    .format(row['DeepSignal_calls'], row['Tombo_calls'], row['Nanopolish_calls'],
                            row['DeepMod_calls'], row['bgTruth'], row['RunPrefix'], row['RunPrefix'],
                            row['parser'], row['minCov'], tool, baseDir, scriptFn)
                print(f"command={command}")
                print(row['RunPrefix'], tool, subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read())
