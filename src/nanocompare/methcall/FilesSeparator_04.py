from sys import argv
import os
from math import ceil
import subprocess


# return all .fast5 file to a list
def recursiveFast5(path):
    fast5files = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".fast5"):
                fast5files.append(os.path.join(root, file))
    return fast5files


# generate a subgroup of fast5 file set into a bash file
def processFiles(inputList, MainOutputDir, identifier, bashPrefix, ignoreSubfolders):
    outfile = open("{}_batch_{}.sh".format(bashPrefix, identifier), 'w')
    outfile.write("mkdir -p {}/{}\n".format(MainOutputDir, identifier))
    for elem in inputList:
        if ignoreSubfolders not in elem:
            outfile.write("cp {} {}/{}/{}\n".format(elem, MainOutputDir, identifier, os.path.basename(elem)))
        else:
            print("ignored:", elem)
    outfile.close()
    print("FileSeperator: ### processFiles: batch {} completed!".format(identifier))


def main():
    path = argv[1]  # "/fastscratch/rosikw/APL/readsAlbacore"
    TargetNum = int(argv[2])  # 500
    bashPrefix = "FilesSeparator"
    fast5files = recursiveFast5(path)  # [1,2,3,4,5,6,7,8,9,0]
    ignoreSubfolders = "notUsed__dffdsddsdfsf"  # "fail"
    MainOutputDir = argv[3]  # "/fastscratch/rosikw/APL_newSept"
    email = argv[4] if len(argv) >= 5 else 'yang.liu@jax.org'

    command = "mkdir -p {}".format(MainOutputDir)
    # print(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read())
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()

    i = 0
    identifier = 0
    tmp = []
    # partialLen = ceil(len(fast5files)/TargetNum)
    step = len(fast5files) // TargetNum
    print("FileSeperator: Each of {} folders will contain {} out of {} files.".format(TargetNum, step, len(fast5files)))
    for k in range(TargetNum - 1):
        tmp = fast5files[k * step:(k + 1) * step]
        # Generate N-1 subgroup copy bash scripts
        processFiles(tmp, MainOutputDir, k, bashPrefix, ignoreSubfolders)

    # last group
    tmp = fast5files[(TargetNum - 1) * step:]
    # Generate Nth subgroup copy bash script
    processFiles(tmp, MainOutputDir, TargetNum - 1, bashPrefix, ignoreSubfolders)

    #	for f in fast5files:
    #		if i < partialLen:
    #			tmp.append(f)
    #			i += 1
    #		else:
    #			processFiles(tmp, MainOutputDir, identifier, bashPrefix, ignoreSubfolders)
    #	#         print(tmp, identifier)
    #
    #			tmp = []
    #			tmp.append(f)
    #			i = 1
    #			identifier += 1
    #	processFiles(tmp, MainOutputDir, identifier, bashPrefix, ignoreSubfolders)
    # print(tmp, identifier)

    # Generate master separator bash
    outfile = open("{}_MasterArray.sh".format(bashPrefix, identifier), 'w')
    outfile.write("""#!/bin/bash
#SBATCH --job-name=sept
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 300g # memory pool for all cores
#SBATCH -t 5-14:00 # time (D-HH:MM)
#SBATCH -o log/%x.%j.out # STDOUT
#SBATCH -e log/%x.%j.err # STDERR
#SBATCH --array=1-{}

########################
###PBS -q batch # batch - up to 72h, long - up to 504h. 
##PBS -l nodes=1:ppn=1   # if I want 32 cores, using machines that has 16 cores, I should specify nodes=2:ppn=16 and be sure 
##PBS -l mem=1000mb # this will ask for 1000mb of RAM; 300gb
##PBS -l walltime=36:00:00
##PBS -M 
##PBS -m abe # be emailed at the beginning, end and error
##PBS -o /home/rosikw/qsubReports
##PBS -e /home/rosikw/qsubReports
###PBS -t 1- # array 

set -x

echo $SLURM_ARRAY_TASK_ID
bash FilesSeparator_batch_$((SLURM_ARRAY_TASK_ID-1)).sh

	""".format(TargetNum, bashPrefix))

    outfile.close()


# Samples: python /projects/liuya/workspace/long_read/utils/FilesSeparator_04.py /fastscratch/liuya/AB.Nanop/input 1 /fastscratch/liuya/AB.Nanop/sept yang.liu@jax.org
if __name__ == "__main__":
    main()
