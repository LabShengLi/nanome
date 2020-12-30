"""
This file will seperate all fast5 file into some subfolders
"""

import os
import shutil
from sys import argv

# return all .fast5 file to a list
from tqdm import tqdm


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
    MainOutputDir = argv[3]  # "/fastscratch/rosikw/APL_newSept"

    fast5files = recursiveFast5(path)  # [1,2,3,4,5,6,7,8,9,0]
    print(f'total files: {len(fast5files)}, indir={path}, TargetNum={TargetNum}, outdir={MainOutputDir}')

    for k in range(TargetNum):
        fdir = os.path.join(MainOutputDir, f'{k}')
        if not os.path.exists(fdir):
            os.umask(0)
            os.makedirs(fdir, exist_ok=True)

    for k, fn in tqdm(enumerate(fast5files)):
        destdir = os.path.join(MainOutputDir, f'{k % TargetNum}')
        shutil.move(fn, destdir)
        #
        # if k % 1000 == 0:
        #     print(f'Processed #{k} files.')


# Samples: python /projects/liuya/workspace/long_read/utils/FilesSeparator_04_deprecated.py /fastscratch/liuya/AB.Nanop/input 1 /fastscratch/liuya/AB.Nanop/sept yang.liu@jax.org
if __name__ == "__main__":
    main()
    # shutil.move('test.txt', "/projects/liuya")
