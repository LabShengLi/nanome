#!/usr/bin/env python3

"""
Seperate all fast5 file into many subfolders 1--N
"""

import os
import shutil
from sys import argv
from tqdm import tqdm


def recursiveFast5(path):
    fast5files = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".fast5"):
                fast5files.append(os.path.join(root, file))
    return fast5files


def main():
    path = argv[1]  # "/fastscratch/rosikw/APL/readsAlbacore"
    TargetNum = int(argv[2])  # 500
    MainOutputDir = argv[3]  # "/fastscratch/rosikw/APL_newSept"
    if len(argv) >= 5:  # prefixStr, such as M1, M2, etc.
        prefixStr = argv[4]
    else:
        prefixStr = ""

    fast5files = recursiveFast5(path)  # [1,2,3,4,5,6,7,8,9,0]
    print(f'total files: {len(fast5files)}, indir={path}, TargetNum={TargetNum}, outdir={MainOutputDir}')

    for k in range(TargetNum):
        fdir = os.path.join(MainOutputDir, f'{prefixStr}{k + 1}')
        if not os.path.exists(fdir):
            os.umask(0)
            os.makedirs(fdir, exist_ok=True)

    for k, fn in tqdm(enumerate(fast5files)):
        destdir = os.path.join(MainOutputDir, f'{prefixStr}{(k % TargetNum) + 1}')
        shutil.move(fn, destdir)


if __name__ == "__main__":
    main()
