#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : FilesSeparator.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/liuyangzzu/nanome

"""
Seperate all fast5 file into more subfolders 1--N

Usage:
    python FileSeparator.py   <original-folder>  N  <output-base-folder>
"""

import os
import shutil
from sys import argv

from tqdm import tqdm


def recursive_find_fast5(basedir):
    """
    Find all fast5 files recursively in folder
    """
    fast5files = []
    for root, dirs, files in os.walk(basedir):
        for file in files:
            if file.endswith(".fast5"):
                fast5files.append(os.path.join(root, file))
    return fast5files


def main():
    indir = argv[1]  # "/fastscratch/rosikw/APL/readsAlbacore"
    target_num = int(argv[2])  # 500
    outdir = argv[3]  # "/fastscratch/rosikw/APL_newSept"
    if len(argv) >= 5:  # prefixStr, such as M1, M2, etc.
        prefix_str = argv[4]
    else:
        prefix_str = ""

    fast5files = recursive_find_fast5(indir)  # [1,2,3,4,5,6,7,8,9,0]
    print(f'total files: {len(fast5files)}, indir={indir}, TargetNum={target_num}, outdir={outdir}')

    for k in range(target_num):
        fdir = os.path.join(outdir, f'{prefix_str}{k + 1}')
        if not os.path.exists(fdir):
            os.umask(0)
            os.makedirs(fdir, exist_ok=True)

    for k, fn in tqdm(enumerate(fast5files)):
        destdir = os.path.join(outdir, f'{prefix_str}{(k % target_num) + 1}')
        shutil.move(fn, destdir)


if __name__ == "__main__":
    main()
