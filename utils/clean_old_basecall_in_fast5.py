#!/usr/bin/env python3

"""
This tool is to clean old basecall analyses for NA12878, etc.
"""
import argparse
import glob
import os.path
from multiprocessing import Pool

import numpy as np
import h5py

# Group name need to be cleaned in fast5 files
groupName = "Analyses"


def parse_arguments():
    parser = argparse.ArgumentParser(description='Clean old basecall info in fast5 files.')
    parser.add_argument('-i', nargs='+', help='list of input fast5 files', default=[])
    parser.add_argument('--processor', type=int, help='number of processors/threads', default=8)
    parser.add_argument('--verbose', action='store_true', help="True if print debug info")
    parser.add_argument('--is-indir', action='store_true', help="True if input is folder")

    args = parser.parse_args()
    return args


def clean_filelist(fnlist):
    """
    Clean for a filelist, will write to them
    :param fnlist:
    :return:
    """
    cntClean = 0
    for fn in fnlist:
        try:
            with h5py.File(fn, 'r+') as handle:  # ref: https://docs.h5py.org/en/stable/high/file.html?highlight=h5py.File#h5py.File
                # print(f"Detect old {groupName}, we clean it.")
                if groupName in list(handle.keys()):  # clean if found any group named 'Analyses'
                    del handle[groupName]
                    cntClean += 1
        except: ## avoid corrupted fast5 files
            pass
    return cntClean


if __name__ == '__main__':
    args = parse_arguments()
    # print(args)
    cntClean = 0

    if args.is_indir:
        fnlist = []
        for bdir in args.i:
            fnlist += glob.glob(os.path.join(bdir, "*.fast5"))
        pass
    else:
        fnlist = args.i

    if args.verbose or args.is_indir:
        print(f"Total fast5 files: {len(fnlist)}")

    # Split to multiple processing
    fnlistArgs = np.array_split(fnlist, args.processor)
    with Pool(processes=args.processor) as pool:
        ret = pool.map(clean_filelist, fnlistArgs)
    totalCleanedFiles = sum(ret)

    if args.verbose or args.is_indir:
        print(f"Total cleaned files: {totalCleanedFiles}")
