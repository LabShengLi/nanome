#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : clean_old_basecall_in_fast5.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
This tool is to clean old basecall analyses for downloaded nanopore reads.
"""
import argparse
import glob
import os.path
from multiprocessing import Pool

import h5py
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(prog='CleanAnalysis (NANOME)',
                                     description='clean old basecall info in downloaded fast5 files')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v1.0')
    parser.add_argument('-i', nargs='+', help='list of input fast5 files/directories', required=True)
    parser.add_argument('--processor', type=int, help='number of processors/threads, default is set to 8', default=8)
    parser.add_argument('--verbose', action='store_true', help="if print debug/verbose info")
    parser.add_argument('--is-indir', action='store_true', help="if input is folder")
    parser.add_argument('--groupName', type=str,
                        help="the group name in h5 data need to be cleaned, default is 'Analyses'", default="Analyses")
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
            with h5py.File(fn,
                           'r+') as handle:  # ref: https://docs.h5py.org/en/stable/high/file.html?highlight=h5py.File#h5py.File
                if args.groupName in list(handle.keys()):  # clean if found any group named 'Analyses'
                    del handle[args.groupName]
                    cntClean += 1
        except:  ## avoid corrupted fast5 files
            pass
    return cntClean


if __name__ == '__main__':
    args = parse_arguments()

    if args.is_indir:
        fnlist = []
        for bdir in args.i:
            ## No recursive search here, assume all files are flat into the folder
            fnlist += glob.glob(os.path.join(bdir, "*.fast5"))
        pass
    else:
        fnlist = args.i

    if args.verbose:
        print(f"Total fast5 files: {len(fnlist)}")

    # Split to multiple processing
    fnlistArgs = np.array_split(fnlist, args.processor)
    with Pool(processes=args.processor) as pool:
        ret = pool.map(clean_filelist, fnlistArgs)
    total_cleaned = sum(ret)

    if args.verbose:
        print(f"Total cleaned files: {total_cleaned}")
