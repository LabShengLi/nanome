"""
This tool is to clean old basecall analyses for NA12878, etc.
"""
import argparse

import h5py

# Group name need to be cleaned in fast5 files
groupName = "Analyses"


def parse_arguments():
    parser = argparse.ArgumentParser(description='Clean old basecall info in fast5 files.')
    parser.add_argument('-i', nargs='+', help='list of input fast5 files', default=[])
    # parser.add_argument('-o', type=str, help='output dir', default=None)
    parser.add_argument('--verbose', action='store_true', help="True if print debug info")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    # print(args)
    cntClean = 0
    for fn in args.i:
        with h5py.File(fn, 'r+') as handle:  # ref: https://docs.h5py.org/en/stable/high/file.html?highlight=h5py.File#h5py.File
            # print(f"Detect old {groupName}, we clean it.")
            if groupName in list(handle.keys()):  # clean if found any group named 'Analyses'
                del handle[groupName]
                cntClean += 1
    if args.verbose:
        print(f"Total cleaned files: {cntClean}")
