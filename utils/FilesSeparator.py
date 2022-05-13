#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : FilesSeparator.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Seperate all fast5 file into more subfolders 1--N,
or a flexible number of subfolders that each with N files

Usage:
    sbatch FilesSeparator.sh \
        -i /fastscratch/liuya/nanocompare/hl60_input/HL60-Nanopore_GT18-07373.fast5.tar \
        -o /fastscratch/liuya/nanocompare/hl60_output_test \
        -n 10 --tar-file

    bash FilesSeparator.sh \
       -i ../test_data/demo1.fast5.reads.tar.gz  \
       -o /fastscratch/liuya/nanocompare/demo1_output \
       -n 50  --tar-file --fix-folder-size

"""

import argparse
import os
import shutil
import subprocess
from glob import glob
from multiprocessing import Pool

import math
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
    target_num = args.n
    outdir = os.path.join(args.o, 'sept_dir')
    os.makedirs(outdir, exist_ok=True)
    prefix_str = args.prefix_name
    fast5files = recursive_find_fast5(indir)  # [1,2,3,4,5,6,7,8,9,0]
    print(f'total files: {len(fast5files)}, indir={indir}, TargetNum={target_num}, outdir={outdir}', flush=True)

    if args.fix_folder_size:
        num_folders = math.ceil(len(fast5files) / args.n)
        print(f"Fixed folder size={args.n}, total folders={num_folders}", flush=True)

        for k in range(num_folders):
            destdir = os.path.join(outdir, f'{prefix_str}{(k % target_num) + 1}')
            os.makedirs(destdir, exist_ok=True)
            for fn in fast5files[k * args.n: (k + 1) * args.n]:
                if args.copy:
                    shutil.copy(fn, destdir)
                else:
                    shutil.move(fn, destdir)
    else:
        print(f"Fixed number of folders={args.n}", flush=True)
        for k in range(target_num):
            fdir = os.path.join(outdir, f'{prefix_str}{k + 1}')
            if not os.path.exists(fdir):
                os.umask(0)
                os.makedirs(fdir, exist_ok=True)

        print(f"Start seperate into N={target_num}", flush=True)
        for k, fn in tqdm(enumerate(fast5files)):
            destdir = os.path.join(outdir, f'{prefix_str}{(k % target_num) + 1}')
            if args.copy:
                shutil.copy(fn, destdir)
            else:
                shutil.move(fn, destdir)


def untar_file(fn, outdir):
    """
    Untar fn into outdir
    Args:
        fn:
        outdir:

    Returns:

    """
    if fn.endswith('.tar'):
        commandStr = f"tar -xf {fn} -C {outdir}"
    elif fn.endswith('.tar.gz'):
        commandStr = f"tar -xzf {fn} -C {outdir}"
    else:
        raise Exception(f"{fn} type is not support")
    process = subprocess.Popen(commandStr,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    process.communicate()
    print(f"### Untar finished for {fn}", flush=True)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='FileSeperator (NANOME)', description='File seperator for fast5/tar/tar.gz')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v1.0')
    parser.add_argument('-i', type=str, help="input dir/file", required=True)
    parser.add_argument('-o', type=str, help="output dir", required=True)
    parser.add_argument('-n', type=int, help="number of seperated folders", required=True)
    parser.add_argument('-p', type=int, help="number of processors", default=8)
    parser.add_argument('--prefix-name', type=str, help="subfolder prefix string", default="")
    parser.add_argument('--copy', action='store_true', help="true if using copy, else using move")
    parser.add_argument('--tar-file', action='store_true',
                        help="true if deal with tar/tar.gz files, else are fast5 files")
    parser.add_argument('--fix-folder-size', action='store_true',
                        help="true if fix each subfolder contain at most N files, else fix number of folders as N")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    print(f"args={args}")

    if args.tar_file:  # tar files
        flist = []
        if os.path.isdir(args.i):
            flist += glob(os.path.join(args.i, '**', '*.tar'), recursive=True)
            flist += glob(os.path.join(args.i, '**', '*.tar.gz'), recursive=True)
        elif os.path.isfile(args.i):
            flist += [args.i]
        print(f"### Find tar file list={flist}", flush=True)

        untar_dir = os.path.join(args.o, 'untar_dir')
        os.makedirs(untar_dir, exist_ok=True)

        arg_list = []
        for fn in flist:
            # untar_file(fn, untar_dir)
            arg_list.append((fn, untar_dir,))
        with Pool(args.p) as p:
            p.starmap(untar_file, arg_list)
        indir = untar_dir
    else:  # fast5 files
        indir = args.i

    main()
    print("### FilesSeparator DONE")
