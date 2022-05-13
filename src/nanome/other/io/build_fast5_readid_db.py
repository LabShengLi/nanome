#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : get_fast5_readid.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Generate a table for readid and filename.
"""

import os.path

import glob
import pandas as pd
from ont_fast5_api.fast5_interface import get_fast5_file
from tqdm import tqdm

from nanome.common.global_config import set_log_debug_level, logger, pic_base_dir

dsname = "NA19240"
basedir = "/projects/li-lab/Nanopore_compare/nanopore_fast5/NA19240-N300-sept"

if __name__ == '__main__':
    set_log_debug_level()

    fnlist = glob.glob(os.path.join(basedir, '*', '*.fast5'), recursive=False)
    logger.debug(f"number of files: {len(fnlist)}")

    readid_list = []
    fast5_fn_list = []

    for infn in tqdm(fnlist[:]):
        try:
            with get_fast5_file(infn, mode="r") as f5:
                for read in f5.get_reads():
                    # print(read.read_id, raw_data)
                    readid_list.append(read.read_id)
                    fast5_fn_list.append(infn)
        except:  # skip failed fast5 files
            pass

    df = pd.DataFrame({'ReadID': readid_list, 'Fast5FileName': fast5_fn_list})
    logger.debug(f"df={df}")
    outfn = os.path.join(pic_base_dir, f'{dsname}_readid_filename_db.csv.gz')
    df.to_csv(outfn)
    logger.info(f"save to {outfn}")
    logger.info("### find fast5 readid DONE")
