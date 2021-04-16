"""
Construct the benchmarking data for resource usage evaluation
"""
import glob
import os
import random
from shutil import copy

import numpy as np

from nanocompare.global_config import set_log_debug_level, logger, pic_base_dir

base_input = '/projects/li-lab/Nanopore_compare/nanopore_fast5/NA19240-N300-sept/10'
out_dir = os.path.join(pic_base_dir, 'BenchData8000')

# sizes = [int(x) for x in np.linspace(1000, 10000, 10)]

sizes = [int(x) for x in np.linspace(800, 8000, 10)]


random.seed(688)

if __name__ == '__main__':
    set_log_debug_level()

    fnlist = glob.glob(os.path.join(base_input, '*.fast5'))
    logger.info(len(fnlist))

    logger.info(f'sizes={sizes}')


    os.makedirs(out_dir, exist_ok=True)

    for t in sizes:
        retfnlist = random.choices(fnlist, k=t)

        benchDir = os.path.join(out_dir, f'MB{t / 1000:02n}K')
        os.makedirs(benchDir, exist_ok=True)

        for fn in retfnlist:
            # copy fn into benchDir
            copy(fn, benchDir)
        logger.info(f'Copy done for t={t}, benchDir={benchDir}')

    pass
