import glob
import os
from collections import defaultdict

import pandas as pd

from nanome.common.global_config import pic_base_dir

baseDir = '/projects/li-lab/Nanopore_compare/suppdata/basecalls.qc'

flist = glob.glob(os.path.join(baseDir, '*basecall.sequencing_summary.txt'))

cutoff = 10000
print(f"We define long reads with length > {cutoff:,}")

dataset = defaultdict(list)

for fn in flist:
    basename = os.path.basename(fn)
    dsname = basename[0:basename.find('-')]
    df = pd.read_csv(fn, sep='\t')
    total_longreads = sum(df['sequence_length_template'] > cutoff)
    portion_longreads = total_longreads / len(df)
    print(
        f"Dataset={dsname}, total_reads={len(df):,}, total_longreads={total_longreads:,}, portion_longreads={portion_longreads * 100:.2f}%")
    dataset['Dataset'].append(dsname)
    dataset['total_reads'].append(len(df))
    dataset['total_longreads'].append(total_longreads)
    dataset['portion_longreads'].append(portion_longreads)

outdf = pd.DataFrame.from_dict(dataset)
outfn = os.path.join(pic_base_dir, f'datasets_long_reads_cutoff_{cutoff}_table.csv')
outdf.to_csv(outfn)
print(f"save to {outfn}")

"""
dsname=NA19240, total_reads=6,219,931, total_longreads=2,286,093, portion_longreads=36.75%
dsname=NA12878, total_reads=14,791,249, total_longreads=4,776,148, portion_longreads=32.29%
dsname=APL, total_reads=1,154,349, total_longreads=204,778, portion_longreads=17.74%
dsname=K562, total_reads=263,019, total_longreads=82,727, portion_longreads=31.45%
"""
