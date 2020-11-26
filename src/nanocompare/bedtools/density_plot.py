#!/home/liuya/anaconda3/envs/nmf/bin/python


import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

nanocompare_prj = "/projects/li-lab/yang/workspace/nano-compare/src"
sys.path.append(nanocompare_prj)

from global_config import set_log_debug_level, pic_base_dir

colnames = ['chr-1', 'start-1', 'end-1', 'meth-freq-1', 'cov-1'] + \
           ['chr-2', 'start-2', 'end-2', 'meth-freq-2', 'cov-2'] + ['nearest-dist']

dtypes = [str, np.int64, np.int64, np.float64, np.int64] + \
         [str, np.int64, np.int64, np.float64, np.int64] + [np.int64]

dtype_dict = {col: dt for col, dt in zip(colnames, dtypes)}


def load_df(infn='K562-WGBS-BGTruth-Tombo-closest.bed'):
    infn = os.path.join('/projects/li-lab/yang/workspace/nano-compare/src/nanocompare/bedtools', infn)
    df = pd.read_csv(infn, sep='\t', header=None, dtype=dtype_dict)
    df.columns = colnames
    df['meth-freq-2'] = df['meth-freq-2'].replace('.', '0.0').astype(np.float64)
    logging.debug(df)
    logging.info(df.dtypes)
    return df


def plot_hist_of_df():
    infn = 'K562-WGBS-Tombo-BGTruth-closest.bed'

    df = load_df(infn=infn)

    data = df['nearest-dist']
    data = data[data <= 10]

    plt.figure(figsize=(5, 4))
    # Density Plot and Histogram of all arrival delays
    ax = sns.distplot(data, hist=True, kde=False, norm_hist=True,
                      bins=10, color='darkblue',
                      hist_kws={'edgecolor': 'black'},
                      kde_kws={'linewidth': 2})

    ax.set_xlim(0, 10)
    ax.set_title('Tombo with BG-Truth on K562')

    # sns.displot(df, x="nearest-dist")
    outfn = os.path.join(pic_base_dir, f'dist-{os.path.basename(infn)}.png')
    plt.savefig(outfn, dpi=600, bbox_inches='tight', format="png")
    logging.info(f"save to {outfn}")

    plt.show()
    plt.close()


if __name__ == '__main__':
    set_log_debug_level()
    plot_hist_of_df()
