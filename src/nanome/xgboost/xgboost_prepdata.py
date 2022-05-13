#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : xgboost_prepdata.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome
"""
Prepare data for training, join all tools log-likelyhood (log meth_prob/unmeth_prob) score, with bs-seq
Note:   no-coverage filter for tools and bs-seq, guppy score is from gcf52ref for reference purpose, not certified by us
        contain NA value for each columns
        run faster for splitting chromsomes
"""
import argparse
from functools import reduce

import pandas as pd

from nanome.common.eval_common import load_tool_read_level_unified_as_df
from nanome.common.global_settings import NANOME_VERSION

parser = argparse.ArgumentParser(prog='xgboost_prepdata (NANOME)', description='XGBoost preparation of data')
parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
parser.add_argument('--dsname', type=str, required=True,
                    help='Dataset name')
parser.add_argument('--methodsfile', '-i', type=str, required=True,
                    help='TSV file containing name and path of the method output tsv file. The output tsv file from the method should be in the format [ID, Pos, Strand, Score]. Can be compressed in gz.')
parser.add_argument('--bsseq', type=str, required=True,
                    help='BS-seq data of sorted bed file')
parser.add_argument('-o', type=str, required=True,
                    help='output combine file name')
parser.add_argument('--chrs', nargs='+', help='chromosome list used, default is None, use all chrs', default=None)
parser.add_argument('--contain-na', help="if allow merge with NA values", action='store_true')
parser.add_argument('--verbose', help="if output verbose info", action='store_true')

args = parser.parse_args()
print(f"args={args}", flush=True)

if __name__ == '__main__':
    df_tsvfile = pd.read_csv(args.methodsfile, header=None, sep='\t')
    tool_list = '_'.join(list(df_tsvfile[0]))
    dflist = []
    for index, row in df_tsvfile.iterrows():
        df = load_tool_read_level_unified_as_df(row[1], row[0], filterChrSet=args.chrs)
        dflist.append(df)
    if args.contain_na:
        pred_combine_df = reduce(
            lambda left, right: pd.merge(left, right, how='outer', on=["ID", "Chr", "Pos", "Strand"]),
            dflist)
    else:
        pred_combine_df = reduce(
            lambda left, right: pd.merge(left, right, how='inner', on=["ID", "Chr", "Pos", "Strand"]),
            dflist)
    ## Del dflist for save memory
    dflist = None
    pred_combine_df.drop_duplicates(subset=["ID", "Chr", "Pos", "Strand"], inplace=True)
    pred_combine_df.reset_index(inplace=True, drop=True)

    if len(pred_combine_df) <= 0:
        raise Exception(f"The combined results are empty, for tool_list={tool_list}, fn_list={df_tsvfile[1]}")

    ## Load bsseq, note site level input is 0-start, 1-end
    print(f"Load bsseq data:{args.bsseq}")
    bsseq_df = pd.read_csv(args.bsseq, sep='\t', header=None, index_col=False).iloc[:, [0, 2, 5, 6, 7]]
    bsseq_df.columns = ['Chr', 'Pos', 'Strand', 'Freq', 'Coverage']
    bsseq_df.drop_duplicates(subset=["Chr", "Pos", "Strand"], inplace=True)
    if args.chrs is not None:
        bsseq_df = bsseq_df[bsseq_df['Chr'].isin(args.chrs)]

    print(f"bsseq_df={bsseq_df}", flush=True)

    ## Combine bsseq to pred_df
    if args.contain_na:
        pred_bsseq_combine_df = pred_combine_df.merge(bsseq_df, how='outer', on=['Chr', 'Pos', 'Strand'])
    else:
        pred_bsseq_combine_df = pred_combine_df.merge(bsseq_df, how='inner', on=['Chr', 'Pos', 'Strand'])
    print(f"pred_bsseq_combine_df={pred_bsseq_combine_df}")
    pred_bsseq_combine_df.to_csv(args.o, index=False)
    print(f"Save to {args.o}", flush=True)

    print(f"### xgboost_prepdata Done")
