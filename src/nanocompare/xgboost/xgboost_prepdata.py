#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : xgboost_prepdata.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome
import argparse
import os.path
from functools import reduce

import pandas as pd

from nanocompare.eval_common import load_tool_read_level_unified_as_df
from nanocompare.global_config import pic_base_dir
from nanocompare.global_settings import NANOME_VERSION

parser = argparse.ArgumentParser(prog='xgboost_prepdata (NANOME)', description='XGBoost preparation of data')
parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
parser.add_argument('--dsname', type=str, required=True,
                    help='Dataset name')
parser.add_argument('--methodsfile', '-i', type=str, required=True,
                    help='TSV file containing name and path of the method output tsv file. The output tsv file from the method should be in the format [ID, Pos, Strand, Score]. Can be compressed in gz.')
parser.add_argument('--bsseq', type=str, required=True,
                    help='BS-seq data of sorted bed file')
parser.add_argument('--contain-na', help="if allow merge with NA values", action='store_true')
parser.add_argument('--verbose', help="if output verbose info", action='store_true')

parser.add_argument('-o', type=str, default=pic_base_dir,
                    help='Where to store the outputs')
args = parser.parse_args()
print(f"args={args}", flush=True)

if __name__ == '__main__':
    df_tsvfile = pd.read_csv(args.methodsfile, header=None, sep='\t')
    tool_list = '_'.join(list(df_tsvfile[0]))
    dflist = []
    for index, row in df_tsvfile.iterrows():
        df = load_tool_read_level_unified_as_df(row[1], row[0])
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
    ## Save combined of prediction by tools
    print(f"combine_df={pred_combine_df}")
    outfn = os.path.join(args.o, f'{args.dsname}_{tool_list}_pred{"_containNA" if args.contain_na else ""}.tsv.gz')
    pred_combine_df.to_csv(outfn, sep='\t', index=False)
    print(f"Save to {outfn}", flush=True)

    ## Load bsseq, note site level input is 0-start, 1-end
    print(f"Load bsseq data:{args.bsseq}")
    bsseq_df = pd.read_csv(args.bsseq, sep='\t', header=None, index_col=False).iloc[:, [0, 2, 5, 6, 7]]
    bsseq_df.columns = ['Chr', 'Pos', 'Strand', 'Freq', 'Coverage']
    bsseq_df.drop_duplicates(subset=["Chr", "Pos", "Strand"], inplace=True)
    print(f"bsseq_df={bsseq_df}", flush=True)

    ## Combine bsseq to pred_df
    if args.contain_na:
        pred_bsseq_combine_df = pred_combine_df.merge(bsseq_df, how='left', on=['Chr', 'Pos', 'Strand'])
    else:
        pred_bsseq_combine_df = pred_combine_df.merge(bsseq_df, how='inner', on=['Chr', 'Pos', 'Strand'])
    print(f"pred_bsseq_combine_df={pred_bsseq_combine_df}")
    outfn = os.path.join(args.o,
                         f'{args.dsname}_{tool_list}_pred_bsseq_combine{"_containNA" if args.contain_na else ""}.tsv.gz')
    pred_bsseq_combine_df.to_csv(outfn, sep='\t', index=False)
    print(f"Save to {outfn}", flush=True)
