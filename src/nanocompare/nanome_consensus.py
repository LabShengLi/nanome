#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : nanome_consensus.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

"""
Consensus top tools's site level results
"""
import argparse
import os
from functools import reduce

import pandas as pd

from nanocompare.eval_common import sort_bed_file
from nanocompare.global_config import logger, set_log_debug_level, set_log_info_level
from nanocompare.global_settings import nanome_version


def consensus(df):
    """
    Consensus all tools site level results:
        cov -  max of tools cov
        freq - average of tools freq
    Args:
        df:

    Returns:

    """
    cov_df = df.filter(regex=("Coverage_*"))
    max_cov_df = cov_df.max(axis=1).rename('Coverage_NANOME').astype(int)

    freq_df = df.filter(regex=("Freq_*"))
    logger.debug(freq_df)
    if args.method == 'average':
        consensus_df = freq_df.mean(axis=1).rename('Freq_NANOME')
        logger.debug(consensus_df)
    else:
        raise Exception(f"Not support consensus method={args.method}")
    outdf = pd.concat([df.iloc[:, [0, 1, 2, 3]], consensus_df, max_cov_df], axis=1)
    outdf['c3'] = '.'
    outdf['c4'] = '.'
    outdf = outdf[['Chr', 'Start', 'End', 'c3', 'c4', 'Strand', 'Freq_NANOME', 'Coverage_NANOME']]
    logger.debug(outdf)
    return outdf


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='consensus (NANOME)',
                                     description='Consensus site level methylation results of top performers')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{nanome_version}')
    parser.add_argument('--site-reports', nargs='+', help='all site level reports need to be consensus', required=True)
    parser.add_argument('--tools', nargs='+', help="tools' names list", default=None)
    parser.add_argument('--join',
                        help="if join all site level results",
                        action='store_true')
    parser.add_argument('--union',
                        help="if union all site level results",
                        action='store_true')
    parser.add_argument('--method',
                        help="consensus method, default is average",
                        default='average')
    parser.add_argument('-o', type=str, help=f"output file name, output format is gzipped file", required=True)
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()
    logger.debug(f"args={args}")

    df_list = []
    for k, fn in enumerate(args.site_reports):
        if args.tools:
            tagname = args.tools[k]
        else:
            tagname = str(k)
        df = pd.read_csv(fn, sep='\t', header=None).iloc[:, [0, 1, 2, 5, 6, 7]]
        df.columns = ['Chr', 'Start', 'End', 'Strand', f'Freq_{tagname}', f'Coverage_{tagname}']
        df_list.append(df)

    if args.join:
        comb_df = reduce(lambda df1, df2: pd.merge(df1, df2, on=['Chr', 'Start', 'End', 'Strand']), df_list)
        logger.debug(comb_df)
        logger.debug(comb_df.columns)
    elif args.union:
        comb_df = reduce(lambda df1, df2: pd.merge(df1, df2, on=['Chr', 'Start', 'End', 'Strand'], how='outer'),
                         df_list)
        logger.debug(comb_df)
        logger.debug(comb_df.columns)
    else:
        raise Exception(f"No actions yet, please use --join or --union")

    outdf = consensus(comb_df)
    ## Save to a gzip output file
    outdf.to_csv(f"{args.o}.tmp.gz", sep='\t', index=False, header=False)
    ## Sort into output file
    sort_bed_file(f"{args.o}.tmp.gz", args.o)
    os.remove(f"{args.o}.tmp.gz")
    logger.info(f"Save to {args.o}")
    logger.info("NANOME consensus top performers DONE")
