#!/usr/bin/env python3
"""
Input is start 0-base, end 1-base, from Nanomethphase phased results meth freq tsv file (with or without header)
Output is base 1-base, used for methylkit input

Note: methylkit prefer 1-based input, since the output DMC will be correct, default input tsv is no header
"""

import argparse

import pandas as pd

from nanome.common.global_config import logger, set_log_debug_level, set_log_info_level


def parse_arguments():
    parser = argparse.ArgumentParser(prog='Convert nanomethphase MethylFreq to methylkit format')

    parser.add_argument('-i', type=str, required=True,
                        help='input file')
    parser.add_argument('-o', type=str, required=True,
                        help='output methylkit formated file')
    parser.add_argument('--has_header', help="if input contains header", action='store_true')
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()
    logger.debug(f"args={args}")

    df = pd.read_csv(args.i, sep='\t', index_col=None, header=True if args.has_header else None)
    df.columns = ['chromosome', 'start', 'end', 'strand', 'NumOfAllCalls', 'NumOfModCalls', 'MethylFreq']
    logger.debug(f"df={df}")
    df['chrBase'] = df['chromosome'] + '.' + df['end'].astype(str)
    df['chr'] = df['chromosome']
    df['base'] = df['end']  # output base column is 1-based
    df['strand'] = df['strand'].apply(lambda x: 'F' if x == '+' else 'R')

    df['coverage'] = df['NumOfAllCalls']
    df['freqC'] = df['MethylFreq'] * 100.0
    df['freqT'] = (1 - df['MethylFreq']) * 100.0

    df = df[['chrBase', 'chr', 'base', 'strand', 'coverage', 'freqC', 'freqT']]
    logger.debug(f"new df={df}")

    df.to_csv(args.o, sep='\t', index=False)
    logger.info(f"save to {args.o}")
