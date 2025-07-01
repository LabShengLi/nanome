#!/usr/bin/env python3
"""
Convert nanome site level output to methylkit format

Nanome site level input:

zcat GT23-09911_Guppy-perSite-cov1.sort.bed.gz| head
chr1	2714	2715	.	.	+	1.0	1
chr1	2715	2716	.	.	-	0.25	4
chr1	2722	2723	.	.	+	1.0	1

Output is base 1-base, used for methylkit's input data format

zcat GM24385_guppy_NanoMethPhase_HP1_methylkit_format.CpG.txt.gz| head
chrBase	chr	base	strand	coverage	freqC	freqT
chr1.2715	chr1	2715	F	1	100.0	0.0
chr1.2723	chr1	2723	F	1	100.0	0.0

Note: methylkit prefer 1-based input, since the output DMC will be correct, default input tsv is no header
"""

import argparse

import pandas as pd
import logging

# interested chrs
chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
        'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX', 'chrY']


def parse_arguments():
    parser = argparse.ArgumentParser(prog='Convert nanome site level format to methylkit format')
    parser.add_argument('-i', type=str, required=True,
                        help='input file')
    parser.add_argument('-o', type=str, required=True,
                        help='output methylkit formated file')
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    parser.add_argument('--chrSet', nargs='+', help='chromosome list, default is human chromosome chr1-22, X and Y',
                        default=chrs)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    else:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    logging.debug(f"args={args}")

    df = pd.read_csv(args.i, sep='\t', index_col=None, header=None)
    df.columns = ['chromosome', 'pos0', 'pos1',
                  'place_hold1', 'place_hold2', 'strand',
                  'meth_freq', 'coverage']
    df = df[df.chromosome.isin(args.chrSet)]

    logging.debug(f"df={df}")

    df['chrBase'] = df['chromosome'] + '.' + df['pos1'].astype(str)
    df['chr'] = df['chromosome']
    df['base'] = df['pos1']  # output base column is 1-based
    df['strand'] = df['strand'].apply(lambda x: 'F' if x == '+' else 'R')
    df = df[df['coverage'] > 0]

    df['freqC'] = df['meth_freq'] * 100.0
    df['freqT'] = 100.0 - df['freqC']

    df = df[['chrBase', 'chr', 'base', 'strand', 'coverage', 'freqC', 'freqT']]
    logging.debug(f"new df={df}")

    df.to_csv(args.o, sep='\t', index=False)
    logging.info(f"save to {args.o}")
