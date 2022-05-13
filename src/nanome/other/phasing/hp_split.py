#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : hp_split.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Split megalodon results based on HP types (H1 and H2) of readid list for each chromosome
For megalodon 5hmc_5mc results, also support splitting into:
    1. H1 and H2 for 5mC and 5C
    2. H1_5hmc and H2_5hmc for 5hmc and 5C
"""
import argparse
import os

import pandas as pd

from nanome.common.global_config import set_log_debug_level, logger, pic_base_dir
from nanome.common.global_settings import NANOME_VERSION
from nanome.other.phasing.mega_parser import import_megalodon_per_read_file, agg_read_to_site, to_read_preds_file, \
    to_site_freq_file, \
    import_nanome_per_read_file

hp_set = ['H1', 'H2']


def parse_arguments():
    parser = argparse.ArgumentParser(prog='hp_split (NANOME)', description='Split megalodon results based on HP types.')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', type=str, help='input file for megalodon per-read', required=True)
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--region', nargs='+',
                        help="chr list for region filtered, default is None, means all chr, no filter",
                        default=None)
    parser.add_argument('--haplotype-list', type=str, help="haplotype list file", required=True)
    parser.add_argument('--num-class', type=int,
                        help="number of class for input file, 2 for 5mc/5c, 3 for 5c/5mc/5hmc, default is 2",
                        default=2)
    parser.add_argument('--tool', type=str,
                        help="input methylation file from tool-name, such as megalodon, nanome3t, nanome2t, megalodon_5hmc_5mc etc., default is megalodon",
                        default='megalodon')
    parser.add_argument('--encode', type=str,
                        help="input methylation file format encode: megalodon/megalodon.ziwei, nanome, etc, default is megalodon",
                        default='megalodon')
    parser.add_argument('--only-test', type=int, help="only for test import few lines",
                        default=None)
    parser.add_argument('-o', type=str, help="output dir",
                        default=None)
    parser.add_argument('--save-unified-read', help="if save a unified read, when read raw input", action='store_true')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    set_log_debug_level()
    args = parse_arguments()
    logger.debug(args)

    if args.o is None:
        outdir = os.path.join(pic_base_dir, f'hp_split_{args.dsname}_{args.tool}')
    else:
        outdir = os.path.join(args.o, f'hp_split_{args.dsname}_{args.tool}')
    os.makedirs(outdir, exist_ok=True)
    logger.debug(f"outdir={outdir}")

    if args.haplotype_list is not None:
        hpdf = pd.read_csv(args.haplotype_list, sep='\t', header=0, index_col=None)
        hpdf = hpdf[hpdf['haplotype'].isin(hp_set)]
        logger.debug(hpdf)
        logger.debug(hpdf['haplotype'].value_counts())
    else:
        hpdf = None

    if args.region is None:
        chrFilter = None
        chrTagname = ""
    else:
        chrFilter = set(args.region)
        chrTagname = "_" + '_'.join(chrFilter)

    for hp_type in hpdf['haplotype'].unique():
        logger.debug(f"hp_type={hp_type}")
        readids = hpdf[hpdf['haplotype'] == hp_type]['#readname']
        logger.debug(f"readids={readids}")

        if len(readids) <= 0:
            continue

        logger.debug(f"Load {str(args.encode).upper()}:{args.i}")
        outfn = os.path.join(outdir,
                             f"{args.dsname}_{str(args.tool).lower()}_perReadScore{chrTagname}_{hp_type}.tsv.gz")

        if str(args.encode).lower() == 'megalodon':
            predDict = import_megalodon_per_read_file(args.i, readid_filter=set(readids), chr_filter=chrFilter,
                                                      only_test=args.only_test,
                                                      save_unified_format=args.save_unified_read, outfn=outfn,
                                                      num_class=args.num_class)
        elif str(args.encode).lower() == 'megalodon.ziwei':
            predDict = import_megalodon_per_read_file(args.i, readid_filter=set(readids), chr_filter=chrFilter,
                                                      chr_col=0, readid_col=2, pos_col=1, strand_col=3,
                                                      only_test=args.only_test,
                                                      save_unified_format=args.save_unified_read, outfn=outfn)
        elif str(args.encode).lower() in ['nanome']:
            predDict = import_nanome_per_read_file(args.i, readid_filter=set(readids), chr_filter=chrFilter,
                                                   only_test=args.only_test, save_unified_format=args.save_unified_read,
                                                   outfn=outfn)
        else:
            raise Exception(f"Not support param tool, args.tool={args.tool}")

        siteDict = agg_read_to_site(predDict, num_class=args.num_class)

        outfn = os.path.join(outdir,
                             f'{args.dsname}_{args.tool.lower()}_read_pred{chrTagname}_{hp_type}.tsv.gz')
        to_read_preds_file(predDict, outfn)
        logger.info(f"save to {outfn}")

        outfn = os.path.join(outdir,
                             f'{args.dsname}_{args.tool.lower()}_site_freq{chrTagname}_{hp_type}.tsv.gz')
        to_site_freq_file(siteDict, outfn, num_class=args.num_class)
        logger.info(f"save to {outfn}")
        logger.info(f"### hp_split DONE")
