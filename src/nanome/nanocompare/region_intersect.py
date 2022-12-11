#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : region_intersect.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome
"""
Tool for CpGs intersect with regions in nanome paper
"""
import argparse
import os

import pandas as pd
import pybedtools

from nanome.common.eval_common import bedtool_convert_0_to_1
from nanome.common.global_config import set_log_debug_level, set_log_info_level, logger
from nanome.common.global_settings import NANOME_VERSION


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='region_intersect (NANOME)',
                                     description='Intersect data with region file')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', type=str, help="input data file", required=True)
    parser.add_argument('-r', type=str, help="region BED file", required=True)
    parser.add_argument('-o', type=str, help="output file", required=True)
    parser.add_argument('--sep', type=str, help="file seperator, default is TAB", default='\t')
    parser.add_argument('--bed-format', type=int, help="BED file start format 0/1, default is 1", default=1)
    parser.add_argument('--strand-intersect', help="if BED file intersect using strand sensitive mode",
                        action='store_true')
    parser.add_argument('--header', type=int, help="if input file contain header, set to 0, default is None",
                        default=None)
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()
    logger.debug(args)

    indf = pd.read_csv(args.i, sep=args.sep, index_col=False, header=args.header)
    logger.debug(f"indf={indf}")

    header_names = indf.columns
    beddf = pybedtools.BedTool.from_dataframe(indf).sort()

    bed_region = pybedtools.BedTool(args.r)
    if args.bed_format == 0:
        ## convert 0 start to 1 start before intersect
        bed_region = bedtool_convert_0_to_1(bed_region)

    if os.path.basename(args.r).startswith('hg38.repetitive') or args.strand_intersect:
        intersectBed = beddf.intersect(bed_region, u=True, wa=True, s=True)
    else:
        intersectBed = beddf.intersect(bed_region, u=True, wa=True)

    logger.debug(f"intersectBed={len(intersectBed):,}")
    outdf = intersectBed.to_dataframe(names=header_names)
    logger.debug(outdf)

    outdf.to_csv(args.o, header=True if args.header is not None else False, index=False, sep=args.sep)
    logger.info(f"save to {args.o}")
    logger.info("### Done")
