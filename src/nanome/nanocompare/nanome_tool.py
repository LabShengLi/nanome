#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : nanome_tool.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome


import argparse
import gzip

from tqdm import tqdm

from nanome.common.eval_common import open_file_gz_or_txt
from nanome.common.global_config import set_log_debug_level, logger, set_log_info_level
from nanome.common.global_settings import NANOME_VERSION


def parse_arguments():
    parser = argparse.ArgumentParser(prog='nanome_tool (NANOME)',
                                     description='Plot and export data for Nanocompare paper.')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument("cmd", help="name of command: to_ucsc_bed, etc.")
    parser.add_argument('-i', type=str, help='input file name', required=True)
    parser.add_argument('-o', type=str, help="output file name",
                        required=True)
    parser.add_argument('--track-name', type=str, help="track name in BED file",
                        default='.')
    parser.add_argument('--scale', type=float, help="track name in BED file",
                        default=100)
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    parser.add_argument('--score-type', type=str, help="track name in BED file",
                        default='int')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()

    logger.debug(args)

    if args.cmd == 'to_ucsc_bed':
        ## Convert sorted bed.gz (site level) into UCSC BED format
        outfile = gzip.open(args.o, 'wt')
        infile, nlines = open_file_gz_or_txt(args.i)
        for row in tqdm(infile, total=nlines, desc=f"{args.cmd}"):
            tmp = row.strip().split('\t')
            if args.score_type == 'int':
                score = int(float(tmp[6]) * args.scale)
                outfile.write(
                    f"{tmp[0]}\t{tmp[1]}\t{tmp[2]}\t{args.track_name}\t{score}\t{tmp[5]}\n")
            else:
                # float results
                score = float(float(tmp[6]) * args.scale)
                outfile.write(
                    f"{tmp[0]}\t{tmp[1]}\t{tmp[2]}\t{args.track_name}\t{score:.3f}\t{tmp[5]}\n")

        infile.close()
        outfile.close()
    else:
        raise Exception(f'Command={args.cmd} is not support')

    logger.info('Nanome tool DONE')
