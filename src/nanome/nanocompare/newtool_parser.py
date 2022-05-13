#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : newtool_parser.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Parse new tool's results
"""
import argparse

from nanome.common.eval_common import *
from nanome.common.global_settings import NANOME_VERSION


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='newtool_parser (NANOME)',
                                     description='Parse new comming tool outputs')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', type=str, help="input file name", required=True)
    parser.add_argument('--read-out', type=str, help="read level unified output file name", required=True)
    parser.add_argument('--site-out', type=str, help="read level unified output file name", required=True)
    parser.add_argument('--column-order', nargs='+', type=int,
                        help='the columns of READID, CHR, POS, STRAND and SCORE',
                        default=[0, 1, 2, 3, 4])
    parser.add_argument('--score-cols', nargs='+', type=int,
                        help='the SCORE column(s), if len=2, [0] is 5mC, [1] is 5C',
                        default=[5])
    parser.add_argument('--log-score', help="if input is log score", action='store_true')
    parser.add_argument('--baseFormat', type=str, help="POS base: 0/1, default is 0", default=0)
    parser.add_argument('--sep', type=str, help="seperator for output csv file, default is tab character", default='\t')
    parser.add_argument('--chrSet', nargs='+', help='chromosome list, default None for no filter',
                        default=None)
    parser.add_argument('--sort', help="if sort output", action='store_true')
    parser.add_argument('--deduplicate', help="if deduplicate read level output", action='store_true')
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    logger.debug(f"args={args}")
    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()

    readid_col = args.column_order[0]
    chr_col = args.column_order[1]
    pos_col = args.column_order[2]
    strand_col = args.column_order[3]

    cpgDict = {}

    infile, lines = open_file_gz_or_txt(args.i)
    read_outf = gzip.open(args.read_out, 'wt')
    read_outf.write(f"ID\tChr\tPos\tStrand\tScore\n")

    for row in tqdm(infile, total=lines):
        tmp = row.strip().split(args.sep)

        if args.chrSet is not None and tmp[chr_col] not in args.chrSet:
            continue
        try:
            pos = int(tmp[pos_col]) + (1 - args.baseFormat)
            if not args.log_score:
                if len(args.score_cols) != 1:
                    raise Exception(f"not support score_cols={args.score_cols}")
                score_col = args.score_cols[0]
                log_score = prob_to_log(float(tmp[score_col]))
            else:
                if len(args.score_cols) == 2:
                    log_score = float(tmp[args.score_cols[0]]) - float(tmp[args.score_cols[1]])
                elif len(args.score_cols) == 1:
                    log_score = float(tmp[args.score_cols[0]])
                else:
                    raise Exception(f"not support score_cols={args.score_cols}")
            read_outf.write(f"{tmp[readid_col]}\t{tmp[chr_col]}\t{pos}\t{tmp[strand_col]}\t{log_score}\n")

            key = (tmp[chr_col], pos, tmp[strand_col])
            if key not in cpgDict:
                cpgDict[key] = [0, 0]
            if log_score > 0:
                cpgDict[key][1] = cpgDict[key][1] + 1
            else:
                cpgDict[key][0] = cpgDict[key][0] + 1
        except Exception as err:
            logger.error(f"Parse error for: tmp={tmp}, err={err}")
    infile.close()
    read_outf.close()

    if args.sort or args.deduplicate:
        logger.debug(f"Sort or deduplicat read level")
        sort_read_fn = args.read_out + '.sort.gz'
        sort_per_read_tsv_file(args.read_out, sort_read_fn, deduplicate=args.deduplicate)
        os.remove(args.read_out)
        os.rename(sort_read_fn, args.read_out)

        ## Reparse after deduplicate and sort
        chr_col = 1
        pos_col = 2
        strand_col = 3
        score_col = 4
        cpgDict = {}
        with gzip.open(args.read_out, 'rt') as infile:
            for row in infile:
                if row.startswith("ID\tChr"):  # skim header
                    continue
                tmp = row.strip().split(args.sep)

                if args.chrSet is not None and tmp[chr_col] not in args.chrSet:
                    continue
                try:
                    pos = int(tmp[pos_col])
                    score = float(tmp[score_col])

                    key = (tmp[chr_col], pos, tmp[strand_col])
                    if key not in cpgDict:
                        cpgDict[key] = [0, 0]
                    if score > 0:
                        cpgDict[key][1] = cpgDict[key][1] + 1
                    else:
                        cpgDict[key][0] = cpgDict[key][0] + 1
                except Exception as err:
                    logger.error(f"Parse error for: tmp={tmp}, err={err}")

    logger.debug(f"Site level output: cpgDict={len(cpgDict)}")
    with gzip.open(args.site_out, 'wt') as outf:
        for key in cpgDict:
            coverage = cpgDict[key][0] + cpgDict[key][1]
            freq = cpgDict[key][1] / coverage
            strlist = [key[0], str(key[1] - 1), str(key[1]), '.', '.', key[2], str(freq),
                       str(coverage)]
            outf.write('\t'.join(strlist) + '\n')

    if args.sort or args.deduplicate:
        logger.debug(f"Sort or deduplicat site level")
        sort_site_fn = args.site_out + '.sort.gz'
        sort_bed_file(args.site_out, sort_site_fn, deduplicate=args.deduplicate)
        os.remove(args.site_out)
        os.rename(sort_site_fn, args.site_out)

    logger.info("### newtool parser DONE")
