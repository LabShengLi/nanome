#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : methcall2bed.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Converts any tool's perReadScore methylation results file to a bed
format and also splits multi-group CpG sites to single group by read-id.

Input is 1-based per-read score (perReadScore unified) file

Output format (methcall2bed):
chrY    10624514        10633058        +       52a2067f-7588-4190-8c15-85887a0f671b    meth_log_ratio  unmeth_log_ratio  meth_sites  unmeth_sites

Note:   start is 0-based, end is 1 based
        for + strand, sites point to CG's C
        for - strand, sites point to CG's C also. (Later will be aligned to CG's G by bam2bis function)
"""
import argparse
import gzip
from collections import defaultdict

from tqdm import tqdm

from nanome.common.eval_common import open_file_gz_or_txt
from nanome.common.global_config import set_log_debug_level, logger, set_log_info_level
from nanome.common.global_settings import NANOME_VERSION


def process_read_score_by_readid(infn, outfn, sep='\t', readid_col=0, chr_col=1, pos_col=2, strand_col=3, score_col=4,
                                 label_col=5):
    """
    Convert meth per read score results into bed for IGV view

    Sample format:
    ID      Chr     Pos     Strand  Score   Label
    0a0a8a8b-3671-4265-8b59-bc30f633e42b    chr11   104837006       -       -3.2771419743285506     c
    0a0a8a8b-3671-4265-8b59-bc30f633e42b    chr11   104837082       -       4.224159732002238       m
    0a0a8a8b-3671-4265-8b59-bc30f633e42b    chr11   104837370       -       3.7890667822438067      m
    Args:
        sep:
        readid_col:

    Returns:

    """
    logger.debug(f"process file:{infn}")
    infile, lines = open_file_gz_or_txt(infn)
    retDict = {}
    for row in tqdm(infile, total=lines, desc="Import-perReadScore"):
        if row.startswith('ID\tChr'):
            continue
        tmp = row.strip().split(sep)
        readid = tmp[readid_col]
        chr = tmp[chr_col]
        strand = tmp[strand_col]
        # input is 1-based, in-memory will be 0-based, for - strand, also point to CG's C, this position will be modified by NanomethPhase later
        pos = int(tmp[pos_col]) - (1 if strand == '+' else 2)
        score = float(tmp[score_col])

        key = (readid, chr, strand)
        if key not in retDict:
            retDict[key] = defaultdict(list)
        if score > 0.0:
            retDict[key]['methylate_sites'].append(pos)
            retDict[key]['llr_score_meth'].append(score)
        else:
            retDict[key]['unmethylate_sites'].append(pos)
            retDict[key]['llr_score_unmeth'].append(score)
    infile.close()
    logger.debug(f"### read in DONE")

    outf = gzip.open(outfn, 'wt')
    # outf.write(f"ID\tChr\tPos\tStrand\tScore\tLabel\n")
    for key in tqdm(retDict, desc="Save-file"):
        readid = key[0]
        chr = key[1]
        strand = key[2]
        all_positions = sorted(retDict[key]['methylate_sites'] + retDict[key]['unmethylate_sites'])
        start = str(all_positions[0])  # output is 0-based start, 1-based end
        end = str(all_positions[-1] + 1)
        if len(retDict[key]['methylate_sites']) == 0:
            retDict[key]['methylate_sites'].append('NA')
            retDict[key]['llr_score_meth'].append('NA')

        if len(retDict[key]['unmethylate_sites']) == 0:
            retDict[key]['unmethylate_sites'].append('NA')
            retDict[key]['llr_score_unmeth'].append('NA')

        llr_meth_str = ','.join(str(x) for x in retDict[key]['llr_score_meth'])
        llr_unmeth_str = ','.join(str(x) for x in retDict[key]['llr_score_unmeth'])
        meth_sites_str = ','.join(str(x) for x in retDict[key]['methylate_sites'])
        unmeth_sites_str = ','.join(str(x) for x in retDict[key]['unmethylate_sites'])

        out_line = '\t'.join(
            (chr, start, end, strand, readid, llr_meth_str, llr_unmeth_str, meth_sites_str, unmeth_sites_str))

        outf.write(f"{out_line}\n")
    outf.close()
    logger.info(f"save to {outfn}")


def parse_arguments():
    parser = argparse.ArgumentParser(prog='process_per_read_score (NANOME)',
                                     description='Convert per read score results into bed format group by read-id.')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', type=str, help='input file for sorted per-read score file', required=True)
    parser.add_argument('-o', type=str, help="output file for bed.gz file", required=True)
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()
    logger.debug(args)

    process_read_score_by_readid(args.i, args.o)
    logger.info(f"### process_per_read_score DONE")
