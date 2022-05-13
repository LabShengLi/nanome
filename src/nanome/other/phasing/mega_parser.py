#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : mega_parser.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Read and parse megalodon/nanome per-read results
Output files:
    1. perReadScore         1-based
    2. read_pred            1-based
    3. site_pred            1-based
"""
import argparse
import gzip
import os
from collections import defaultdict

import math
import numpy as np
from tqdm import tqdm

from nanome.common.eval_common import open_file_gz_or_txt
from nanome.common.global_config import set_log_debug_level, logger, pic_base_dir
from nanome.common.global_settings import EPSLONG, NANOME_VERSION

indicator_5c = 0
indicator_5mc = 1
indicator_5hmc = 2

class_label = {
    indicator_5c: 'c',
    indicator_5mc: 'm',
    indicator_5hmc: 'h',
}


def import_megalodon_per_read_file(infn, chr_filter=None, readid_filter=None, sep='\t', readid_col=0, chr_col=1,
                                   strand_col=2, pos_col=3,
                                   meth_log_prob_col=4,
                                   unmeth_log_prob_col=5, meth_type_col=6, prob_cutoff=0.80, include_score=False,
                                   outBase=1, only_test=None, save_unified_format=False, outfn=None, num_class=2):
    """
    Import read-level megalodon, compatible with 5mc or 5hmc_5mc input
    Args:
        infn:
        readid_filter: must be set type for fast running!!!
        sep:
        readid_col:

    Returns: CPG -> [0 1 0 1 2], or [ (0, 0.88), (1, 0.99) ], in which 0-5c, 1-5mc, 2-5hmc
    Note: no need to check if 5hmc_5mc or 5mc

    Input samples:
    2d39b6b0-75b4-4ef8-8c6b-75871c924fd2	chr1	+	412238	-7.546572208404541	-8.03528033495954	h
    2d39b6b0-75b4-4ef8-8c6b-75871c924fd2	chr1	+	412238	-0.0008521132986061275	-8.03528033495954	m

    We validated 0-based start for Megalodon
    """
    if outBase not in [0, 1]:
        raise Exception(f"outBase={outBase} is not allowed")

    if num_class not in [2, 3]:
        raise Exception(f"num_class={num_class} is not allowed")

    if save_unified_format:
        if not outfn.endswith('tsv.gz'):
            raise Exception(f"outfn is not correct named, outfn={outfn}")
        logger.debug(f"Save unified format to {outfn}")
        outf1 = gzip.open(outfn, 'wt')
        outf1.write(f"ID\tChr\tPos\tStrand\tScore\tLabel\n")
        if num_class > 2:
            # seperate 5mc and 5hmc
            outf2 = gzip.open(outfn.replace(".tsv.gz", "_5hmc.tsv.gz"), 'wt')
            outf2.write(f"ID\tChr\tPos\tStrand\tScore\tLabel\n")

    if chr_filter is not None:  # set is fast than list
        chr_filter = set(chr_filter)

    if readid_filter is not None:  # set is fast than list
        readid_filter = set(readid_filter)

    cpgDict = defaultdict(list)
    infile, lines = open_file_gz_or_txt(infn)
    nreads = 0  # count pred prob > cutoff
    for row in tqdm(infile, total=lines, desc="Import-Megalodon"):
        if row.startswith('read_id\t'):
            continue
        if only_test and nreads >= only_test:
            break
        tmp = row.strip().split(sep)
        readid = tmp[readid_col]
        if readid_filter is not None and readid not in readid_filter:
            continue

        chr = tmp[chr_col]
        if chr_filter is not None and chr not in chr_filter:
            continue

        strand = tmp[strand_col]
        pos = int(tmp[pos_col]) + outBase
        meth_prob = float(np.e ** float(tmp[meth_log_prob_col]))
        unmeth_prob = float(np.e ** float(tmp[unmeth_log_prob_col]))
        meth_type = tmp[meth_type_col]

        if meth_prob + unmeth_prob > 1.0 + EPSLONG:
            raise Exception(
                f"Assert prob sum <= 1 failed, please check log transform functions or file correctness, meth_prob={meth_prob:.3f}, unmeth_prob={unmeth_prob:.3f}")

        if meth_prob >= prob_cutoff or (unmeth_prob >= prob_cutoff and meth_type == 'm'):
            key = (chr, pos, strand)
            if meth_prob >= prob_cutoff:
                prob_value = meth_prob
                if meth_type == 'm':
                    meth_indicator = indicator_5mc
                elif meth_type == 'h':
                    meth_indicator = indicator_5hmc
                else:
                    raise Exception(f"Not correct meth_type={meth_type}")
            elif unmeth_prob >= prob_cutoff:
                prob_value = unmeth_prob
                if meth_type == 'm':
                    meth_indicator = indicator_5c
                else:
                    raise Exception(f"Code bug, meth_type={meth_type} is not allowed here")

            if include_score:
                cpgDict[key].append((meth_indicator, prob_value))
            else:
                cpgDict[key].append(meth_indicator)
            nreads += 1

            ## keep only above cutoff reads
            if save_unified_format:
                score = math.log((meth_prob + EPSLONG) / (unmeth_prob + EPSLONG))
                if class_label[meth_indicator] == 'c':  # keep 5c to both files
                    outf1.write(
                        f"{readid}\t{chr}\t{pos}\t{strand}\t{score}\t{class_label[meth_indicator]}\n")
                    if num_class > 2:
                        outf2.write(
                            f"{readid}\t{chr}\t{pos}\t{strand}\t{score}\t{class_label[meth_indicator]}\n")
                elif class_label[meth_indicator] == 'm':  # keep 5mc to 5mc file
                    outf1.write(
                        f"{readid}\t{chr}\t{pos}\t{strand}\t{score}\t{class_label[meth_indicator]}\n")
                elif class_label[meth_indicator] == 'h' and num_class > 2:  # keep 5hmc to 5hmc file
                    outf2.write(
                        f"{readid}\t{chr}\t{pos}\t{strand}\t{score}\t{class_label[meth_indicator]}\n")

    if save_unified_format:
        outf1.close()
        logger.debug(f'Save unified output format to {outfn}')
        if num_class > 2:
            outf2.close()
            logger.debug(f'Save unified output format for both 5mc and 5hmc')

    logger.debug(f"import Megalodon: cpgs={len(cpgDict):,}\t read_level_preds={nreads:,}")
    return cpgDict


def import_nanome_per_read_file(infn, chr_filter=None, readid_filter=None, sep='\t', readid_col=0, chr_col=1,
                                strand_col=3, pos_col=2,
                                meth_prob_col=-1,
                                meth_indicator_col=-2, include_score=False,
                                outBase=1, only_test=None, save_unified_format=False, outfn=None):
    """
    Import read-level nanome input
    Args:
        infn:
        readid_filter: must be set type for fast running!!!
        sep:
        readid_col:

    Returns: CPG -> [0 1 0 1 2], or [ (0, 0.88), (1, 0.99) ], in which 0-5c, 1-5mc, 2-5hmc

    Input samples:
    ID	Chr	Pos	Strand	nanopolish	megalodon	deepsignal	Prediction	Prob_methylation
    1f3bc0e1-294c-4f1e-9986-a460d904ab49	chr12	10320	-	-0.59	1.5610911471496811		1	0.6527313
    38f9b801-1ceb-41a5-995b-7fe1797e1257	chr12	10320	-	-0.11	1.2781343948949806		1	0.6100183
    24367ea9-6b72-4f4f-847f-18c184331fb3	chr12	10320	-	-0.06		00.49295574

    We checked input as 1-based format for start col
    """
    if outBase not in [0, 1]:
        raise Exception(f"outBase={outBase} is not allowed")

    if save_unified_format:
        logger.debug(f"Save unified format to {outfn}")
        outf = gzip.open(outfn, 'wt')
        outf.write(f"ID\tChr\tPos\tStrand\tScore\tLabel\n")

    if chr_filter is not None:  # set is fast than list
        chr_filter = set(chr_filter)

    if readid_filter is not None:  # set is fast than list
        readid_filter = set(readid_filter)

    cpgDict = defaultdict(list)
    infile, lines = open_file_gz_or_txt(infn)
    nreads = 0  # count pred prob > cutoff
    for row in tqdm(infile, total=lines, desc="Import-NANOME"):
        if row.startswith('ID\tChr'):
            continue
        if only_test and nreads >= only_test:
            break
        tmp = row.strip().split(sep)
        readid = tmp[readid_col]
        if readid_filter is not None and readid not in readid_filter:
            continue

        chr = tmp[chr_col]
        if chr_filter is not None and chr not in chr_filter:
            continue

        strand = tmp[strand_col]
        pos = int(tmp[pos_col]) + (outBase - 1)
        meth_prob = float(tmp[meth_prob_col])
        meth_indicator = int(tmp[meth_indicator_col])

        key = (chr, pos, strand)
        if include_score:
            cpgDict[key].append((meth_indicator, meth_prob))
        else:
            cpgDict[key].append(meth_indicator)
        nreads += 1

        ## keep only above cutoff reads
        if save_unified_format:
            # output to 1-based for meteore, ref: https://github.com/comprna/METEORE/blob/master/script/format_megalodon.R#L16
            score = math.log((meth_prob + EPSLONG) / (1 - meth_prob + EPSLONG))
            outf.write(
                f"{readid}\t{chr}\t{pos}\t{strand}\t{score}\t{class_label[meth_indicator]}\n")

    if save_unified_format:
        outf.close()
        logger.debug(f'Save unified output format to {outfn}')

    logger.debug(f"import NANOME: cpgs={len(cpgDict):,}\t read_level_preds={nreads:,}")
    return cpgDict


def count_class_occurs(alist, class_label):
    """
    Count how many preds are for class_label (0, 1 or 2)
    Args:
        alist: may be int list, or tuple list for (class_label, read_pred_prob)
        class_label:

    Returns:

    """
    if type(alist[0]) == tuple:
        newlist = [k[0] for k in alist]
    elif type(alist[0]) == int:
        newlist = alist
    else:
        raise Exception(f"Element in alist is not support, type(alist[0])={type(alist[0])}, alist[0]={alist[0]}")
    return newlist.count(class_label)


def agg_read_to_site(cpgDict, num_class=3):
    """
    Aggregate read cpgDict into site freq dict, must know if 5hmc_5mc or 5mc input before convert to site freq
    Args:
        cpgDict:  CPG -> [0 1 0 1 2]
        num_class: 2 or 3

    Returns:
        siteDict:   CPG -> (0.3, 0.7, cov) / (0.1, 0.2, 0.7, cov)
    """
    siteDict = {}
    for key in cpgDict:
        freq_5c = count_class_occurs(cpgDict[key], indicator_5c) / len(cpgDict[key])
        freq_5mc = count_class_occurs(cpgDict[key], indicator_5mc) / len(cpgDict[key])
        freq_5hmc = count_class_occurs(cpgDict[key], indicator_5hmc) / len(cpgDict[key])
        if num_class == 3:
            siteDict[key] = (freq_5c, freq_5mc, freq_5hmc, len(cpgDict[key]))
        elif num_class == 2:
            siteDict[key] = (freq_5c, freq_5mc, len(cpgDict[key]))
        else:
            raise Exception(f"params num_class is not allowed, num_class={num_class}")
    return siteDict


def to_read_preds_file(cpgDict, outfn, inBase=1, is_tuple_element=False):
    outf = gzip.open(outfn, 'wt')
    outf.write(f"Chr\tPos_base{inBase}\tStrand\tPreds\tCoverage\n")

    for key in cpgDict:
        if is_tuple_element:
            preds_str = ','.join(f"{x[0]}:{x[1]:.3f}" for x in cpgDict[key])
        else:
            preds_str = ','.join(str(x) for x in cpgDict[key])
        outf.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{preds_str}\t{len(cpgDict[key])}\n")
    outf.close()


def to_site_freq_file(siteDict, outfn, inBase=1, num_class=3):
    outf = gzip.open(outfn, 'wt')

    if num_class == 3:
        outf.write(f"Chr\tPos_base{inBase}\tStrand\tP5c\tP5mc\tP5hmc\tCoverage\n")
    elif num_class == 2:
        outf.write(f"Chr\tPos_base{inBase}\tStrand\tP5c\tP5mc\tCoverage\n")
    else:
        raise Exception(f"param num_class is not correct, num_class={num_class}")

    for key in siteDict:
        probs_str = '\t'.join(f"{x:.3f}" for x in siteDict[key][:-1])
        probs_str += f"\t{siteDict[key][-1]}"
        outf.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{probs_str}\n")
    outf.close()


def parse_arguments():
    parser = argparse.ArgumentParser(prog='mega_parser (NANOME)', description='Parse megalodon per-read results.')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('-i', type=str, help='input file for megalodon per-read', required=True)
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('-o', type=str, help="output dir",
                        default=None)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    set_log_debug_level()
    args = parse_arguments()
    logger.debug(args)

    if args.o is None:
        outdir = pic_base_dir
    else:
        outdir = args.o
        os.makedirs(outdir, exist_ok=True)
    logger.debug(f"outdir={outdir}")

    # dsname = "APL"
    # infn = "/fastscratch/liuya/nanome/APL_analysis/APL_methcall/APL_megalodon_per_read_sort.tsv.gz"

    # dsname = "CD34_Human"
    # infn = "/fastscratch/liuya/nanome/CD34_Human_analysis/CD34_Human_methcall/CD34_Human_megalodon_per_read_sort.tsv.gz"

    dsname = args.dsname
    infn = args.i

    predDict = import_megalodon_per_read_file(infn)
    siteDict = agg_read_to_site(predDict)

    outfn = os.path.join(pic_base_dir, f'{dsname}_megalodon_read_pred.tsv.gz')
    to_read_preds_file(predDict, outfn)

    outfn = os.path.join(pic_base_dir, f'{dsname}_megalodon_site_freq.tsv.gz')
    to_site_freq_file(siteDict, outfn)
