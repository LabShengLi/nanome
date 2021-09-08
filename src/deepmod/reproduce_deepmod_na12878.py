#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : reproduce_deepmod.na12878.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/liuyangzzu/nanome

"""
Sanity check NA12878 for DeepMod paper.
"""
import argparse
import gzip
import logging
import os
import sys
from collections import defaultdict

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix

# chr 1-22 X and Y
humanChrSet = [f'chr{k}' for k in range(1, 23)] + ['chrX', 'chrY']


def open_file_gz_or_txt(infn):
    """
    Open a txt or gz file, based on suffix
    :param infn:
    :return:
    """
    logger.debug(f"open file: {infn}")
    if infn.endswith('.gz'):
        infile = gzip.open(infn, 'rt')  # using rt option, no need to convert bytearray
    else:
        infile = open(infn, 'r')
    return infile


def import_encode_bsseq(infileName, chr_col=0, start_col=1, methfreq_col=-1, cov_col=4, strand_col=5,
                        covCutoff=1, baseFormat=1, filterChr=humanChrSet, includeCov=True):
    """
    ### Description of the columns in this format (https://www.encodeproject.org/data-standards/wgbs/):

    1. Reference chromosome or scaffold
    2. Start position in chromosome (0-based position)
    3. End position in chromosome (1-based position)
    4. Name of item
    5. Score from 0-1000. Capped number of reads
    6. Strandedness, plus (+), minus (-), or unknown (.)
    7. Start of where display should be thick (start codon)
    8. End of where display should be thick (stop codon)
    9. Color value (RGB)
    10. Coverage, or number of reads
    11. Percentage of reads that show methylation at this position in the genome

    We check the input is 0-based.

    e.g.:
    chr17   115761  115762  .       0       -       115761  115762  0,255,0 0       0
    chr17   116083  116084  .       3       +       116083  116084  0,255,0 3       0
    chr17   116084  116085  .       7       -       116084  116085  0,255,0 7       0
    chr17   116353  116354  .       9       +       116353  116354  0,255,0 9       0
    chr17   116354  116355  .       1       -       116354  116355  0,255,0 1       0
    chr17   116445  116446  .       12      +       116445  116446  255,255,0       12      50
    chr17   116446  116447  .       1       -       116446  116447  0,255,0 1       0
    chr17   116703  116704  .       17      +       116703  116704  0,255,0 17      0
    chr17   116704  116705  .       4       -       116704  116705  0,255,0 4       0

    Note, that the first row above have coverage=0, so they list all CpGs (this is from WGBS data).
    This is not the case for RRBS, where they only list covered sites:

    e.g.:
    chr1    9943211 9943212 K562_Rep3_RRBS  168     -       10003269        10003270        0,255,0 168     0
    chr1    9943228 9943229 K562_Rep3_RRBS  168     -       10003286        10003287        0,255,0 168     1
    chr1    9943239 9943240 K562_Rep3_RRBS  1       +       10003297        10003298        0,255,0 1       0
    chr1    9943240 9943241 K562_Rep3_RRBS  168     -       10003298        10003299        0,255,0 168     4

    Return:
        key -> [meth_freq, meth_cov] in which,
                key:    (chr, start, strand)
                value:  meth_frep in [0,1], meth_cov is int
    """
    cpg_dict = {}

    infile = open_file_gz_or_txt(infileName)

    nrow = 0
    for row in infile:
        nrow += 1

        tmp = row.strip().split("\t")

        if filterChr is not None and tmp[chr_col] not in filterChr:  # Filter out non-human chrs
            continue

        strand = tmp[strand_col]

        if strand not in ['+', '-']:
            raise Exception(f'input file format error when parsing: {row}, strand={strand}')

        try:
            if baseFormat == 1:
                start = int(tmp[start_col]) + 1
            elif baseFormat == 0:
                start = int(tmp[start_col])
            else:
                logger.error(
                    "###\timportGroundTruth_BedMethyl_from_Encode InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(
                        baseFormat))
                sys.exit(-1)
            cov_cnt = int(tmp[cov_col])
            meth_freq = int(tmp[methfreq_col]) / 100.0  # [0,1]
        except:
            logger.error(f" ### Error parse gbTruth row = {row}")
            continue

        if cov_cnt < covCutoff:
            continue

        key = (tmp[chr_col], start, strand)

        if key in cpg_dict:
            raise Exception(f"Found duplicate CpG sites for {row}")

        if includeCov:
            cpg_dict[key] = (meth_freq, cov_cnt)
        else:
            cpg_dict[key] = meth_freq

    infile.close()
    logger.debug(
        f"###\timportGroundTruth_BedMethyl_from_Encode: loaded information for {len(cpg_dict):,} CpGs, with cutoff={covCutoff} ({nrow:,} rows)")

    return cpg_dict


def combine_bsseq_replicates(bgtruth_list, freqCutoff=0.9, covCutoff=1, filterChrs=humanChrSet):
    """
    Combine two replicates by DeepMod, >90% in both as methylated, =0% in both as unmethylated, remove others
    :param bgtruth_list:
    :return:
    """
    jointBGTruth = {}  # intersection of CpG sites
    if len(bgtruth_list) == 2:  # sites must in both replicates, and
        unionSet = set(bgtruth_list[0].keys()).union(set(bgtruth_list[1].keys()))
        jointSet = set(bgtruth_list[0].keys()).intersection(set(bgtruth_list[1].keys()))
        logger.debug(
            f'unionSet={len(unionSet):,}, jointSet={len(jointSet):,}, cov_cutoff={covCutoff}, for chr={filterChrs}')
        methCnt = unmethCnt = 0

        totalCnt = 0

        for key in jointSet:
            if key[0] not in filterChrs:
                continue
            totalCnt += 1
            cov1 = bgtruth_list[0][key][1]
            methfreq1 = bgtruth_list[0][key][0]

            cov2 = bgtruth_list[1][key][1]
            methfreq2 = bgtruth_list[1][key][0]

            if cov1 < covCutoff or cov2 < covCutoff:
                continue  # ensure both >= cov_cutoff

            if 1e-5 <= methfreq1 < freqCutoff or 1e-5 <= methfreq2 < freqCutoff:
                continue  # ensure both >90% or =0%

            if methfreq1 < 1e-5 and methfreq2 < 1e-5:
                meth_indicator = 0
                unmethCnt += 1
            elif methfreq1 >= freqCutoff and methfreq2 >= freqCutoff:
                meth_indicator = 1
                methCnt += 1
            else:
                continue

            jointBGTruth[key] = (meth_indicator, min(cov1, cov2))
        logger.debug(
            f'chrSet={filterChrs} with cov-cutoff={covCutoff}, jointBGTruth={len(jointBGTruth):,}, unmethCnt={unmethCnt:,}, methCnt={methCnt:,}, not-used={totalCnt - methCnt - unmethCnt:,}')
    else:
        raise Exception(f'len={len(bgtruth_list)}, is not support now.')

    ret = (methCnt, unmethCnt, totalCnt - methCnt - unmethCnt)
    return jointBGTruth, ret


def import_deepmod(infileName, chr_col=0, start_col=1, strand_col=5, coverage_col=-4,
                   meth_cov_col=-2, clustered_meth_freq_col=-1, baseFormat=1,
                   output_first=False, siteLevel=True, include_score=False, filterChr=humanChrSet,
                   total_cols=13):
    """
    DeepMod RNN+Cluster results format for human genome
    Note that DeepMod only outputs site level stats, we use DeepMod clustered results for site level evaluation only.
    This function have three outputs:
    1. read level ouput of [0 0 1 1 1 ]
    2. read level output of [(0, 0.0), (1, 1.0)]
    3. site level output of [(0.25, 25), (0.8, 12)], i.e.,  [methFrequency, coverage] such as key -> values [50 (freq 0-100), 10 (cov)]

    ### Example input format from DeepMod (clustered - following Step 4 from "Example 3: Detect 5mC on Na12878" section; https://github.com/WGLab/DeepMod/blob/master/docs/Reproducibility.md):
    chr2 241991445 241991446 C 3 -  241991445 241991446 0,0,0 3 100 3 69
    chr2 241991475 241991476 C 3 -  241991475 241991476 0,0,0 3 33 1 75
    chr2 241991481 241991482 C 2 -  241991481 241991482 0,0,0 2 50 1 76
    Note: it is white space separated, not tab-separated file
    ============
    """
    infile = open_file_gz_or_txt(infileName)

    cpg_dict = defaultdict()
    count_calls = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        tmp = row.strip().split()

        if len(tmp) != total_cols:
            raise Exception(f"DeepMod RNN+clustered output format error: tmp={tmp}")

        if tmp[chr_col] not in filterChr:
            continue

        if output_first:
            logger.debug(f'row = {list(enumerate(tmp))}')
            output_first = False

        if baseFormat == 1:
            start = int(tmp[start_col]) + 1
            strand = tmp[strand_col]
        elif baseFormat == 0:
            start = int(tmp[start_col])
            strand = tmp[strand_col]
        else:
            logger.error(
                "###\timportPredictions_DeepMod InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(
                    baseFormat))
            sys.exit(-1)

        if strand not in ['+', '-']:
            raise Exception(f'Strand info={strand} is not correctly recognized, row={row}, in file {infileName}')

        key = (tmp[chr_col], start, strand)

        meth_freq = int(tmp[clustered_meth_freq_col]) / 100.0
        coverage = int(tmp[coverage_col])
        meth_cov = int(tmp[meth_cov_col])

        if key in cpg_dict:
            raise Exception(f'In DeepMod_Cluster results, we found duplicate key={key}, this is not correct')

        if siteLevel:  # in dict
            cpg_dict[key] = (meth_freq, coverage)  # {'freq': meth_freq, 'cov': coverage}
        elif include_score:  # For read-level include scores
            cpg_dict[key] = [(1, 1.0)] * meth_cov + [(0, 0.0)] * (coverage - meth_cov)
        else:  # For read-level no scores
            cpg_dict[key] = [1] * meth_cov + [0] * (coverage - meth_cov)

        count_calls += coverage
        meth_cnt += meth_cov
        unmeth_cnt += coverage - meth_cov

    infile.close()

    logger.info(
        f"###\timportDeepMod_clustered Parsing SUCCESS: {count_calls:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpg_dict):,} CpGs from {infileName} file")

    return cpg_dict


def parse_arguments():
    """
    Sample usage:

    python /projects/li-lab/yang/workspace/nano-compare/src/deepmod/reproduce_deepmod_na12878.py --chr chr20 --deepmod-input /projects/li-lab/Nanopore_compare/suppdata/deepmod-albacore/na12878_pred/ecoli_org/cpredecoli_org_clusterCpG.chr20.C.bed --bgtruth '/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz;/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz' --pred-threshold 0.5

    :return:
    """
    parser = argparse.ArgumentParser(description='Sanity check deepmod paper')
    parser.add_argument('--deepmod-input', type=str, help='DeepMod format output results', required=True)
    parser.add_argument('--bgtruth', type=str, help="BS-seq two replicates files: <file-name1>;<file-name1>",
                        required=True)
    parser.add_argument('--chr', type=str, help="chromosome to evaluate", required=True)
    parser.add_argument('--processors', type=int, help="running processors", default=4)
    parser.add_argument('--pred-threshold', type=float, help="prediction threshold for methylation level", default=0.5)
    parser.add_argument('-o', type=str, help="output dir", default='.')
    parser.add_argument('--report-joined', help="if also report joined performance", action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - [%(filename)s:%(lineno)d] - %(levelname)s: %(message)s')
    logger = logging.getLogger('SanityDeepMod')

    args = parse_arguments()
    logger.debug(args)

    # out_dir = os.path.join(args.o, "DeepMod_paper_sanity_check")
    # os.makedirs(out_dir, exist_ok=True)
    # logger.info(f'Output to dir:{out_dir}')
    logger.info(f'\n\n####################\n\n')

    # we import multiple (1 or 2) replicates and join them
    fnlist = args.bgtruth.split(';')

    if len(fnlist) > 2:
        raise Exception(f'Currently only support bgtruth with upto two, but found more: {fnlist}')

    logger.debug(f'BGTruth fnlist={fnlist}')

    bgtruth_list = []
    for fn in fnlist:
        if len(fn) == 0:  # incase of input like 'bismark:/a/b/c;'
            continue
        # import if cov >= 1 firstly, then after join two replicates step, remove low coverage
        bgtruth = import_encode_bsseq(fn)
        bgtruth_list.append(bgtruth)

    logger.info(f'\n\n####################\n\n')

    ## Evaluate for chr
    logger.debug(f'Start check paper results for chr={args.chr}')
    chrSet = [args.chr]
    # Combine one/two replicates using DeepMod methods
    combined_bgtruth_cov1, ret1 = combine_bsseq_replicates(bgtruth_list, covCutoff=1, filterChrs=chrSet)
    # value is (meth_indicator: int, cov)
    combined_bgtruth_cov5 = {key: combined_bgtruth_cov1[key] for key in combined_bgtruth_cov1 if
                             combined_bgtruth_cov1[key][1] >= 5}

    # value is (meth_freq, cov)
    deepmod_calls = import_deepmod(args.deepmod_input)
    calls_set = set(deepmod_calls.keys())
    bsseq_set = set(combined_bgtruth_cov5.keys())
    joined_set = calls_set.intersection(bsseq_set)
    logger.info(
        f"BS-seq(cov>=5)={len(combined_bgtruth_cov5):,}, DeepMod={len(deepmod_calls):,}, intersect={len(joined_set):,}")

    if args.report_joined:
        ## joined sets performance report
        y_truth = []
        y_pred = []
        for key in joined_set:
            bs_label = combined_bgtruth_cov5[key][0]
            pred_label = 1 if deepmod_calls[key][0] >= args.pred_threshold else 0
            y_truth.append(bs_label)
            y_pred.append(pred_label)
        ## Calculate precision and recall
        prec = precision_score(y_truth, y_pred)
        recall = recall_score(y_truth, y_pred)
        f1 = f1_score(y_truth, y_pred)
        accuracy = accuracy_score(y_truth, y_pred)
        logger.info(
            f"Joined report: chr={args.chr} (Meth={sum(y_truth):,}, Unmeth={len(y_truth) - sum(y_truth)}), precision={prec}, recall={recall}, f1={f1}, accuracy={accuracy}")
        print(
            f"Type\tChr\tMeth\tUnmeth\tPrecision\tRecall\tF1-score\tAccuracy")
        print(
            f"Joined report\t{args.chr}\t{sum(y_truth)}\t{len(y_truth) - sum(y_truth)}\t{prec}\t{recall}\t{f1}\t{accuracy}")

    y_truth = []
    y_pred = []
    for key in combined_bgtruth_cov5:
        bs_label = combined_bgtruth_cov5[key][0]
        if key in deepmod_calls:
            pred_label = 1 if deepmod_calls[key][0] >= args.pred_threshold else 0
        else:
            pred_label = 0
        y_truth.append(bs_label)
        y_pred.append(pred_label)
    ## Calculate precision and recall
    prec = precision_score(y_truth, y_pred)
    recall = recall_score(y_truth, y_pred)
    f1 = f1_score(y_truth, y_pred)
    accuracy = accuracy_score(y_truth, y_pred)
    conf_matrix = confusion_matrix(y_truth, y_pred, labels=[1,0])
    logger.info(
        f"BS-seq report: chr={args.chr} (Meth={sum(y_truth):,}, Unmeth={len(y_truth) - sum(y_truth):,}), precision={prec}, recall={recall}, f1={f1}, accuracy={accuracy}")
    print(
        f"Type\tChr\tMeth\tUnmeth\tPrecision\tRecall\tF1-score\tAccuracy")
    print(
        f"Bs-seq report\t{args.chr}\t{sum(y_truth)}\t{len(y_truth) - sum(y_truth)}\t{prec}\t{recall}\t{f1}\t{accuracy}")
    print(f"Confusion matrix=\n{conf_matrix}")

    logger.info("### Sanity check DONE")
