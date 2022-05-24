#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : eval_common.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Common functions used by read level and site level evaluations in nanome paper.

Such as import_DeepSignal, import_BGTruth, etc.
"""

import glob
import gzip
import os.path
import pickle
import re
import subprocess
import sys
import warnings
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, wait
from itertools import combinations
from multiprocessing import Pool

import math
import numpy as np
import pandas as pd
import psutil
import pysam
from Bio import SeqIO
from pybedtools import BedTool
from sklearn.metrics import roc_curve, auc, average_precision_score, f1_score, precision_score, recall_score
from tqdm import tqdm

from nanome.common.global_config import *
from nanome.common.global_settings import HUMAN_CHR_SET, ToolEncodeList, BGTruthEncodeList, reference_genome_hg38_fn, \
    enable_base_detection_bedfile, region_filename_dict, genome_wide_tagname, EPSLONG, CHUNKSIZE


def importPredictions_Nanopolish(infileName, chr_col=0, start_col=2, strand_col=1, readid_col=4, log_lik_ratio_col=5,
                                 sequence_col=-1, num_motifs_col=-2, baseFormat=1, score_cutoff=(-2.0, 2.0),
                                 output_first=False,
                                 include_score=False, filterChr=HUMAN_CHR_SET, save_unified_format=False, outfn=None,
                                 stringent_cutoff=True):
    """
    We checked the input is 0-based for the start col
    Return dict of key='chr1\t123\t123\t+', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.

    ### Example default input format from Nanopolish (pre-processed to containing strand-info):
    head /projects/li-lab/yang/results/12-11/K562.methylation_calls-nanopolish-strand-info.tsv
    chromosome	start	end	read_name	log_lik_ratio	log_lik_methylated	log_lik_unmethylated	num_calling_strands	num_cpgs	sequence	strand-info
    chr1	24450	24450	23cf5aac-6664-4fdb-9334-5b7e55f33335	-8.02	-117.7	-109.68	1	1	GAAAACGTGAA	-
    chr1	24553	24553	23cf5aac-6664-4fdb-9334-5b7e55f33335	4.4	-142.26	-146.66	1	1	GCTCTCGGACT	-
    chr1	24637	24637	23cf5aac-6664-4fdb-9334-5b7e55f33335	-8.97	-162.65	-153.68	1	1	AGGACCGGGAT	-
    chr1	24784	24784	23cf5aac-6664-4fdb-9334-5b7e55f33335	1.78	-137.14	-138.91	1	1	GCATCCGCCAT	-
    chr1	24809	24812	23cf5aac-6664-4fdb-9334-5b7e55f33335	8.49	-164.69	-173.18	1	2	CCTCTCGCCGCAGG	-
    chr1	24837	24837	23cf5aac-6664-4fdb-9334-5b7e55f33335	12.53	-182.95	-195.47	1	1	GGGCACGGCAT	-
    chr1	24899	24910	23cf5aac-6664-4fdb-9334-5b7e55f33335	4.19	-215.96	-220.15	1	3	TGGGTCGGAGCCGGAGCGTCAG	-
    chr1	24925	24937	23cf5aac-6664-4fdb-9334-5b7e55f33335	25.63	-213.44	-239.06	1	3	ACCCACGACCACCGGCACGCCCC	-
    ###############
    Referenced by Nanopolish script for handling with conversion of calls to frequencies:
    https://github.com/jts/nanopolish/blob/master/scripts/calculate_methylation_frequency.py
    """
    filterChr = set(filterChr)
    if score_cutoff is None:
        score_cutoff = (-2.0, 2.0)
    llr_cutoff = abs(score_cutoff[1])

    cpgDict = defaultdict(list)
    call_cnt = 0
    meth_cnt = 0
    unmeth_cnt = 0

    infile, lines = open_file_gz_or_txt(infileName)

    if save_unified_format:
        # outf = open(outfn, 'w')
        outf = gzip.open(outfn, 'wt')
        outf.write(f"ID\tChr\tPos\tStrand\tScore\n")

    for row in tqdm(infile, total=lines, desc="Import-Nanopolish"):
        tmp = row.strip().split("\t")

        if tmp[chr_col] != "chromosome":
            if output_first:
                logger.debug(list(enumerate(tmp)))
                output_first = False

            if tmp[chr_col] not in filterChr:
                continue

            try:  # try to find if these columns are interpretable
                start = int(tmp[start_col])
                num_sites = int(tmp[num_motifs_col])
                llr = float(tmp[log_lik_ratio_col])

                meth_score = llr
                if llr > 0:
                    meth_indicator = 1

                elif llr < 0:
                    meth_indicator = 0

                strand_info = tmp[strand_col]
                if strand_info == '-':  # - strand, point to positive sequence C, so need +1 to point to G
                    start = start + 1

                if strand_info not in ['-', '+']:
                    raise Exception(
                        f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')
            except:
                logger.error(f'###\tError when parsing row=[{row}] in {infileName}')
                continue

            if stringent_cutoff:
                final_cutoff = llr_cutoff * num_sites
            else:
                final_cutoff = llr_cutoff

            if num_sites == 1:  # we have singleton, i.e. only one CpG within the area
                if baseFormat == 0:
                    key = (tmp[chr_col], start, strand_info)
                elif baseFormat == 1:
                    key = (tmp[chr_col], start + 1, strand_info)
                else:
                    logger.error(
                        "###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(
                            baseFormat))
                    sys.exit(-1)

                if save_unified_format:
                    # output to 1-based for meteore, ref: https://github.com/comprna/METEORE/blob/master/script_in_snakemake/format_nanopolish.R#L14
                    outf.write(f"{tmp[readid_col]}\t{tmp[chr_col]}\t{start + 1}\t{tmp[strand_col]}\t{meth_score}\n")

                if abs(llr) < final_cutoff:  # Consider all sites as a group when there are multiple sites
                    continue

                if include_score:
                    cpgDict[key].append((meth_indicator, meth_score))
                else:
                    cpgDict[key].append(meth_indicator)

                call_cnt += 1
                if meth_indicator == 1:
                    meth_cnt += 1
                else:
                    unmeth_cnt += 1
            else:  # we deal with non-singleton
                firstCpgLoc = int(tmp[start_col]) - 5
                sequence = tmp[sequence_col]
                for cpg in re.finditer("CG", sequence):
                    cpgStart = cpg.start() + firstCpgLoc
                    if strand_info == '-':
                        cpgStart = cpgStart + 1
                    elif strand_info != '+':
                        raise Exception(f'The file [{infileName}] contains no strand-info, please check it')

                    if baseFormat == 0:
                        key = (tmp[chr_col], cpgStart, strand_info)
                    elif baseFormat == 1:
                        key = (tmp[chr_col], cpgStart + 1, strand_info)
                    else:
                        logger.error(
                            "###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(
                                baseFormat))
                        sys.exit(-1)

                    if save_unified_format:
                        # output to 1-based for meteore, ref: https://github.com/comprna/METEORE/blob/master/script_in_snakemake/format_nanopolish.R#L14
                        outf.write(
                            f"{tmp[readid_col]}\t{tmp[chr_col]}\t{cpgStart + 1}\t{tmp[strand_col]}\t{meth_score}\n")

                    if abs(llr) < final_cutoff:  # Consider all sites as a group when there are multiple sites
                        continue

                    if include_score:
                        cpgDict[key].append((meth_indicator, meth_score))
                    else:
                        cpgDict[key].append(meth_indicator)
                    call_cnt += 1
                    if meth_indicator == 1:
                        meth_cnt += 1
                    else:
                        unmeth_cnt += 1

    infile.close()

    if save_unified_format:
        outf.close()
        logger.debug(f'Save METEORE output format to {outfn}')

    logger.debug(
        f"###\timportPredictions_Nanopolish SUCCESS: {call_cnt:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file with score_cutoff={score_cutoff}")
    return cpgDict


def importPredictions_DeepSignal(infileName, chr_col=0, start_col=1, strand_col=2, readid_col=4, meth_prob_col=7,
                                 meth_col=8, baseFormat=1, score_cutoff=(0.5, 0.5), include_score=False,
                                 filterChr=HUMAN_CHR_SET,
                                 save_unified_format=False, outfn=None):
    """
    We checked input as 0-based format for start col.
    Return dict of key='chr1\t123\t123\t+', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.
    Note that the function requires per read stats, not frequencies of methylation.
    ### Example input format from DeepSignal:
    head /projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/K562/K562.deepsignal_MethCalls.tsv
    chr9	124877811	+	124877811	fd575bfa-96d2-41f6-852e-b3cd3b67a8e4	t	0.7913327	0.2086673	0	CCTGGTCACGTCTCCTG
    chr9	124877883	+	124877883	fd575bfa-96d2-41f6-852e-b3cd3b67a8e4	t	0.20028716	0.79971284	1	GAACTAAACGTCAGAAA
    chr9	124878198	+	124878198	fd575bfa-96d2-41f6-852e-b3cd3b67a8e4	t	0.15723002	0.84277	1	TTAAATTACGTATATTT
    chr22	22922449	+	22922449	38af944b-5bbf-402f-8fc9-90ba4f3392f0	t	0.6276733	0.37232664	0	CCACTCACCGCTGACCT
    chr22	22922479	+	22922479	38af944b-5bbf-402f-8fc9-90ba4f3392f0	t	0.84463745	0.15536247	0	CAAGGGTCCGGCCTGAG
    chr22	22922588	+	22922588	38af944b-5bbf-402f-8fc9-90ba4f3392f0	t	0.46641168	0.5335883	1	ATCCACCCCGCAGGTCA
    chr10	71503608	-	62293813	6579588c-6785-4cd0-ada8-c6408302aaa1	t	0.69284683	0.3071532	0	GCCCATCACGCAGCACA
    chr10	71503427	-	62293994	6579588c-6785-4cd0-ada8-c6408302aaa1	t	0.8272412	0.17275886	0	AGCCACAACGGGAAGAG
    ### Input file format description:
    - chrom: the chromosome name
    - pos: 0-based position of the targeted base in the chromosome
    - strand: +/-, the aligned strand of the read to the reference
    - pos_in_strand: 0-based position of the targeted base in the aligned strand of the chromosome
    - readname: the read name
    - read_strand: t/c, template or complement
    - prob_0: [0, 1], the probability of the targeted base predicted as 0 (unmethylated)
    - prob_1: [0, 1], the probability of the targeted base predicted as 1 (methylated)
    - called_label: 0/1, unmethylated/methylated
    - k_mer: the kmer around the targeted base
    ** by default if this probability will be higher than > 0.5, DeepSignal will tell that this is methylated site, or else is unmethylated
    """
    filterChr = set(filterChr)
    if score_cutoff is None:
        score_cutoff = (0.5, 0.5)
    infile, lines = open_file_gz_or_txt(infileName)

    if save_unified_format:
        outf = gzip.open(outfn, 'wt')
        outf.write(f"ID\tChr\tPos\tStrand\tScore\n")

    cpgDict = defaultdict(list)
    row_count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in tqdm(infile, total=lines, desc="Import-DeepSignal"):
        tmp = row.strip().split("\t")

        if tmp[chr_col] not in filterChr:
            continue

        if baseFormat == 1:
            start = int(tmp[start_col]) + 1
            strand = tmp[strand_col]
        elif baseFormat == 0:
            start = int(tmp[start_col])
            strand = tmp[strand_col]
        else:
            logger.error(
                "###\timportPredictions_DeepSignal InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(
                    baseFormat))
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')

        if save_unified_format:
            # output to 1-based for meteore, ref: https://github.com/comprna/METEORE/blob/master/script_in_snakemake/format_deepsignal.R#L19
            try:
                meteore_score = math.log2((float(tmp[meth_prob_col]) + 1e-4) / (float(tmp[meth_prob_col - 1]) + 1e-4))
                outf.write(f"{tmp[readid_col]}\t{tmp[chr_col]}\t{start}\t{tmp[strand_col]}\t{meteore_score}\n")
            except:
                logger.warn(f"METEORE score calculate for DeepSignal error for row:{row}")

        key = (tmp[chr_col], start, strand)

        if include_score:
            the_score = float(tmp[meth_prob_col])
            if np.isnan(the_score):
                the_score = 0.0
            cpgDict[key].append((int(tmp[meth_col]), the_score,))
        else:
            cpgDict[key].append(int(tmp[meth_col]))

        row_count += 1
        if int(tmp[meth_col]) == 1:
            meth_cnt += 1
        else:
            unmeth_cnt += 1

    infile.close()

    if save_unified_format:
        outf.close()
        logger.debug(f'Save METEORE output format to {outfn}')

    logger.debug(
        f"###\timportPredictions_DeepSignal SUCCESS: {row_count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file with score_cutoff={score_cutoff}")

    return cpgDict


def importPredictions_Tombo(infileName, chr_col=0, start_col=1, readid_col=3, strand_col=5, meth_col=4, baseFormat=1,
                            score_cutoff=(-1.5, 2.5), output_first=False, include_score=False, filterChr=HUMAN_CHR_SET,
                            save_unified_format=False, outfn=None):
    """
    We checked input as 0-based start format.
    Return dict of key='chr1 123 +', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.
    Note that the function requires per read stats, not frequencies of methylation.

    ### Example input format from Tombo
    chr1    48020    48020    3526811b-6958-49f8-b78c-a205c1b5fc6e    1.185219591257949    +    TATTACACCCG
    chr1    48022    48022    3526811b-6958-49f8-b78c-a205c1b5fc6e    1.6267354150537658    +    TTACACCCGTT
    chr1    48023    48023    3526811b-6958-49f8-b78c-a205c1b5fc6e    2.6122662196889728    +    TACACCCGTTA
    chr1    48024    48024    3526811b-6958-49f8-b78c-a205c1b5fc6e    2.771131774766473    +    ACACCCGTTAA
    chr1    48041    48041    3526811b-6958-49f8-b78c-a205c1b5fc6e    6.524775544143312    +    GATTTCTAAAT
    chr1    48048    48048    3526811b-6958-49f8-b78c-a205c1b5fc6e    1.9142728191641216    +    AAATGCATTGA
    chr1    48054    48054    3526811b-6958-49f8-b78c-a205c1b5fc6e    1.8675210090110548    +    ATTGACATTTG
    ......
    chr1    8447736    8447736    c9339e26-1898-4483-a312-b78c3fafc6a9    8.073560995614967    -    CTGTGCTGTGT
    chr1    8447745    8447745    c9339e26-1898-4483-a312-b78c3fafc6a9    2.4467964154940858    -    GTTGACCGTGT
    chr1    8447746    8447746    c9339e26-1898-4483-a312-b78c3fafc6a9    1.966921521322515    -    TTGACCGTGTA
    chr1    8447754    8447754    c9339e26-1898-4483-a312-b78c3fafc6a9    5.387457000225035    -    GTATGCAATGG
    chr1    8447761    8447761    c9339e26-1898-4483-a312-b78c3fafc6a9    -0.8580941645036908    -    ATGGACACAGA
    ============
    """
    filterChr = set(filterChr)
    if score_cutoff is None:
        score_cutoff = (-1.5, 2.5)
    infile, lines = open_file_gz_or_txt(infileName)

    if save_unified_format:
        outf = gzip.open(outfn, 'wt')
        outf.write(f"ID\tChr\tPos\tStrand\tScore\n")

    cpgDict = defaultdict(list)
    row_count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in tqdm(infile, total=lines, desc="Import-Tombo"):
        tmp = row.strip().split("\t")

        if tmp[chr_col] not in filterChr:
            continue

        if output_first:
            logger.debug(f'row = {list(enumerate(tmp))}')
            output_first = False

        if baseFormat == 1:
            try:
                start = int(tmp[start_col]) + 1
                strand = tmp[strand_col]
                if strand == '-':
                    start = start + 1
            except:
                logger.error(f" ####Tombo parse error at row={row}")
                continue
        elif baseFormat == 0:
            try:
                start = int(tmp[start_col])
                strand = tmp[strand_col]

                if strand == '-':
                    start = start + 1
            except Exception as e:
                logger.error(f" ####Tombo parse error at row={row}, exception={e}")
                continue
        else:
            logger.error(
                f"###\timportPredictions_Tombo InputValueError: baseCount value set to '{baseFormat}'. It should be equal to 0 or 1")
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')

        try:
            methCallTombo = float(tmp[meth_col])
        except Exception as e:
            logger.error(f" ####Tombo parse error at row={row}, exception={e}")
            continue

        meth_score = -methCallTombo

        if save_unified_format:
            # output to 1-based for meteore, ref: https://github.com/comprna/METEORE/blob/master/script_in_snakemake/extract_tombo_per_read_results.py
            outf.write(f"{tmp[readid_col]}\t{tmp[chr_col]}\t{start}\t{tmp[strand_col]}\t{methCallTombo}\n")

        key = (tmp[chr_col], start, strand)
        if methCallTombo < score_cutoff[0]:  # below -1.5 is methylated by default
            meth_indicator = 1
            meth_cnt += 1
        elif methCallTombo > score_cutoff[1]:  # above 2.5 is methylated by default
            meth_indicator = 0
            unmeth_cnt += 1
        else:
            continue

        if include_score:
            cpgDict[key].append((meth_indicator, meth_score))
        else:
            cpgDict[key].append(meth_indicator)
        row_count += 1

    infile.close()

    if save_unified_format:
        outf.close()
        logger.debug(f'Save METEORE output format to {outfn}')

    logger.debug(
        f"###\timportPredictions_Tombo SUCCESS: {row_count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-call={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs with score_cutoff={score_cutoff} from {infileName} file")
    return cpgDict


def importPredictions_DeepMod_C(infileName, chr_col=0, start_col=1, strand_col=5, coverage_col=-3, meth_freq_col=-2,
                                meth_cov_col=-1, baseFormat=1, output_first=False, include_score=False,
                                siteLevel=False, filterChr=HUMAN_CHR_SET, total_cols=12):
    """
    DeepMod RNN results format
    We treate input as 0-based format for start col.
    Return dict of key='chr1  123  +', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.
    Note that DeepMod only generate genome-level stats, we use meth cov and total coverage columns for read level evaluation.

    DeepMod BED format ref: https://github.com/WGLab/DeepMod/blob/master/docs/Results_explanation.md#2-format-of-output
    The output is in a BED format like below. The first six columns are Chr, Start pos, End pos, Base, Capped coverage, and Strand, and the last three columns are Real coverage, Mehylation percentage and Methylation coverage.

    Description (https://github.com/WGLab/DeepMod/blob/master/docs/Results_explanation.md):
    The output is in a BED format like below. The first six columns are Chr,
    Start pos, End pos, Base, Capped coverage, and Strand, and the last three
    columns are Real coverage, Mehylation percentage and Methylation coverage.

   The output is in a BED format like below. The first six columns are Chr, Start pos, End pos, Base, Capped coverage, and Strand, and the last three columns are Real coverage, Mehylation percentage and Methylation coverage.
    chr6 148655 148656 C 10 -  148655 148656 0,0,0 10 10 1
    chr6 148657 148658 C 12 +  148657 148658 0,0,0 12 8 1
    chr6 148674 148675 C 14 -  148674 148675 0,0,0 14 7 1
    chr6 148675 148676 C 15 -  148675 148676 0,0,0 15 6 1
    chr6 148676 148677 C 14 -  148676 148677 0,0,0 14 7 1
    chr6 148684 148685 C 12 -  148684 148685 0,0,0 12 25 3
    chr6 148685 148686 C 16 -  148685 148686 0,0,0 16 6 1
    chr6 148689 148690 C 11 +  148689 148690 0,0,0 11 72 8
    chr6 148691 148692 C 10 +  148691 148692 0,0,0 10 50 5
    chr6 148693 148694 C 8 +  148693 148694 0,0,0 8 100 8
    chr6 148694 148695 C 11 -  148694 148695 0,0,0 11 54 6
    chr6 148695 148696 C 10 +  148695 148696 0,0,0 10 90 9
    chr6 148697 148698 C 12 +  148697 148698 0,0,0 12 50 6
    chr6 148699 148700 C 9 +  148699 148700 0,0,0 9 22 2
    chr6 148701 148702 C 13 -  148701 148702 0,0,0 13 7 1
    chr6 148703 148704 C 13 -  148703 148704 0,0,0 13 15 2
    chr6 148706 148707 C 9 -  148706 148707 0,0,0 9 22 2
    Note: it is space-separated in original result file, not tab-separated file
    ============
    """
    filterChr = set(filterChr)
    infile, lines = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict(list)
    count_calls = 0
    meth_cnt = 0
    unmeth_cnt = 0

    first_row_checked = False

    for row in tqdm(infile, total=lines, desc="Import-DeepMod.C"):
        tmp = row.strip().split()

        if len(tmp) != total_cols:
            raise Exception(f"DeepMod RNN output format error: tmp={tmp}")

        if tmp[chr_col] not in filterChr:
            continue

        if output_first:
            logger.debug(f'row = {list(enumerate(tmp))}')
            output_first = False

        if not first_row_checked:
            # We check if this is DeepMod.C format, not DeepMod.Cluster Mode, ensure correct file used
            cov_num = int(tmp[coverage_col])
            meth_num = int(tmp[meth_cov_col])
            meth_freq_num = int(tmp[meth_freq_col])
            if int(100.0 * meth_num / cov_num) != meth_freq_num:
                raise Exception(f'Found not DeepMod.C format in file {infileName}, with line: {row}')
            first_row_checked = True

        if baseFormat == 1:
            start = int(tmp[start_col]) + 1
            strand = tmp[strand_col]
        elif baseFormat == 0:
            start = int(tmp[start_col])
            strand = tmp[strand_col]
        else:
            logger.debug(
                "###\timportPredictions_DeepMod InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(
                    baseFormat))
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')

        key = (tmp[chr_col], start, strand)

        methReads = int(tmp[meth_cov_col])
        coverage = int(tmp[coverage_col])

        meth_freq = methReads / coverage

        if siteLevel:  ## Site level return
            methCallsList = (methReads / coverage, coverage)
        elif include_score:  ## Read level return with score
            methCallsList = [(1, 1.0)] * methReads + [(0, 0.0)] * (coverage - methReads)
        else:  ## Read level return only class predictions
            methCallsList = [1] * methReads + [0] * (coverage - methReads)

        if key in cpgDict:
            raise Exception(f'In DeepMod.C results, we found duplicate key={key}, this is not correct')

        cpgDict[key] = methCallsList

        count_calls += len(methCallsList)
        meth_cnt += methReads
        unmeth_cnt += coverage - methReads

    infile.close()

    logger.debug(
        f"###\timportPredictions_DeepMod SUCCESS: {count_calls:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")
    return cpgDict


def importPredictions_DeepMod_Clustered(infileName, chr_col=0, start_col=1, strand_col=5, coverage_col=-4,
                                        meth_cov_col=-2, clustered_meth_freq_col=-1, baseFormat=1,
                                        output_first=False, siteLevel=True, include_score=False,
                                        filterChr=HUMAN_CHR_SET,
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
    filterChr = set(filterChr)
    infile, lines = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict()
    count_calls = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in tqdm(infile, total=lines, desc="Import-DeepMod.Cluster"):
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

        if key in cpgDict:
            raise Exception(f'In DeepMod_Cluster results, we found duplicate key={key}, this is not correct')

        if siteLevel:  # in dict
            cpgDict[key] = (meth_freq, coverage)  # {'freq': meth_freq, 'cov': coverage}
        elif include_score:  # For read-level include scores
            cpgDict[key] = [(1, 1.0)] * meth_cov + [(0, 0.0)] * (coverage - meth_cov)
        else:  # For read-level no scores
            cpgDict[key] = [1] * meth_cov + [0] * (coverage - meth_cov)

        count_calls += coverage
        meth_cnt += meth_cov
        unmeth_cnt += coverage - meth_cov

    infile.close()

    logger.debug(
        f"###\timportDeepMod_clustered Parsing SUCCESS: {count_calls:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")

    return cpgDict


def importPredictions_DeepMod(infileName, chr_col=0, start_col=1, strand_col=5, coverage_col=-4, meth_freq_col=-3,
                              meth_cov_col=-2, clustered_meth_freq_col=-1, baseFormat=1,
                              output_first=False, siteLevel=True, include_score=False, filterChr=HUMAN_CHR_SET):
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
    filterChr = set(filterChr)
    infile, lines = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict()
    count_calls = 0
    meth_cnt = 0
    unmeth_cnt = 0
    deepmod_format = None
    for row in tqdm(infile, total=lines, desc="Import-DeepMod"):
        tmp = row.strip().split()
        if deepmod_format is None:
            if len(tmp) == 12:
                deepmod_format = "DeepMod.C"
                # redefine DeepMod.C column
                coverage_col = -3
                meth_freq_col = -2
                meth_cov_col = -1
                clustered_meth_freq_col = None
            elif len(tmp) == 13:
                deepmod_format = "DeepMod.Cluster"
            else:
                raise Exception(f"Detect DeepMod output format error: tmp={tmp}")
            logger.debug(f"Detect DeepMod format is:{deepmod_format}")

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

        if clustered_meth_freq_col is not None:
            meth_freq = int(tmp[clustered_meth_freq_col]) / 100.0
        else:
            meth_freq = int(tmp[meth_freq_col]) / 100.0
        coverage = int(tmp[coverage_col])
        meth_cov = int(tmp[meth_cov_col])

        if key in cpgDict:
            raise Exception(f'In DeepMod_Cluster results, we found duplicate key={key}, this is not correct')

        if siteLevel:  # For site-level return
            cpgDict[key] = (meth_freq, coverage)
        elif include_score:  # For read-level include scores
            cpgDict[key] = [(1, 1.0)] * meth_cov + [(0, 0.0)] * (coverage - meth_cov)
        else:  # For read-level no scores
            cpgDict[key] = [1] * meth_cov + [0] * (coverage - meth_cov)

        count_calls += coverage
        meth_cnt += meth_cov
        unmeth_cnt += coverage - meth_cov

    infile.close()

    logger.debug(
        f"###\timportDeepMod_clustered Parsing SUCCESS: {count_calls:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")

    return cpgDict


def importPredictions_Megalodon(infileName, readid_col=0, chr_col=1, start_col=3, strand_col=2, mod_log_prob_col=4,
                                can_log_prob_col=5, meth_type_col=6, baseFormat=1, score_cutoff=(0.2, 0.8), sep='\t',
                                output_first=False,
                                include_score=False, filterChr=HUMAN_CHR_SET, save_unified_format=False, outfn=None):
    """
    0-based start for Magelodon：
        1.  baseFormat=0， start=Megalondon start；
        2.  baseFormat=1， start=Megalondon start +1
    Return dict of key='chr1 123 +', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.

    Note that the function requires per-read file, not per-site methlyation level.

    ### Example default input format from Megalondon (pre-processed to containing strand-info):
    head /fastscratch/c-panz/K562/megalodon/per_read_modified_base_calls.merged.sorted.bed
    chr1	10468	a0b5e934-3784-4163-b46c-f575ac1015bf	+	-1.6103846163370363	-0.2229070153110119	    m
    chr1	10470	a0b5e934-3784-4163-b46c-f575ac1015bf	+	-2.035311286540776	-0.13999775139928522	m
    chr1	10483	a0b5e934-3784-4163-b46c-f575ac1015bf	+	-1.5477196338381982	-0.2391872270542014 	m

    ** By default setting in Megalodon, for a modified base, mod_prob > 0.8; for a canonical base, can_prob > 0.8 (In other words, mod_prob < 0.2)
    For any position where no probability is greater than 0.8 neither canonical not modified bases get a count.
    (https://github.com/nanoporetech/megalodon/issues/47#issuecomment-673742805)
    ============
    """
    filterChr = set(filterChr)

    if score_cutoff is None:
        score_cutoff = (0.2, 0.8)
    infile, lines = open_file_gz_or_txt(infileName)
    if save_unified_format:
        outf = gzip.open(outfn, 'wt')
        outf.write(f"ID\tChr\tPos\tStrand\tScore\n")

    cpgDict = defaultdict(list)
    call_cnt = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in tqdm(infile, total=lines, desc="Import-Megalodon"):
        tmp = row.strip().split(sep)

        if tmp[chr_col] not in filterChr:
            continue

        if tmp[meth_type_col] != 'm':  # neglect 5hmc row
            continue

        if output_first:
            logger.debug(f'row = {list(enumerate(tmp))}')
            output_first = False

        if baseFormat == 1:
            try:
                start = int(tmp[start_col]) + 1
                strand = tmp[strand_col]
            except Exception as e:
                logger.error(f" ####Megalodon parse error at row={row}, exception={e}")
                continue
        elif baseFormat == 0:
            try:
                start = int(tmp[start_col])
                strand = tmp[strand_col]
            except Exception as e:
                logger.error(f" ####Megalodon parse error at row={row}, exception={e}")
                continue
        else:
            logger.error(
                f"###\timportPredictions_Megalodon_Read_Level InputValueError: baseFormat value set to '{baseFormat}'. It should be equal to 0 or 1")
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')

        if save_unified_format:
            # output to 1-based for meteore, ref: https://github.com/comprna/METEORE/blob/master/script/format_megalodon.R#L16
            meteore_score = float(tmp[mod_log_prob_col]) - float(tmp[can_log_prob_col])
            outf.write(f"{tmp[readid_col]}\t{tmp[chr_col]}\t{start}\t{tmp[strand_col]}\t{meteore_score}\n")

        key = (tmp[chr_col], start, strand)

        try:
            meth_prob = float(np.e ** float(tmp[mod_log_prob_col]))  # Calculate mod_prob
        except Exception as e:
            logger.error(f" ####Megalodon parse error at row={row}, exception={e}")
            continue

        # logger.debug(f"meth_prob={meth_prob:.4f}, unmeth_prob={float(np.e ** float(tmp[can_log_prob_col])):.4f}, score_cutoff={score_cutoff}")
        if meth_prob > score_cutoff[1]:  ##Keep methylated reads
            meth_indicator = 1
            meth_cnt += 1
        elif meth_prob < score_cutoff[0]:  ##Count unmethylated reads
            meth_indicator = 0
            unmeth_cnt += 1
        else:  ## Neglect other cases 0.2<= prob <=0.8
            continue

        if include_score:
            cpgDict[key].append((meth_indicator, meth_prob))
        else:
            cpgDict[key].append(meth_indicator)
        call_cnt += 1

    infile.close()

    if save_unified_format:
        outf.close()
        logger.debug(f'Save METEORE output format to {outfn}')

    logger.debug(
        f"###\timportPredictions_Megalodon SUCCESS: {call_cnt:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-call={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs with score_cutoff={score_cutoff} from {infileName} file")
    return cpgDict


def importPredictions_Guppy(infileName, baseFormat=1, sep='\t', output_first=False, include_score=False,
                            siteLevel=False, filterChr=HUMAN_CHR_SET, formatSource="raw"):
    """
    Import Guppy results by Fast5mod from ONT developed tools, it can report read-level or site-level results.
    Start is 0-based, and combined + and - CpGs together.
    :param infileName:
    :param sep:
    :param output_first:
    :param include_score:
    :param filterChr:
    :return:
    """
    filterChr = set(filterChr)
    if formatSource == 'raw':  # raw results of fast5mod
        chr_col = 0
        start_col = 1
        fw_meth_col = 3
        rv_meth_col = 4
        fw_unmeth_col = 5
        rv_unmeth_col = 6

        cpgDict = defaultdict(list)
        call_cnt = methcall_cnt = unmethcall_cnt = 0
        infile, lines = open_file_gz_or_txt(infileName)

        for row in tqdm(infile, total=lines, desc="Import-Guppy"):
            tmp = row.strip().split(sep)
            if tmp[chr_col] not in filterChr:
                continue

            try:
                start = int(tmp[start_col])
                fw_meth = int(tmp[fw_meth_col])
                rv_meth = int(tmp[rv_meth_col])
                fw_unmeth = int(tmp[fw_unmeth_col])
                rv_unmeth = int(tmp[rv_unmeth_col])

                fw_key = (tmp[chr_col], start + baseFormat, '+')
                rv_key = (tmp[chr_col], start + 1 + baseFormat, '-')

                if fw_meth + fw_unmeth > 0:  # + strand, cov>0, report it
                    coverage = fw_meth + fw_unmeth
                    methfreq = fw_meth / coverage

                    call_cnt += coverage
                    methcall_cnt += fw_meth
                    unmethcall_cnt += fw_unmeth

                    if siteLevel:  # in dict
                        cpgDict[fw_key] = (methfreq, coverage)  # {'freq': meth_freq, 'cov': coverage}
                    elif include_score:  # For read-level include scores
                        cpgDict[fw_key] = [(1, 1.0)] * fw_meth + [(0, 0.0)] * fw_unmeth
                    else:  # For read-level no scores
                        cpgDict[fw_key] = [1] * fw_meth + [0] * fw_unmeth
                if rv_meth + rv_unmeth > 0:  # - strand, cov>0, report it
                    coverage = rv_meth + rv_unmeth
                    methfreq = rv_meth / coverage

                    call_cnt += coverage
                    methcall_cnt += rv_meth
                    unmethcall_cnt += rv_unmeth

                    if siteLevel:  # in dict
                        cpgDict[rv_key] = (methfreq, coverage)  # {'freq': meth_freq, 'cov': coverage}
                    elif include_score:  # For read-level include scores
                        cpgDict[rv_key] = [(1, 1.0)] * rv_meth + [(0, 0.0)] * rv_unmeth
                    else:  # For read-level no scores
                        cpgDict[rv_key] = [1] * rv_meth + [0] * rv_unmeth
            except:
                logger.error(f"Parse Guppy fast5mod raw line failed: [{row}]")
                continue

        logger.debug(
            f"###\timportPredictions_Guppy (fast5mod) SUCCESS: {call_cnt:,} methylation calls (meth-calls={methcall_cnt:,}, unmeth-calls={unmethcall_cnt:,}) mapped to {len(cpgDict):,} CpGs, using default Fast5mod cutoff from {infileName} file")

        return cpgDict
    elif formatSource == 'raw-correct':  # raw results of fast5mod, using pandas
        the_dtype = {"chrom": str, "position": int, "motif": str,
                     "fwd.meth.count": int, "rev.meth.count": int,
                     "fwd.canon.count": int, "rev.canon.count": int}
        methdata = pd.read_csv(infileName, sep=sep, header=None, dtype=the_dtype, error_bad_lines=False,
                               warn_bad_lines=False, engine='c',
                               names=["chrom", "position", "motif",
                                      "fwd.meth.count", "rev.meth.count",
                                      "fwd.canon.count", "rev.canon.count"])

        # Forward strand
        fwd_methdata = methdata[["chrom", "position", "motif", "fwd.meth.count", "fwd.canon.count"]].copy()
        # Convert 0-based into 1-based
        fwd_methdata['position'] = fwd_methdata['position'] + baseFormat
        fwd_methdata['fwd.coverage'] = fwd_methdata['fwd.meth.count'] + fwd_methdata['fwd.canon.count']
        # Keep >0 coverage
        fwd_methdata = fwd_methdata[fwd_methdata['fwd.coverage'] > 0]
        fwd_methdata['fwd.methfreq'] = fwd_methdata['fwd.meth.count'] / fwd_methdata['fwd.coverage']
        fwd_methdata['strand'] = ['+'] * fwd_methdata.shape[0]

        # Reversed strand
        rev_methdata = methdata[["chrom", "position", "motif", "rev.meth.count", "rev.canon.count"]].copy()
        # Convert 0-based into 1-based, indicate G position
        rev_methdata['position'] = rev_methdata['position'] + 1 + baseFormat
        ## Keep >0 coverage
        rev_methdata['rev.coverage'] = rev_methdata['rev.meth.count'] + rev_methdata['rev.canon.count']
        rev_methdata = rev_methdata[rev_methdata['rev.coverage'] > 0]
        rev_methdata['rev.methfreq'] = rev_methdata['rev.meth.count'] / rev_methdata['rev.coverage']
        rev_methdata['strand'] = ['-'] * rev_methdata.shape[0]

        # Merge result
        fwd_methdata.columns = rev_methdata.columns = ["chrom", "position", "motif", "meth.count", "canon.count",
                                                       "coverage", "meth_freq", "strand"]
        df_combined = pd.concat([fwd_methdata, rev_methdata], ignore_index=True, axis=0)
        df_combined = df_combined[
            ["chrom", "position", "strand", "motif", "meth.count", "canon.count", "coverage", "meth_freq"]]
    else:  # our preprocessed results by ZW
        if baseFormat != 1:
            raise Exception(f"formatSource={formatSource} is not finished for base=0")
        df_combined = pd.read_csv(infileName, sep=sep, header=None,
                                  names=["chrom", "position", "strand", "motif", "meth.count", "canon.count",
                                         "coverage", "meth_freq"])
    logger.debug(df_combined)

    ## extract cpg into dict return object
    cpgDict = defaultdict()
    call_cnt = methcall_cnt = unmethcall_cnt = 0
    for index, row in df_combined.iterrows():
        chr = row['chrom']

        if chr not in filterChr:  # Filter out interested chrs
            continue
        start = int(row['position'])
        strand = row['strand']

        methfreq = float(row["meth_freq"])
        coverage = int(row["coverage"])

        meth_cov = int(row["meth.count"])
        unmeth_cov = int(row["canon.count"])

        call_cnt += coverage
        methcall_cnt += meth_cov
        unmethcall_cnt += unmeth_cov

        if strand not in ['+', '-']:
            raise Exception(
                f'strand={strand}  for row={row} is not acceptable, please check use correct function to parse Guppy output file {infileName}')

        key = (chr, int(start), strand)

        if siteLevel:  # in dict
            cpgDict[key] = (methfreq, coverage)  # {'freq': meth_freq, 'cov': coverage}
        elif include_score:  # For read-level include scores
            cpgDict[key] = [(1, 1.0)] * meth_cov + [(0, 0.0)] * unmeth_cov
        else:  # For read-level no scores
            cpgDict[key] = [1] * meth_cov + [0] * unmeth_cov

        if output_first:
            logger.debug(f'line={row}, key={key}, cpgDict[key]={cpgDict[key]}')
        output_first = False

    logger.debug(
        f"###\timportPredictions_Guppy SUCCESS: {call_cnt:,} methylation calls (meth-calls={methcall_cnt:,}, unmeth-calls={unmethcall_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")
    return cpgDict


def importPredictions_Guppy_gcf52ref(infileName, baseFormat=1, chr_col=0, strand_col=1, start_col=2, readid_col=4,
                                     log_ratio_col=5, log_lik_methylated_col=6, score_cutoff=(0.25, 0.5), header=None,
                                     sep='\t', output_first=False, include_score=False, filterChr=HUMAN_CHR_SET,
                                     save_unified_format=False, outfn=None):
    """
    Parse read level gcf52ref format, start position is 0-base, end is 1-base
    Sample file:
    #chromosome     strand  start   end     read_name       log_lik_ratio   log_lik_methylated      log_lik_unmethylated    num_calling_strands     num_motifs      sequence
    chr1    -       11343976        11343977        9d47a371-d2af-4104-8097-c5c159035f1e    0.847   -0.0561 -0.903  1       1
           CG
    chr1    -       11343987        11343988        9d47a371-d2af-4104-8097-c5c159035f1e    -1.8    -1.81   -0.00512        1
           1       CG
    chr1    -       11344167        11344168        9d47a371-d2af-4104-8097-c5c159035f1e    2.41    0.0     -2.41   1       1
           CG
    chr1    -       11344218        11344219        9d47a371-d2af-4104-8097-c5c159035f1e    1.8     -0.00512        -1.81   1
           1       CG
    :param infileName:
    :param baseFormat:
    :param chr_col:
    :param strand_col:
    :param start_col:
    :param log_lik_methylated_col:
    :param score_cutoff:
    :param sep:
    :param output_first:
    :param include_score:
    :param siteLevel:
    :param filterChr:
    :return:
    """
    filterChr = set(filterChr)
    ### Using scanner way
    if score_cutoff is None:
        score_cutoff = (0.25, 0.5)
    cpgDict = defaultdict(list)
    call_cnt = methcall_cnt = 0
    infile, lines = open_file_gz_or_txt(infileName)

    if save_unified_format:
        outf = gzip.open(outfn, 'wt')
        outf.write(f"ID\tChr\tPos\tStrand\tScore\n")

    for row in tqdm(infile, total=lines, desc="Import-Guppy.gcf52ref"):
        tmp = row.strip().split(sep)
        if tmp[chr_col] not in filterChr:
            continue
        chr = tmp[chr_col]
        strand = tmp[strand_col]
        start = int(tmp[start_col])
        log_lik_methylated = float(tmp[log_lik_methylated_col])
        log_ratio = float(tmp[log_ratio_col])
        readid = tmp[readid_col]
        if strand not in ['+', '-']:
            raise Exception(f"gcf52ref format strand parse error, row={row}")
        # Correct the coordinates into 1-based
        if strand == "+":
            start = start + baseFormat
        elif strand == "-":
            start = start + baseFormat + 1

        prob_methylated = 10 ** log_lik_methylated

        if save_unified_format:
            # output to 1-based for meteore, ref: https://github.com/comprna/METEORE/blob/master/script_in_snakemake/format_guppy.R
            outf.write(f"{readid}\t{chr}\t{start}\t{strand}\t{log_ratio}\n")
        if (prob_methylated > score_cutoff[0]) and (prob_methylated < score_cutoff[1]):
            continue

        if prob_methylated <= score_cutoff[0]:
            meth_indicator = 0
        else:
            meth_indicator = 1
        call_cnt += 1
        methcall_cnt += meth_indicator

        key = (chr, int(start), strand)

        if include_score:  # For read-level include scores
            cpgDict[key].append((meth_indicator, prob_methylated))
        else:  # For read-level no scores
            cpgDict[key].append(meth_indicator)
    infile.close()
    if save_unified_format:
        outf.close()
        logger.debug(f'Save METEORE output format to {outfn}')

    unmethcall_cnt = call_cnt - methcall_cnt
    logger.debug(
        f"###\timportPredictions_Guppy_gcf52ref SUCCESS: {call_cnt:,} methylation calls (meth-calls={methcall_cnt:,}, unmeth-calls={unmethcall_cnt:,}) mapped to {len(cpgDict):,} CpGs, using score_cutoff={score_cutoff} from {infileName} file")
    return cpgDict


def importPredictions_METEORE(infileName, readid_col=0, chr_col=1, start_col=2, meth_indicator_col=-3, meth_prob_col=-2,
                              strand_col=-1, baseFormat=1, include_score=False, filterChr=HUMAN_CHR_SET,
                              save_unified_format=False, outfn=None, toolname="METEORE"):
    """
    We checked input as 1-based format for start col.
    Return dict of key='chr1\t123\t123\t+', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.
    Note that the function requires per read stats, not frequencies of methylation.
    ### Example input format from DeepSignal:
    zcat HL60.METEORE.megalodon_deepsignal-default-model-perRead.combine.tsv.gz | head
    ID	Chr	Pos	Prediction	Prob_methylation
    ID	Chr	Pos	Prediction	Prob_methylation	Strand
    abcc33b7-2895-4e0e-830c-3d0ed441760d	chr1	103867	0	0.331208615558345	-
    abcc33b7-2895-4e0e-830c-3d0ed441760d	chr1	103984	0	0.331208615558345	-
    abcc33b7-2895-4e0e-830c-3d0ed441760d	chr1	103989	0	0.331208615558345	-
    abcc33b7-2895-4e0e-830c-3d0ed441760d	chr1	104048	0	0.331208615558345	-
    abcc33b7-2895-4e0e-830c-3d0ed441760d	chr1	104243	0	0.3042032634832675	-
    """
    filterChr = set(filterChr)

    infile, lines = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict(list)
    row_count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    if save_unified_format:
        outf = gzip.open(outfn, 'wt')
        outf.write(f"ID\tChr\tPos\tStrand\tScore\n")

    for row in tqdm(infile, total=lines, desc=f"Import-{toolname}"):
        if row.startswith("ID\tChr"):  # skim header
            continue
        tmp = row.strip().split("\t")

        if tmp[chr_col] not in filterChr:
            continue

        readid = tmp[readid_col]
        start_1_base = int(tmp[start_col])
        meth_indicator = int(tmp[meth_indicator_col])
        meth_prob = float(tmp[meth_prob_col])
        strand = tmp[strand_col]

        if save_unified_format:
            # output to 1-based readl-level format
            log_ratio_score = math.log2((meth_prob + EPSLONG) / (1 - meth_prob + EPSLONG))
            outf.write(f"{readid}\t{tmp[chr_col]}\t{start_1_base}\t{strand}\t{log_ratio_score}\n")

        ## Since METEORE report chr pos (1-based), strand
        ## pos is point to CG' C in +, and CG's G in -
        key = (tmp[chr_col], start_1_base - (1 - baseFormat), strand)
        if include_score:
            cpgDict[key].append((meth_indicator, meth_prob))
        else:
            cpgDict[key].append(meth_indicator)

        row_count += 1
        if meth_indicator == 1:
            meth_cnt += 1
        else:
            unmeth_cnt += 1
    infile.close()
    if save_unified_format:
        outf.close()
        logger.debug(f'Save {toolname} output format to {outfn}')
    logger.debug(
        f"###\timportPredictions_{toolname} SUCCESS: rows={row_count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs (include + and -) from {infileName} file")
    return cpgDict


def importPredictions_UNIREAD(infileName, readid_col=0, chr_col=1, start_col=2, strand_col=3, meth_score_col=4,
                              score_cutoff=(0, 0), sep='\t',
                              baseFormat=1, include_score=False, filterChr=HUMAN_CHR_SET,
                              save_unified_format=False, outfn=None, toolname=None):
    """
    Import the unified read level format

    For Nanopolish/Megalodon UNIREAD, may apply score_cutoff=(-2,2), ref prob=0.2,0.8
    For DeepSignal UNIREAD, apply score_cutoff=(0,0), ref prob=0.5

    Sample format:
    ID	Chr	Pos	Strand	Score
    a0b5e934-3784-4163-b46c-f575ac1015bf	chr1	10469	+	-1.3874776010260244
    a0b5e934-3784-4163-b46c-f575ac1015bf	chr1	10471	+	-1.895313535141491
    a0b5e934-3784-4163-b46c-f575ac1015bf	chr1	10484	+	-1.3085324067839967
    a0b5e934-3784-4163-b46c-f575ac1015bf	chr1	10489	+	4.122552265302514
    a0b5e934-3784-4163-b46c-f575ac1015bf	chr1	10493	+	4.4188115095774565
    """
    filterChr = set(filterChr)

    if score_cutoff is None:  # Default if provided None
        score_cutoff = (0, 0)
    infile, lines = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict(list)
    row_count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    if save_unified_format:
        outf = gzip.open(outfn, 'wt')
        outf.write(f"ID\tChr\tPos\tStrand\tScore\n")

    for row in tqdm(infile, total=lines, desc=f"Import-{toolname}"):
        if row.startswith("ID\tChr"):  # skim header
            continue
        tmp = row.strip().split(sep)

        chr = tmp[chr_col]
        if chr not in filterChr:
            continue
        readid = tmp[readid_col]
        start_1_base = int(tmp[start_col])
        strand = tmp[strand_col]
        meth_score = float(tmp[meth_score_col])

        if meth_score <= score_cutoff[0]:
            meth_indicator = 0
        elif meth_score > score_cutoff[1]:
            meth_indicator = 1
        else:
            continue

        if save_unified_format:
            # output to 1-based readl-level format
            outf.write(f"{readid}\t{chr}\t{start_1_base}\t{strand}\t{meth_score}\n")

        ## Since METEORE report chr pos (1-based), strand
        ## pos is point to CG' C in +, and CG's G in -
        key = (tmp[chr_col], start_1_base - (1 - baseFormat), strand)
        if include_score:
            cpgDict[key].append((meth_indicator, meth_score))
        else:
            cpgDict[key].append(meth_indicator)

        row_count += 1
        if meth_indicator == 1:
            meth_cnt += 1
        else:
            unmeth_cnt += 1
    infile.close()
    if save_unified_format:
        outf.close()
        logger.debug(f'Save {toolname} output format to {outfn}')
    logger.debug(
        f"###\timportPredictions_{toolname} SUCCESS: rows={row_count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs with score_cutoff={score_cutoff}, from {infileName} file")
    return cpgDict


def importPredictions_UNISITE(infileName, chr_col=0, start_col=1, strand_col=5,
                              meth_freq_col=6, meth_cov_col=7,
                              sep='\t', covCutoff=1, baseFormat=1,
                              include_score=False, siteLevel=True, includeCov=True,
                              filterChr=HUMAN_CHR_SET, toolname=None):
    """
    Import the unified read level format

    Sample format:
    chr10	10738	10739	.	.	+	0.87	1
    chr10	10750	10751	.	.	+	0.9	1
    chr10	10752	10753	.	.	+	0.9	1
    chr10	10767	10768	.	.	+	0.9	1
    chr10	10779	10780	.	.	+	0.9	1
    """
    filterChr = set(filterChr)

    infile, lines = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict()
    count_calls = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in tqdm(infile, total=lines, desc=f"Import-{toolname}"):
        if row.startswith("ID\tChr"):  # skim header
            continue
        tmp = row.strip().split(sep)

        chr = tmp[chr_col]
        if chr not in filterChr:
            continue
        start = int(tmp[start_col]) + baseFormat
        strand = tmp[strand_col]
        meth_freq = float(tmp[meth_freq_col])
        coverage = int(tmp[meth_cov_col])
        meth_cov = int(meth_freq * coverage)

        if coverage < covCutoff:
            continue

        ## Since METEORE report chr pos (1-based), strand
        ## pos is point to CG' C in +, and CG's G in -
        key = (tmp[chr_col], start, strand)
        if siteLevel:
            if includeCov:
                cpgDict[key] = (meth_freq, coverage)
            else:
                cpgDict[key] = meth_freq
        elif include_score:
            cpgDict[key] = [(1, 1.0)] * meth_cov + [(0, 0.0)] * (coverage - meth_cov)
        else:
            cpgDict[key] = [1] * meth_cov + [0] * (coverage - meth_cov)
        count_calls += coverage
        meth_cnt += meth_cov
        unmeth_cnt += coverage - meth_cov

    infile.close()
    logger.debug(
        f"###\timportPredictions_{toolname} SUCCESS: {count_calls:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs (include + and -) from {infileName} file")
    return cpgDict


def readLevelToSiteLevelWithCov(ontDict, minCov=1, toolname="Tool"):
    """
    Convert read level dict into site level dict.

    Filtering low coverage CpG sites in ONT call dict object, keep only sites with cov >= minCov. Convert to unified format of [[freq, cov], ......]
    1. input is key-value, value like [000011100] or like [{'freq':0.7, 'cov':8}, ...]
    2. output will be converted to [[freq, cov], ......]

    Convert orignial call object from dict[cpg] = {cpg: [1 0 11 1 1 0, ..., meth_freq_readn]} to dict[cpg] = {cpg:[meth_freq, coverage_number]}

    meth_freq   in [0.0,1.0]
    cov_num     in int
    Read-level -> Site-level

    :param ontDict:
    :param minCov:
    :param byLength: if False, will deal with DeepMod_cluster results
    :return:
    """
    result = defaultdict()
    for cpg in tqdm(ontDict, desc=f"read_to_site_{toolname}"):
        if type(ontDict[cpg]) == list:  # value is [0 0 1 1 0 ...]
            if len(ontDict[cpg]) >= minCov:
                result[cpg] = (sum(ontDict[cpg]) / float(len(ontDict[cpg])), len(ontDict[cpg]))
        elif type(ontDict[cpg]) == tuple:  # Used by DeepMod_cluster results, value is (freq, cov)
            if ontDict[cpg][1] >= minCov:  # no change for site level format, e.g., DeepModCluster
                result[cpg] = ontDict[cpg]
        else:
            raise Exception(
                f'Not support type of value, type(ontDict[cpg])={type(ontDict[cpg])} for toolname={toolname}')

    logger.debug(
        f"###\treadLevelToSiteLevelWithCov: completed filtering with minCov={minCov}: leave {len(result):,} CpG sites left for {toolname}")
    return result


def open_file_gz_or_txt(infn, return_lines=True):
    """
    Open a txt or gz file, based on suffix
    :param infn:
    :return:
    """
    logger.debug(f"open file: {infn}")

    if infn.endswith('.gz'):
        infile = gzip.open(infn, 'rt')  # using rt option, no need to convert bytearray
        command = f"zcat {infn} | wc -l"
    else:
        infile = open(infn, 'r')
        command = f"wc -l  {infn}"

    if return_lines:
        ret = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")
        num_lines = int(ret.strip())
        return (infile, num_lines)

    return infile


# encode format
def importGroundTruth_from_Encode(infileName, chr_col=0, start_col=1, methfreq_col=-1, cov_col=4, strand_col=5,
                                  covCutoff=1, baseFormat=1, filterChr=HUMAN_CHR_SET, includeCov=True):
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
    filterChr = set(filterChr)
    cpgDict = {}

    infile, lines = open_file_gz_or_txt(infileName)

    nrow = 0
    for row in tqdm(infile, total=lines, desc="Import-Encode"):
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

        if key in cpgDict:
            raise Exception(f"Found duplicate CpG sites for {row}")

        if includeCov:
            cpgDict[key] = (meth_freq, cov_cnt)
        else:
            cpgDict[key] = meth_freq

    infile.close()
    logger.debug(
        f"###\timportGroundTruth_BedMethyl_from_Encode: loaded information for {len(cpgDict):,} CpGs, with cutoff={covCutoff} ({nrow:,} rows)")

    return cpgDict


# bismark format, Deal with file name like 'bismark_bt2.CpG_report.txt.gz'
def importGroundTruth_from_Bismark(infn, chr_col=0, start_col=1, strand_col=2, meth_col=3, unmeth_col=4, covCutoff=1,
                                   baseFormat=1, filterChr=HUMAN_CHR_SET, includeCov=True, sep='\t', print_first=False):
    """
    We checked the input file is start using 1-based format.
    We use this format than others in Bismark due to it contains strand info.
    Ensure that in + strand, position is pointed to CG's C
                in - strand, position is pointed to CG's G, (for positive strand sequence)
            in our imported into programs.
    We design this due to other bismark output contains no strand-info.

    The genome-wide cytosine methylation output file is tab-delimited in the following format:
    ==========================================================================================
    <chromosome>  <position>  <strand>  <count methylated>  <count non-methylated>  <C-context>  <trinucleotide context>

    Sample input file:
    gunzip -cd HL60_RRBS_ENCFF000MDA.Read_R1.Rep_1_trimmed_bismark_bt2.CpG_report.txt.gz| head
    chr4	10164	+	0	0	CG	CGC
    chr4	10165	-	0	0	CG	CGT
    chr4	10207	+	0	0	CG	CGG
    chr4	10208	-	0	0	CG	CGA
    chr4	10233	+	0	0	CG	CGT
    chr4	10234	-	0	0	CG	CGG
    chr4	10279	+	0	0	CG	CGT
    chr4	10280	-	0	0	CG	CGT
    chr4	10297	+	0	0	CG	CGC
    chr4	10298	-	0	0	CG	CGC

    :param infn:
    :param chr_col:
    :param start_col:
    :param strand_col:
    :param meth_col:
    :param unmeth_col:
    :param ccontect_col:
    :param covCutoff:
    :param baseFormat:
    :param includeCov:
    :return:
    """
    filterChr = set(filterChr)
    infile, lines = open_file_gz_or_txt(infn)

    cpgDict = {}
    row_cnt = 0
    for row in tqdm(infile, total=lines, desc="Import-Bismark"):
        tmp = row.strip().split(sep)

        chr = tmp[chr_col]
        start = int(tmp[start_col])
        strand = tmp[strand_col]
        meth_cov = int(tmp[meth_col])
        unmeth_cov = int(tmp[unmeth_col])
        coverage = meth_cov + unmeth_cov

        if coverage < 1:  # Filter out all 0 coverage
            continue
        row_cnt += 1
        meth_freq = meth_cov / coverage

        if strand not in ['+', '-']:
            raise Exception(f"Error strand={strand}, for row={row}")

        if filterChr is not None and chr not in filterChr:
            continue

        if coverage < covCutoff:
            continue

        key = (chr, start - (1 - baseFormat), strand)
        if key in cpgDict:
            raise Exception(f"Found duplicate key, in row={row}, key={key}, from File={infn}")

        if includeCov:
            cpgDict[key] = (meth_freq, coverage)
        else:
            cpgDict[key] = meth_freq

        if print_first:
            logger.debug(f"key={key}, cpgDict[key]={cpgDict[key]}")
            print_first = False

    infile.close()
    logger.debug(
        f"###\timportGroundTruth_from_Bismark: loaded information for {len(cpgDict):,} CpGs with cutoff={covCutoff}, before cutoff={row_cnt:,}")
    return cpgDict

    df = pd.read_csv(infn, sep='\t', compression='gzip', header=None)

    df['cov'] = df.iloc[:, meth_col] + df.iloc[:, unmeth_col]
    df = df[df['cov'] >= covCutoff]

    # df = df[df.iloc[:, ccontect_col] == 'CG']

    # Based on import format 0, we need minus 1, due to input format is 1-based start
    if baseFormat == 0:
        df.iloc[:, start_col] = df.iloc[:, start_col] - 1
        df['end'] = df.iloc[:, start_col] + 1
    elif baseFormat == 1:
        df['end'] = df.iloc[:, start_col]

    df['meth-freq'] = (df.iloc[:, meth_col] / df['cov'] * 100).astype(np.int32)

    df = df.iloc[:,
         [chr_col, start_col, df.columns.get_loc("end"), df.columns.get_loc("meth-freq"), df.columns.get_loc("cov"),
          strand_col]]

    df.columns = ['chr', 'start', 'end', 'meth-freq', 'cov', 'strand']

    cpgDict = defaultdict(list)
    for index, row in df.iterrows():
        chr = row['chr']

        if chr not in HUMAN_CHR_SET:  # Filter out non-human chrs
            continue

        start = int(row['start'])
        strand = row['strand']

        if strand not in ['+', '-']:
            raise Exception(
                f'strand={strand}  for row={row} is not acceptable, please check use correct function to parse bgtruth file {infn}')

        key = (chr, int(start), strand)
        if key not in cpgDict:
            if includeCov:
                cpgDict[key] = (row['meth-freq'] / 100.0, row['cov'])
            else:
                cpgDict[key] = row['meth-freq'] / 100.0
        else:
            raise Exception(
                f'In genome-wide, we found duplicate sites: for key={key} in row={row}, please check input file {infn}')

    logger.debug(
        f"###\timportGroundTruth_genome_wide_from_Bismark: loaded {len(cpgDict):,} CpGs with cutoff={covCutoff} from file {infn}")

    return cpgDict


def importGroundTruth_from_5hmc_ziwei(infn, chr_col=0, start_col=1, strand_col=7, meth_freq_col=3,
                                      meth_5hmc_freq_col=4, sep='\t', covCutoff=1, baseFormat=1, includeCov=True,
                                      filterChr=HUMAN_CHR_SET, cov_col=(8, 9), provided_cov=5, print_first=False):
    """
    Input format is below, start is 1-based:
    chr1	10542	10543	1	0	0	0	+	2	1
    chr1	10563	10564	0.75	0	0.25	0	+	3	1
    chr1	10571	10572	1	0	0	0	+	3	1
    chr1	10577	10578	1	0	0	0	+	3	2
    chr1	10579	10580	1	0	0	0	+	3	2
    chr1	10589	10590	1	0	0	0	+	3	2
    chr1	10609	10610	1	0	0	0	+	3	2
    chr1	10617	10618	1	0	0	0	+	3	2
    chr1	10620	10621	1	0	0	0	+	2	2
    chr1	10631	10632	1	0	0	0	+	2	2

    chr1	10542	10543	1	0	0	0	+
    chr1	10563	10564	0.75	0	0.25	0	+
    chr1	10571	10572	1	0	0	0	+
    chr1	10577	10578	1	0	0	0	+
    chr1	10579	10580	1	0	0	0	+
    chr1	10589	10590	1	0	0	0	+
    chr1	10609	10610	1	0	0	0	+
    chr1	10617	10618	1	0	0	0	+
    chr1	10620	10621	1	0	0	0	+
    chr1	10631	10632	1	0	0	0	+

    The columns are chromosome name, start position, end position, 5-mC level, 5-hmC level, unmethylated level and number of conflicts.
    """
    filterChr = set(filterChr)
    infile = open_file_gz_or_txt(infn)

    cpgDict = {}
    row_cnt = 0
    for row in infile:
        tmp = row.strip().split(sep)
        row_cnt += 1

        chr = tmp[chr_col]
        start = int(tmp[start_col])
        strand = tmp[strand_col]
        meth_freq = float(tmp[meth_freq_col])

        if strand not in ['+', '-']:
            raise Exception(f"Error strand={strand}, for row={row}")

        if filterChr is not None and chr not in filterChr:
            continue

        if covCutoff is not None and cov_col is not None:
            # max cov of bs and oxbs
            coverage = max([int(tmp[cov_col_i]) for cov_col_i in cov_col])
            if coverage < covCutoff:
                continue

        key = (chr, start - (1 - baseFormat), strand)
        if key in cpgDict:
            raise Exception(f"Found duplicate key, in row={row}, key={key}")

        if includeCov:
            if cov_col is None:
                cpgDict[key] = (meth_freq, provided_cov)
            else:
                # max cov of bs and oxbs
                coverage = max([int(tmp[cov_col_i]) for cov_col_i in cov_col])
                cpgDict[key] = (meth_freq, coverage)
        else:
            cpgDict[key] = meth_freq

        if print_first:
            logger.debug(f"key={key}, cpgDict[key]={cpgDict[key]}")
            print_first = False

    infile.close()
    logger.debug(
        f"###\timportGroundTruth_5hmc_ziwei: loaded information for {len(cpgDict):,} CpGs with cutoff={covCutoff}, before cutoff={row_cnt:,}")
    return cpgDict


def import_call(infn, encode, baseFormat=1, include_score=False, siteLevel=False, enable_cache=False, using_cache=False,
                filterChr=HUMAN_CHR_SET, save_unified_format=False, outfn=None, cache_dir=None, toolname=None,
                score_cutoff=None):
    """
    General purpose for import any tools methylation calling input files.

    Import fn based on callname and return key=(chr, start, strand) -> value=[0 0 1 1 1 1 ]

    If include_score=True, then value= (0, score)
        call0   -   original results
    :param infn:
    :param encode:
    :return:
    """
    logger.debug("\n\n####################\n\n")
    logger.debug(f"Start load encode={encode}, infn={infn}")

    if enable_cache and using_cache:
        ret = check_cache_available(infn=infn, encode=encode, baseFormat=baseFormat, include_score=include_score,
                                    siteLevel=siteLevel, cache_dir=cache_dir, file_type='ont-call')
        if ret is not None:
            logger.debug(f'Import {encode} finished from cache, CpGs={len(ret):,}\n')
            return ret
        logger.debug(f'Not cached yet, we load from raw file')
    call0 = None
    if encode == 'DeepSignal':
        calls0 = importPredictions_DeepSignal(infn, baseFormat=baseFormat, include_score=include_score,
                                              score_cutoff=score_cutoff,
                                              filterChr=filterChr, save_unified_format=save_unified_format, outfn=outfn)
    elif encode == 'Tombo':
        calls0 = importPredictions_Tombo(infn, baseFormat=baseFormat, include_score=include_score, filterChr=filterChr,
                                         score_cutoff=score_cutoff,
                                         save_unified_format=save_unified_format, outfn=outfn)
    elif encode == 'Nanopolish':
        calls0 = importPredictions_Nanopolish(infn, baseFormat=baseFormat, include_score=include_score,
                                              score_cutoff=score_cutoff,
                                              filterChr=filterChr, save_unified_format=save_unified_format, outfn=outfn)
    elif encode == 'Megalodon':  # Original format
        calls0 = importPredictions_Megalodon(infn, baseFormat=baseFormat, include_score=include_score,
                                             score_cutoff=score_cutoff,
                                             filterChr=filterChr, save_unified_format=save_unified_format, outfn=outfn)
    elif encode == 'Megalodon.ZW':  # Here, ziwei preprocess the raw file to this format: chr_col=0, start_col=1, strand_col=3
        calls0 = importPredictions_Megalodon(infn, baseFormat=baseFormat, include_score=include_score,
                                             score_cutoff=score_cutoff,
                                             readid_col=2, chr_col=0, start_col=1, strand_col=3, filterChr=filterChr,
                                             save_unified_format=save_unified_format, outfn=outfn)
    elif encode == 'DeepMod':  # import DeepMod, support both C and Cluster input format
        calls0 = importPredictions_DeepMod(infn, baseFormat=baseFormat, siteLevel=siteLevel,
                                           include_score=include_score, filterChr=filterChr)
    elif encode == 'DeepMod.Cluster':  # import DeepMod Clustered output for site level
        calls0 = importPredictions_DeepMod_Clustered(infn, baseFormat=baseFormat, siteLevel=siteLevel,
                                                     include_score=include_score, filterChr=filterChr)
    elif encode == 'DeepMod.C':  # import DeepMod itself tool for read level
        calls0 = importPredictions_DeepMod_C(infn, baseFormat=baseFormat, siteLevel=siteLevel,
                                             include_score=include_score, filterChr=filterChr)
    elif encode == 'Guppy':  # import Guppy itself tool results by fast5mod raw results
        calls0 = importPredictions_Guppy(infn, baseFormat=baseFormat, siteLevel=siteLevel,
                                         include_score=include_score, filterChr=filterChr)
    elif encode == 'Guppy.ZW':  # import Guppy itself tool results by fast5mod processed by Ziwei
        calls0 = importPredictions_Guppy(infn, baseFormat=baseFormat, siteLevel=siteLevel,
                                         include_score=include_score, filterChr=filterChr, formatSource="Guppy.ZW")
    elif encode == 'Guppy.gcf52ref':  # import Guppy gcf52ref read level results
        if siteLevel:
            raise Exception(f"Not support site-level import for gcf52ref format file={infn}")
        calls0 = importPredictions_Guppy_gcf52ref(infn, baseFormat=baseFormat, include_score=include_score,
                                                  filterChr=filterChr, header=None, score_cutoff=score_cutoff,
                                                  save_unified_format=save_unified_format,
                                                  outfn=outfn)
    elif encode == 'METEORE':
        calls0 = importPredictions_METEORE(infn,
                                           baseFormat=baseFormat, include_score=include_score,
                                           filterChr=filterChr, save_unified_format=save_unified_format, outfn=outfn)
    elif encode == 'NANOME':
        calls0 = importPredictions_METEORE(infn, strand_col=3, meth_indicator_col=-2, meth_prob_col=-1,
                                           baseFormat=baseFormat, include_score=include_score,
                                           filterChr=filterChr, save_unified_format=save_unified_format, outfn=outfn,
                                           toolname="NANOME")
    elif encode == 'UNIREAD':
        calls0 = importPredictions_UNIREAD(infn,
                                           baseFormat=baseFormat, include_score=include_score,
                                           filterChr=filterChr, save_unified_format=save_unified_format, outfn=outfn,
                                           toolname=toolname, score_cutoff=score_cutoff)
    elif encode == 'UNISITE':
        calls0 = importPredictions_UNISITE(
            infn, baseFormat=baseFormat, siteLevel=siteLevel, include_score=include_score,
            filterChr=filterChr, toolname=toolname)
    else:
        raise Exception(f'Not support {encode} for file {infn} now')

    if enable_cache:
        save_to_cache(infn, calls0, encode=encode, baseFormat=baseFormat, include_score=include_score,
                      siteLevel=siteLevel, cache_dir=cache_dir, file_type='ont-call')

    logger.debug(f'Import {encode} finished!\n')
    return calls0


def import_bgtruth(infn, encode, covCutoff=1, baseFormat=1, includeCov=True, filterChr=HUMAN_CHR_SET,
                   enable_cache=False,
                   using_cache=False,
                   cache_dir=None):
    """
    General purpose for import BG-Truth input files.
    Import bgtruth from file fn using encode, when use new dataset, MUST check input file start baseFormat and import functions are consistent!!!
    :param infn:
    :param encode:
    :param includeCov:  if true return (freq, cov) as value of each key=(chr, (int)start, strand), or just value=freq. Note: freq is in range of [0,1]
    :return:
    """
    if enable_cache and using_cache:
        ret = check_cache_available(infn, encode=encode, cov=covCutoff, baseFormat=baseFormat, includeCov=includeCov,
                                    cache_dir=cache_dir, file_type='bs-seq')
        if ret is not None:
            logger.debug(f'Import BG-Truth using encode={encode} finished from cache, CpGs={len(ret):,}\n')
            return ret
        logger.debug(f'Not cached yet, we load from raw file')

    if encode == "encode":  # Encode BS-seq data, 0-based start
        bgTruth = importGroundTruth_from_Encode(infn, covCutoff=covCutoff, filterChr=filterChr,
                                                baseFormat=baseFormat, includeCov=includeCov)
    elif encode == "bismark":  # for genome-wide Bismark Report ouput format, 1-based start input
        bgTruth = importGroundTruth_from_Bismark(infn, covCutoff=covCutoff, filterChr=filterChr,
                                                 baseFormat=baseFormat, includeCov=includeCov)
    elif encode == 'UNISITE':
        bgTruth = importPredictions_UNISITE(infn, covCutoff=covCutoff, filterChr=filterChr, siteLevel=True,
                                            baseFormat=baseFormat, includeCov=includeCov)
    elif encode == "5hmc_ziwei":  # for sanity check 5hmc
        bgTruth = importGroundTruth_from_5hmc_ziwei(infn, covCutoff=covCutoff,
                                                    baseFormat=baseFormat, includeCov=includeCov)
    else:
        raise Exception(f"encode={encode} is not supported yet, for inputfile={infn}")

    if enable_cache:
        save_to_cache(infn, bgTruth, encode=encode, cov=covCutoff, baseFormat=baseFormat, includeCov=includeCov,
                      cache_dir=cache_dir, file_type='bs-seq')

    logger.debug(f'Import BG-Truth using encode={encode} finished!\n')
    return bgTruth


def calldict2bed(inputDict, null_str='.'):
    """
    Directly convert dict into bed object, assume key is (chr1  123   +),
    and output is bed like chr1  123  123 . .  +
    Note columns name:  chrom  start  end name score strand
    :param inputDict:
    :return:
    """
    dataset = defaultdict(list)
    for key in inputDict:
        dataset['chrom'].append(key[0])
        dataset['start'].append(int(key[1]))
        dataset['end'].append(int(key[1]))
        dataset['name'].append(null_str)
        dataset['score'].append(null_str)
        dataset['strand'].append(key[2])
    df = pd.DataFrame.from_dict(dataset)[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    bedret = BedTool.from_dataframe(df).sort()
    return bedret


def calldict2txt(inputDict):
    """
    Convert all keys in dict to a string txt, txt file is:
    chr1  123  123 . .  +

    :param inputDict:
    :return:
    """
    text = ""
    for key in inputDict:
        text += f'{key[0]}\t{key[1]}\t{key[1]}\t.\t.\t{key[2]}\n'
    return text


def bedtxt2dict(pybed, strand_col=5):
    """
    convert bed object to a dict with keys in bed file

    From
    chr123  123  123 +
    To
    ('chr123', 123, '+') -> 1

    :param pybed: bed object
    :return:
    """
    ret_dict = {}
    for t in pybed:
        ret = str(t).strip().split('\t')
        key = (ret[0], int(ret[1]), ret[strand_col])
        ret_dict[key] = 1
    return ret_dict


def load_single_sites_bed_as_set(infn):
    """
    Load the single sites BED files, and return a Set object represents sites.
    :param infn:
    :return:
    """
    if infn.endswith('.gz'):
        infile = gzip.open(infn, 'rt')
    else:
        infile = open(infn, 'r')
    ret = set()
    for row in infile:  # each row: chr123  123   123  .  .  +
        rowsplit = row[:-1].split('\t')
        key = (rowsplit[0], int(rowsplit[1]), rowsplit[5])
        ret.add(key)
    infile.close()
    return ret


def computePerReadPerfStats(ontCalls, bgTruth, title, coordBedFileName=None, secondFilterBedFileName=None,
                            cutoff_fully_meth=1.0, outdir=None, prefix_name=None, save_curve_data=True):
    """
    Compute ontCalls with bgTruth performance results by per-read count.
    coordBedFileName        -   full file name/ bed tuple of coordinate used to eval, this is the really region
    secondFilterBed         -   joined sets of four tools with bg-truth, or None with out joined sets

    Note: in our experiments, bed files of singleton, non-singleton and related files are all start 1-based format.

    90% methlation level if cutoff_meth=0.9.

    bedFile is coordinate: False(GenomeWide), Singleton, etc.
    secondFilterBed is second filter coordinate: False(No filter), Joined file name (four tools joined with BGTruth)

    secondFilterBed_4Corr is all tools at least four calls results joined, currently very confused results

    cutoff_meth is the percentage used in evaluation, 1-fully methylated  0.9- >=90% methylated used

    ontCalls - dictionary of CpG coordinates with their per read methylation call (1 or 0) // Format: {"chr\tstart\tend\n" : [list of methylation calls]}
    bsReference - dictionary of CpG coordinates with their methylation frequencies (range 0 - 1). This list is already prefiltered to meet minimal coverage (e.g. 4x) at this point. // Format: {"chr\tstart\tend\n" : methylation level (format: float (0-1))}
    title - prefix of the analysis, output plots etc. - should be as short as possible, but unique in context of other analyses
    bedFile - BED file which will be used to narrow down the list of CpGs for example to those inside CGIs or promoters etc.. By default "False" - which means no restrictions are done (i.e. genome wide)
    secondFilterBed - these should be CpGs covered in some reference list. Format: BED
    """
    # Number of returned tuple size
    num_tuples_ret = 21

    # Firstly reduce ontCalls with in bgTruth keys
    ontCallsKeySet = set(ontCalls.keys()).intersection(set(bgTruth.keys()))

    ontCalls_narrow_set = None  # Intersection of ontCall with coord, or None if genome-wide
    if coordBedFileName[1] != genome_wide_tagname:
        # Try ontCall intersect with coord (Genomewide, Singletons, etc.)
        # ontCalls_bed = BedTool(calldict2txt(ontCallsKeySet), from_string=True).sort()
        ontCalls_bed = calldict2bed(ontCallsKeySet)

        if isinstance(coordBedFileName, str):  # bed file path
            coordBed = get_region_bed(coordBedFileName)
            bedfn_basename = os.path.basename(coordBedFileName)
        elif isinstance(coordBedFileName, tuple):  # directly get the results
            coordBed = coordBedFileName[2]
            bedfn_basename = os.path.basename(coordBedFileName[0])
        else:
            raise Exception(f"not correct coordBedFileName={coordBedFileName}")

        ontCalls_intersect = intersect_bed_regions(ontCalls_bed, coordBed, bedfn_basename)
        ontCalls_narrow_set = set(bedtxt2dict(ontCalls_intersect).keys())
        tagname = get_region_tagname(bedfn_basename)
    else:  # for genome-wide case
        # bedfn_basename = 'x.x.Genome-wide'
        # region_name = 'Genome-wide'
        bedfn_basename = 'x.x.Genome-wide'  # used for intersect bed regions, strand identify
        tagname = coordBedFileName[1]

    # return none tuple 1. tagname is None, 2. no cpgs in intersections in the region
    ret_none_tuple = tuple([None] * (num_tuples_ret - 1) + [tagname])

    if tagname is None:
        return ret_none_tuple

    ontCalls_narrow_second_set = None  # if using joined sites of all tools, or None for not using joined sites
    if secondFilterBedFileName is not None:
        joined_set = load_single_sites_bed_as_set(secondFilterBedFileName)
        ontCalls_narrow_second_set = set(ontCalls.keys()).intersection(joined_set)

    # Initial evaluation vars
    TP_5mC = FP_5mC = FN_5mC = TN_5mC = TP_5C = FP_5C = FN_5C = TN_5C = 0
    y_of_bgtruth = []
    ypred_of_ont_tool = []
    yscore_of_ont_tool = []

    mCalls = 0  # count how many reads call is methy of a tool
    cCalls = 0  # count how many reads call is unmethy of a tool

    referenceCpGs = 0  # number of sites

    # really CpGs for fully methylation and unmethylation,
    mCsites_BGTruth = 0  # count sites 5mC
    Csites_BGTruth = 0  # count sites 5C

    # We find the narrowed CpG set to evaluate, try to reduce running time
    targetedSet = ontCallsKeySet

    ## ontCalls_narrow_set
    ##          - set: need to intersect
    ##          - None: no need to intersect
    if ontCalls_narrow_set is not None:
        targetedSet = targetedSet.intersection(ontCalls_narrow_set)
    if len(targetedSet) == 0:  # no cpg sites for evaluation
        return ret_none_tuple

    if ontCalls_narrow_second_set is not None:
        targetedSet = targetedSet.intersection(ontCalls_narrow_second_set)
    if len(targetedSet) == 0:  # no cpg sites for evaluation
        return ret_none_tuple

    for cpgKey in targetedSet:  # key = (chr, start, strand)
        ##### for each sites, we perform per read stats:
        if satisfy_fully_meth_or_unmeth(bgTruth[cpgKey][0]):
            referenceCpGs += 1

            if is_fully_meth(bgTruth[cpgKey][0]):
                mCsites_BGTruth += 1
            elif is_fully_unmeth(bgTruth[cpgKey][0]):
                Csites_BGTruth += 1
            else:
                raise Exception(f'We must see all certain sites here, but see meth_freq={bgTruth[cpgKey][0]}')

            for perCall in ontCalls[cpgKey]:  # perCall is a tupple of (pred_class, pred_score)
                if perCall[0] == 1:
                    mCalls += 1
                elif perCall[0] == 0:
                    cCalls += 1
                else:
                    raise Exception(f'Pred_class is only 0 or 1, but is {perCall}')

                ### variables needed to compute precission, recall etc.:
                if perCall[0] == 1 and is_fully_meth(bgTruth[cpgKey][0]):  # true positive
                    TP_5mC += 1
                elif perCall[0] == 1 and is_fully_unmeth(bgTruth[cpgKey][0]):  # false positive
                    FP_5mC += 1
                elif perCall[0] == 0 and is_fully_meth(bgTruth[cpgKey][0]):  # false negative
                    FN_5mC += 1
                elif perCall[0] == 0 and is_fully_unmeth(bgTruth[cpgKey][0]):  # true negative
                    TN_5mC += 1

                if perCall[0] == 0 and is_fully_unmeth(bgTruth[cpgKey][0]):  # true positive
                    TP_5C += 1
                elif perCall[0] == 0 and is_fully_meth(bgTruth[cpgKey][0]):  # false positive
                    FP_5C += 1
                elif perCall[0] == 1 and is_fully_unmeth(bgTruth[cpgKey][0]):  # false negative
                    FN_5C += 1
                elif perCall[0] == 1 and is_fully_meth(bgTruth[cpgKey][0]):  # true negative
                    TN_5C += 1

                ### prediction results, AUC related:
                ypred_of_ont_tool.append(perCall[0])
                if np.isnan(perCall[1]):
                    yscore_of_ont_tool.append(0.0)
                else:
                    yscore_of_ont_tool.append(perCall[1])

                if is_fully_meth(bgTruth[cpgKey][0]):  # BG Truth label
                    y_of_bgtruth.append(1)
                else:
                    y_of_bgtruth.append(0)
        else:
            raise Exception(f'We must see all certain sites here, but see meth_freq={bgTruth[cpgKey][0]}')

    ### compute all per read stats:
    with warnings.catch_warnings(record=True) as w:
        try:
            accuracy = (TP_5mC + TN_5mC) / float(TP_5mC + FP_5mC + FN_5mC + TN_5mC)
        except ZeroDivisionError:
            accuracy = 0

        try:
            predicted_condition_positive_5mC = float(TP_5mC + FP_5mC)
            precision_5mC = TP_5mC / predicted_condition_positive_5mC
        except ZeroDivisionError:
            precision_5mC = 0

        try:
            predicted_condition_positive_5C = float(TP_5C + FP_5C)
            precision_5C = TP_5C / predicted_condition_positive_5C
        except ZeroDivisionError:
            precision_5C = 0

        try:
            recall_5mC = TP_5mC / float(TP_5mC + FN_5mC)
        except ZeroDivisionError:
            recall_5mC = 0

        try:
            recall_5C = TP_5C / float(TP_5C + FN_5C)
        except ZeroDivisionError:
            recall_5C = 0

        # F1 score, precision and recall
        f1_micro = f1_score(y_of_bgtruth, ypred_of_ont_tool, average='micro')
        f1_macro = f1_score(y_of_bgtruth, ypred_of_ont_tool, average='macro')

        precision_micro = precision_score(y_of_bgtruth, ypred_of_ont_tool, average='micro')
        precision_macro = precision_score(y_of_bgtruth, ypred_of_ont_tool, average='macro')

        recall_micro = recall_score(y_of_bgtruth, ypred_of_ont_tool, average='micro')
        recall_macro = recall_score(y_of_bgtruth, ypred_of_ont_tool, average='macro')

        try:
            F1_5mC = 2 * ((precision_5mC * recall_5mC) / (precision_5mC + recall_5mC))
        except ZeroDivisionError:
            F1_5mC = 0

        try:
            F1_5C = 2 * ((precision_5C * recall_5C) / (precision_5C + recall_5C))
        except ZeroDivisionError:
            F1_5C = 0

        fprSwitch = 1
        try:
            fpr, tpr, _ = roc_curve(y_of_bgtruth, yscore_of_ont_tool)
            average_precision = average_precision_score(y_of_bgtruth, yscore_of_ont_tool)
        except ValueError:
            logger.error(
                f"###\tERROR for roc_curve: y(Truth):{len(y_of_bgtruth)}, scores(Call pred):{len(yscore_of_ont_tool)}, \nother settings: {title}, {tagname}, {secondFilterBedFileName}")
            fprSwitch = 0
            roc_auc = 0.0
            average_precision = 0.0

        if fprSwitch == 1:
            roc_auc = auc(fpr, tpr)

    ########################
    if save_curve_data:
        # save y and y-pred and y-score for later plot:
        curve_data = {'yTrue': y_of_bgtruth, 'yPred': ypred_of_ont_tool, 'yScore': yscore_of_ont_tool}

        os.makedirs(os.path.join(outdir, 'curve_data'), exist_ok=True)
        outfn = os.path.join(outdir, 'curve_data', f'{prefix_name}.{tagname.replace(" ", "_")}.curve_data.pkl')
        with open(outfn, 'wb') as handle:
            pickle.dump(curve_data, handle)

    return (accuracy, roc_auc, average_precision, f1_macro, f1_micro, \
            precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, \
            recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, \
            F1_5mC, mCalls, referenceCpGs, Csites_BGTruth, mCsites_BGTruth, tagname,)


def save_keys_to_single_site_bed(keys, outfn, callBaseFormat=1, outBaseFormat=1, nonstr='.'):
    """
    Save all keys in set of ('chr  123  123  .  .  +\n', etc.) to outfn.
    We use non-string like . in 3rd, 4th columns by BED file format.
    :param keys:
    :param outfn:
    :return:
    """
    if outfn.endswith('.gz'):
        outfile = gzip.open(outfn, 'wt')
    else:
        outfile = open(outfn, 'w')
    for key in keys:
        if outBaseFormat == 0:
            outfile.write(
                f'{key[0]}\t{key[1] - callBaseFormat + outBaseFormat}\t{key[1] - callBaseFormat + outBaseFormat + 1}\t{nonstr}\t{nonstr}\t{key[2]}\n')
        else:
            outfile.write(
                f'{key[0]}\t{key[1] - callBaseFormat + outBaseFormat}\t{key[1] - callBaseFormat + outBaseFormat}\t{nonstr}\t{nonstr}\t{key[2]}\n')
    outfile.close()


def do_singleton_nonsingleton_scanner():
    """
    Do generate singleton and nonsingleton BED file
    :return:
    """
    kbp = 5
    singletonFilename = os.path.join(pic_base_dir, f'hg38_singletons_{kbp}bp.bed.gz')
    nonsingletonFilename = os.path.join(pic_base_dir, f'hg38_nonsingletons_{kbp}bp.bed.gz')
    SingletonsAndNonSingletonsScanner(reference_genome_hg38_fn, singletonFilename, nonsingletonFilename, kbp=kbp)

    kbp = 10
    singletonFilename = os.path.join(pic_base_dir, f'hg38_singletons_{kbp}bp.bed.gz')
    nonsingletonFilename = os.path.join(pic_base_dir, f'hg38_nonsingletons_{kbp}bp.bed.gz')
    SingletonsAndNonSingletonsScanner(reference_genome_hg38_fn, singletonFilename, nonsingletonFilename, kbp=kbp)


def SingletonsAndNonSingletonsScanner(referenceGenomeFile, outfileName_s, outfileName_ns, kbp=10):
    """
    Generate singleton and non-singletons BED file, based on Reference Genome and KBP up and down streams.

    The output file is in 1-based at start coordinate system.
    kbp is up and down k-bp regions and evaluated on positive strand.
    Singletons: only one CpG in the region
    Nonsingletons: more than one CpG in the region
    """
    reference = SeqIO.to_dict(SeqIO.parse(referenceGenomeFile, "fasta"))
    logger.debug(
        f"###\tSingletonsAndNonSingletonsScanner: {referenceGenomeFile} reference genome file is parsed, up and down bp={kbp}")

    outfile_s = gzip.open(outfileName_s, "wt")  # "s" stands for Singletons
    outfile_ns = gzip.open(outfileName_ns, "wt")  # "ns" stands for Non-Singletons

    for chromosome in list(reference.keys()):
        if chromosome not in HUMAN_CHR_SET:
            continue
        idxs = re.finditer('CG', str(reference[chromosome].seq).upper())

        singleton = -1  # 1 will stand for yes, 0 for no
        for idx in idxs:
            #             print(chromosome, idx, idx.start(), idx.end())
            if singleton == -1:
                s = idx.start() + 1  # here 8: mock1 <_sre.SRE_Match object; span=(8, 10), match='CG'> 8 10
                end_index = idx.end()  # here 10: mock1 <_sre.SRE_Match object; span=(8, 10), match='CG'> 8 10
                singleton = 1
            else:
                if (idx.start() - end_index) < kbp:
                    # we just found a non-singleton. I.e. accordingly to the Nanopolish approach, CGs closer than 5bp, are considered as non-singletons
                    # Singletons are SR=XXXXXCGXXXXX
                    # Non-singletons are SR=XXXXXCGXXXXCGXXXCGXXCGCGXXXXX  , <5bp for pair of neighbor CGs
                    end_index = idx.end()
                    singleton = 0
                else:
                    # current CG is not part of non-singleton. It might mean that its not part of a big non-singleton or singleton upstream from it. We test which of these options below
                    # The current CG is up to k-bp to previous CG, we store the previous regions into Singletons or Non-singletons
                    if singleton == 1:
                        #                         print(chromosome, s, e, "SINGLETON")
                        outfile_s.write("{}\t{}\t{}\n".format(chromosome, s, end_index))
                    else:
                        #                         print(chromosome, s, e, "NON-SINGLETON")
                        outfile_ns.write("{}\t{}\t{}\n".format(chromosome, s, end_index))
                    s = idx.start() + 1
                    end_index = idx.end()
                    singleton = 1

        if singleton == 1:  # this code repetition takes care of the last instance in the long list of CG indexes
            #             print(chromosome, s, e, "SINGLETON")
            outfile_s.write("{}\t{}\t{}\n".format(chromosome, s, end_index))
        else:
            #             print(chromosome, s, e, "NON-SINGLETON")
            outfile_ns.write("{}\t{}\t{}\n".format(chromosome, s, end_index))

        logger.debug(f"###\tNonSingletonsScanner: chromosome {chromosome} processed")

    outfile_s.close()
    outfile_ns.close()

    logger.debug(
        f"###\tSingletonsAndNonSingletonsScanner: {referenceGenomeFile} file processed, kbp={kbp}, save to Singletons:{outfile_s}, and Nonsingletons:{outfile_ns}")


def eval_concordant_within_kbp_region(cpgList, evalcpg, kbp=10):
    """
    Evaluate the kth elem if is concordant, by checking state of sites in k-bp region meth states

    For example with in 5bp: CGXXXXCG 4 times X is considered in 5bp region
    :param cpgList:
    :param k:
    :return: True if concordant, False is discordant
    """
    # cpgList=list(cpgList)

    for cpg in cpgList:  # cpg is (start,strand, meth_indicator)
        if abs(evalcpg[0] - cpg[
            0]) - 2 >= kbp:  # not within kbp regions, skip (Note: kbp is number of X between CG, ex. CGXXCGXXCG)
            continue
        if evalcpg[1] != cpg[1]:  # found there is a CPG not same state, return Discordant
            return False
    return True  # Concordant


def nonSingletonsPostprocessing(absoluteBGTruth, nsRegionsBedFileName, nsConcordantFileName, nsDisCordantFileName,
                                kbp=10, print_first=False,
                                genome_annotation_dir=os.path.join(data_base_dir, 'genome-annotation')):
    """
    Define concordant and discordant based on BG-Truth.
    Return 1-based Cocordant and Discordant regions in bed file

    Based on only 100% or 0% bg-truth in BS-seq (absoluteBGTruth), we define:
    Concordant: All CpGs in 10-bp up and down region is same states, such as 0000, or 1111.
    Discordant: CpGs in the 10-bp region is mixed with fully-methylated (1) and unmethylated(0) state.

    This kind of preprocessing will have to be done for each studied library separately.

    Output format for Concordant and Discordant file:
    chr   start  end  meth_state    coverage
    chr1  123   124   1             16
    """
    logger.debug(f"nonSingletonsPostprocessing, based on file={nsRegionsBedFileName}, kbp={kbp}")
    bedBGTruth = BedTool(calldict2txt(absoluteBGTruth), from_string=True).sort()

    infn = os.path.join(genome_annotation_dir, nsRegionsBedFileName)
    regionNonsingletons = BedTool(infn).sort()

    regionWithBGTruth = regionNonsingletons.intersect(bedBGTruth, wa=True,
                                                      wb=True)  # chr start end   chr start end .  .  strand

    regionDict = defaultdict(
        list)  # key->value, key=region of (chr, start, end), value=list of [f1,f2,etc.] , suche as {regionCoords : [methylation percentage list]}

    is_print_first = print_first
    cntBedLines = 0

    for ovr in regionWithBGTruth:
        cntBedLines += 1
        if is_print_first:
            logger.debug(f'ovr={ovr}')

        regionKey = (ovr[0], int(ovr[1]), int(ovr[2]))  # chr  start  end
        methKey = (ovr[3], int(ovr[4]), ovr[8])  # chr, start, strand

        tmpStart = int(ovr[4])
        tmpStrand = ovr[8]

        if tmpStrand == '-':  # group + - strand CpG together, for simplicity
            tmpStart -= 1
        meth_freq = absoluteBGTruth[methKey][0]
        coverage = absoluteBGTruth[methKey][1]
        if is_fully_meth(meth_freq):
            methIndicator = 1
        else:
            methIndicator = 0

        regionDict[regionKey].append((tmpStart, methIndicator, coverage))  # save value as (start, 0/1, cov)

        if is_print_first:
            logger.debug(f'regionDict[regionKey]={regionDict[regionKey]}, regionKey={regionKey}')
        is_print_first = False

    ## logger.debug(f'cntBedLines={cntBedLines:,}')

    is_print_first = print_first
    concordantList = {}  # list of (chr, start)
    discordantList = {}
    meth_cnt_dict = defaultdict(int)  # count meth states in concordant and discordant
    unmeth_cnt_dict = defaultdict(int)

    bar = tqdm(regionDict)
    bar.set_description("Scan non-singletons")
    for region in bar:  # region is (chr, start, end)
        cpgList = regionDict[region]  # get values as the list of  (start, 0/1, cov)
        chrOut = region[0]

        if is_print_first:
            logger.debug(f'cpgList={cpgList}')

        for k in range(len(cpgList)):  # for each CpG
            # Since before we group + - to +, now we assure use start
            startOut = cpgList[k][0]
            methIndicator = cpgList[k][1]
            methCov = cpgList[k][2]

            if eval_concordant_within_kbp_region(cpgList, cpgList[k], kbp=kbp):
                # add this cpg to concordant
                concordantList.update(
                    {(chrOut, startOut): (chrOut, startOut, methIndicator, methCov)})  # list of (chr, start)
                if cpgList[k][1] == 1:
                    meth_cnt_dict['Concordant'] += 1
                else:
                    unmeth_cnt_dict['Concordant'] += 1
            else:
                # add this cpg to discordant
                discordantList.update({(chrOut, startOut): (chrOut, startOut, methIndicator, methCov)})
                if cpgList[k][1] == 1:
                    meth_cnt_dict['Discordant'] += 1
                else:
                    unmeth_cnt_dict['Discordant'] += 1
            pass
        if is_print_first:
            logger.debug(f'concordantList={concordantList}, discordantList={discordantList}')
        is_print_first = False

    concordantSet = list(set(concordantList) - set(discordantList))  # remove duplicate and same in discordant
    # outfile_concordant = open(nsConcordantFileName, "w")
    outfile_concordant = gzip.open(nsConcordantFileName + ".tmp.gz", "wt")
    for cpgKey in concordantSet:  # (chrOut, startOut) -> (chrOut, startOut, methIndicator, methCov)
        cpg = concordantList[cpgKey]
        region_txt = '\t'.join([cpg[0], str(cpg[1]), str(cpg[1] + 1), str(cpg[2]), str(cpg[3])]) + '\n'
        outfile_concordant.write(region_txt)
    outfile_concordant.close()

    discordantSet = list(set(discordantList))  # remove duplicate
    # outfile_discordant = open(nsDisCordantFileName, "w")
    outfile_discordant = gzip.open(nsDisCordantFileName + ".tmp.gz", "wt")
    for cpgKey in discordantSet:
        cpg = discordantList[cpgKey]
        region_txt = '\t'.join([cpg[0], str(cpg[1]), str(cpg[1] + 1), str(cpg[2]), str(cpg[3])]) + '\n'
        outfile_discordant.write(region_txt)
    outfile_discordant.close()

    ## sort bed file
    sort_bed_file(nsConcordantFileName + ".tmp.gz", nsConcordantFileName, deduplicate=True)
    os.remove(nsConcordantFileName + ".tmp.gz")
    sort_bed_file(nsDisCordantFileName + ".tmp.gz", nsDisCordantFileName, deduplicate=True)
    os.remove(nsDisCordantFileName + ".tmp.gz")

    logger.debug(f'save to {[nsConcordantFileName, nsDisCordantFileName]}')
    # logger.debug(f'meth_cnt={meth_cnt_dict}, unmeth_cnt={unmeth_cnt_dict}')
    ret = {'Concordant.5mC': meth_cnt_dict['Concordant'], 'Concordant.5C': unmeth_cnt_dict['Concordant'],
           'Discordant.5mC': meth_cnt_dict['Discordant'], 'Discordant.5C': unmeth_cnt_dict['Discordant']}
    return ret


## Will deprecated
def load_nanopolish_df(
        infn='/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/APL/APL.nanopolish_methylation_calls.tsv'):
    """
    Load the nanopolish original output results tsv into a dataframe

    head /projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/APL/APL.nanopolish_methylation_calls.tsv
    chromosome	start	end	read_name	log_lik_ratio	log_lik_methylated	log_lik_unmethylated	num_calling_strands	num_cpgs	sequence
    chr1	14348	14353	542a58d6-b2f8-4ccf-829f-7a1977876c6c	8.14	-244.99-253.13	1	2	GACCCCGAGACGTTTG
    chr1	14434	14434	542a58d6-b2f8-4ccf-829f-7a1977876c6c	0.79	-126.97-127.77	1	1	TGTGCCGTTTT
    chr1	14468	14468	542a58d6-b2f8-4ccf-829f-7a1977876c6c	-0.97	-186.33-185.36	1	1	AGTGGCGCAGG
    :param infn:
    :return:
    """
    df = pd.read_csv(infn, sep='\t')
    # logger.debug(df)
    # logger.info(df.iloc[:, -3].value_counts())

    return df


## Will deprecated
def load_tombo_df(
        infn='/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/K562/K562.tombo_perReadsStats.bed'):
    """
    Load the nanopolish original output results tsv into a dataframe

    head /projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/APL/APL.nanopolish_methylation_calls.tsv
    chromosome	start	end	read_name	log_lik_ratio	log_lik_methylated	log_lik_unmethylated	num_calling_strands	num_cpgs	sequence
    chr1	14348	14353	542a58d6-b2f8-4ccf-829f-7a1977876c6c	8.14	-244.99-253.13	1	2	GACCCCGAGACGTTTG
    chr1	14434	14434	542a58d6-b2f8-4ccf-829f-7a1977876c6c	0.79	-126.97-127.77	1	1	TGTGCCGTTTT
    chr1	14468	14468	542a58d6-b2f8-4ccf-829f-7a1977876c6c	-0.97	-186.33-185.36	1	1	AGTGGCGCAGG
    :param infn:
    :return:
    """
    df = pd.read_csv(infn, sep='\t', header=None)
    # logger.debug(df)
    # logger.info(df.iloc[:, -3].value_counts())

    return df


## Will deprecated
def load_deepmod_df(infn):
    """
    Load the DeepMod original output results tsv into a dataframe
    :param infn:
    :return:
    """
    df = pd.read_csv(infn, sep=' ', header=None)
    return df


def load_sam_as_strand_info_df(infn='/projects/li-lab/yang/workspace/nano-compare/data/bam-files/K562.sam'):
    """
    Load strand info from SAM files, and make a df of read-name and strand-info as follows:

    return format is followings:
    2020-12-11 12:12:54,582 - [meth_stats_common.py:3066] - INFO:                                    read-name strand-info
    0       dd31ef2a-8826-4f1a-afd9-1d7a9bdea414           +
    1       628585d4-1947-463d-9787-04040600cb65           +
    2       db24bf7d-086f-424f-afc4-e48dca59de80           +
    3       a3ea655a-32ac-4b6b-b9ec-92dae00fc91c           +
    4       be5ad23e-2379-4a6a-87e9-024676bdbbdb           +
    ...                                      ...         ...
    263013  60e6fb8d-0e99-4b02-91ec-562325503640           -
    263014  c22e60d4-285c-4947-a949-a886793d4a1e           -
    263015  635e260a-82c4-4473-894d-952008eff026           -
    263016  39cfbc37-dcc7-4b3f-8e27-af1bd82e9d95           +
    263017  99351046-7f42-4187-b880-3ebfcf33de86           +

    [263018 rows x 2 columns]
    :param infn:
    :return:
    """
    data = {}  # key=read-name, value= + or -
    strand_info = '+'
    samfile = pysam.AlignmentFile(infn, "r")
    for read in samfile.fetch():
        # logger.debug(read)
        # logger.debug(read.flag)

        if (read.flag & 0x10) == 0x10:  # according to spec: https://samtools.github.io/hts-specs/SAMv1.pdf
            strand_info = '-'
        else:
            strand_info = '+'
        # logger.debug(type(read))
        # logger.debug(dir(read))
        # logger.debug(read.query_name)
        data.update({read.query_name: strand_info})
        # break
    samfile.close()

    df = pd.Series(data).to_frame()
    df = df.reset_index()
    df.columns = ['read-name', 'strand-info']
    # df.columns[0] = 'read-name'
    # df.columns[1] = 'strand-info'
    logger.debug(len(data))
    logger.debug(df)
    logger.debug(df['strand-info'].value_counts())
    return df


def get_dna_base_from_reference(chr, start, num_seq=5, ref_fasta=None):
    """
    Get the base of DNA from start-num_seq to start+num_seq, totally 2*num_seq+1
    The center start is 0-based start position

    :param chr:
    :param start:
    :param end:
    :return:
    """

    if ref_fasta is None:
        raise Exception('Please specify params: ref_fasta')

    long_seq = ref_fasta[chr].seq
    short_seq = str(long_seq)[start - num_seq:start + num_seq + 1]

    # logger.info(short_seq)
    return short_seq


def get_dna_seq_from_reference(chr, start, end, ref_fasta=None):
    """
    Get the sequence from [start, end), totally 2*num_seq+1
    The center start is 0-based start position

    :param chr:
    :param start:
    :param end:
    :return:
    """

    if ref_fasta is None:
        raise Exception('Please specify params: ref_fasta')

    long_seq = ref_fasta[chr].seq
    short_seq = str(long_seq)[start:end]

    # logger.info(short_seq)
    return short_seq


def get_ref_fasta(ref_fn=reference_genome_hg38_fn):
    ref_fasta = SeqIO.to_dict(SeqIO.parse(open(ref_fn), 'fasta'))
    logger.debug(f'load ref file from {ref_fn}')
    return ref_fasta


refGenome = None


def sanity_check_dna_sequence(chr='chr10', start_base0=10493):
    """
    start is a 0-based position
    :param chr:
    :param start:
    :return:
    """
    ret_seq = get_dna_base_from_reference(chr, start_base0, ref_fasta=refGenome)
    logger.debug(f'Report is 0-based for the input:\n{chr}:{start_base0}\n{ret_seq}\n-----^-----\n\n')


def get_cache_filename(infn, params):
    """
    Get the file name that encoded in cache based on parameters, the params must have cache_dir for searching cache files.
    :param infn:
    :param params:
    :return:
    """
    if 'cache_dir' not in params:
        raise Exception(f"cache_dir is not in params={params}")
    cache_dir = params['cache_dir']
    if cache_dir is None:
        raise Exception(f"'cache_dir={cache_dir} is not properly set")

    basefn = os.path.basename(infn)
    if params['file_type'] in ['ont-call', 'bs-seq']:  # for ont calls and bs-seq
        cachefn = f'cachefile.{params["file_type"]}.encode.{params["encode"]}.{basefn}.base.{params["baseFormat"]}'
        if params["encode"] in ToolEncodeList:
            cachefn += f'.inscore.{params["include_score"]}'
            if params["encode"] in ['DeepMod.Cluster', 'DeepMod.C', 'DeepMod', 'Guppy', 'Guppy.ZW', 'UNISITE']:
                cachefn += f'.siteLevel.{params["siteLevel"]}'
        elif params["encode"] in BGTruthEncodeList:
            cachefn += f'.cov.{params["cov"]}.incov.{params["includeCov"]}'
        else:
            raise Exception(f'Encode {params["encode"]} is not support now')
    elif params['file_type'] in ['genome-annotation']:
        cachefn = f'cachefile.genome.annotation.{basefn}.base.{params["baseFormat"]}'
    ret_cachefn = os.path.join(cache_dir, cachefn + '.pkl')
    return ret_cachefn


def save_to_cache(infn, data, **params):
    """
    Save the data from program into cache, encoded by parameters
    :param infn:
    :param data:
    :param params:
    :return:
    """
    if not data:
        return
    if 'cache_dir' not in params:
        Exception(f"cache_dir is not in params={params}")
    cache_dir = params['cache_dir']

    if cache_dir is None:
        raise Exception(f"cache_dir is not properly set")

    os.makedirs(cache_dir, exist_ok=True)

    cache_fn = get_cache_filename(infn, params)
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir, exist_ok=True)
    with open(cache_fn, 'wb') as outf:
        pickle.dump(data, outf)
    logger.debug(f'Cached to file:{cache_fn}')


def check_cache_available(infn, **params):
    """
    Check if the input file is in cache, and return if it is exists
    :param infn:
    :param params:
    :return:
    """
    cachefn = get_cache_filename(infn, params)

    if os.path.exists(cachefn):
        # logger.debug(f'Start to get from cache:{cachefn}')
        try:
            with open(cachefn, 'rb') as inf:
                ret = pickle.load(inf)
                logger.debug(f'Get from cache:{cachefn}')
            return ret
        except:
            return None
    return None


def filter_cpg_dict(cpgDict, filterDict, toolname="Tool"):
    """
    Filter and keep only filterDict keys
    :param cpgDict:
    :param filterDict:
    :return:
    """
    retDict = defaultdict()
    joinedKeys = set(filterDict.keys()).intersection(set(cpgDict.keys()))
    for k in tqdm(joinedKeys, desc=f"FilterCPG-{toolname}"):
        retDict[k] = cpgDict[k]
    return retDict


def filter_cpg_dict_by_cov(cpgDict, coverage=1):
    """
    Filter cpg with coverage
    Args:
        cpgDict:
        coverage:

    Returns:

    """
    if coverage <= 1:
        return cpgDict
    retDict = {}
    for key in cpgDict:
        if len(cpgDict[key]) >= coverage:
            retDict[key] = cpgDict[key]
    return retDict


def is_fully_meth(methfreq, eps=1e-5, cutoff_fully_meth=1.0):
    """
    Check if the freq is fully-methylated, can be 1.0, or 0.9
    :param methfreq:
    :param eps:
    :param cutoff_fully_meth:
    :return:
    """
    if methfreq > 1.0 + eps or methfreq < 0:
        raise Exception(f'detect error value for freq={methfreq}')

    if methfreq > cutoff_fully_meth - eps:  # near 1
        return True
    return False


def is_fully_unmeth(methfreq, eps=1e-5):
    """
    Check if the freq is unmethylated, means 0.0 (almost)
    :param methfreq:
    :param eps:
    :return:
    """
    if methfreq > 1.0 + eps or methfreq < 0:
        raise Exception(f'detect error value for freq={methfreq}')

    if methfreq < eps:  # near 0
        return True
    return False


def satisfy_fully_meth_or_unmeth(methfreq, eps=1e-5, cutoff_fully_meth=1.0):
    """
    Return true if fully meth or unmeth, also known as certain sites
    :param methfreq:
    :return:
    """
    if is_fully_meth(methfreq, eps=eps, cutoff_fully_meth=cutoff_fully_meth) or is_fully_unmeth(methfreq, eps=eps):
        return True
    return False


def combineBGTruthList(bgTruthList, covCutoff=1):
    """
    Combine two replicates together, we unioned two replicates together as one bgtruth, and retain only cov >= covCutoff sites.
    Combined coverage = cov1 + cov2.
    :param bgTruthList:
    :return:
    """
    unionBGTruth = {}  # used for singleton and non-singletons detect, used for performance eval
    jointBGTruth = {}  # intersection of CpG sites
    if len(bgTruthList) == 2:  # sites must in both replicates, and
        logger.debug(f'Start study union of multiple BG-Truth, with coverage={covCutoff}')
        unionSet = set(bgTruthList[0].keys()).union(set(bgTruthList[1].keys()))
        jointSet = set(bgTruthList[0].keys()).intersection(set(bgTruthList[1].keys()))

        for key in tqdm(unionSet, desc="Combine-BSseq"):
            meth1 = meth2 = 0
            cov1 = cov2 = 0
            if key in bgTruthList[0]:
                cov1 = bgTruthList[0][key][1]
                meth1 = round(bgTruthList[0][key][0] * cov1)

            if key in bgTruthList[1]:
                cov2 = bgTruthList[1][key][1]
                meth2 = round(bgTruthList[1][key][0] * cov2)

            # Joined of two replicates
            meth_freq = (meth1 + meth2) / float(cov1 + cov2)

            if meth_freq > 1.0:
                raise Exception(
                    f"Compute joined meth_freq={meth_freq} error, using meth1={meth1}, cov1={cov1}, meth2={meth2}, cov2={cov2}, in sites key={key}, and bgTruthList[0][key]={bgTruthList[0][key]}, bgTruthList[1][key]={bgTruthList[1][key]}")

            cov = int(cov1 + cov2)

            if cov < covCutoff:
                continue

            unionBGTruth[key] = (meth_freq, cov)
            if key in jointSet:
                jointBGTruth[key] = (meth_freq, cov)
        logger.debug(
            f'unionBGTruth = {len(unionBGTruth):,}, jointBGTruth={len(jointBGTruth):,}, with cov-cutoff={covCutoff}')
    elif len(bgTruthList) == 1:
        for key in tqdm(bgTruthList[0], desc="Combine-BSseq"):
            if bgTruthList[0][key][1] >= covCutoff:
                unionBGTruth[key] = bgTruthList[0][key]
        logger.debug(f'Only 1 replicates, BGTruth = {len(unionBGTruth):,}, with cov-cutoff={covCutoff}')
    else:
        raise Exception(f'len={len(bgTruthList)}, is not support now.')

    return unionBGTruth


def combineBGTruthList_by_DeepModPaper(bgTruthList, freqCutoff=0.9, covCutoff=1, filterChrs=HUMAN_CHR_SET):
    """
    Combine two replicates by DeepMod, >90% in both as methylated, =0% in both as unmethylated, remove others
    :param bgTruthList:
    :return:
    """
    jointBGTruth = {}  # intersection of CpG sites
    if len(bgTruthList) == 2:  # sites must in both replicates, and
        unionSet = set(bgTruthList[0].keys()).union(set(bgTruthList[1].keys()))
        jointSet = set(bgTruthList[0].keys()).intersection(set(bgTruthList[1].keys()))
        logger.debug(
            f'unionSet={len(unionSet):,}, jointSet={len(jointSet):,}, cov_cutoff={covCutoff}, for chr={filterChrs}')
        methCnt = unmethCnt = 0

        totalCnt = 0

        for key in jointSet:
            if key[0] not in filterChrs:
                continue
            totalCnt += 1
            cov1 = bgTruthList[0][key][1]
            methfreq1 = bgTruthList[0][key][0]

            cov2 = bgTruthList[1][key][1]
            methfreq2 = bgTruthList[1][key][0]

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
        raise Exception(f'len={len(bgTruthList)}, is not support now.')

    ret = (methCnt, unmethCnt, totalCnt - methCnt - unmethCnt)
    return jointBGTruth, ret


def filter_cpgkeys_using_bedfile(cpgKeys, bedFileName):
    """
    Keep only cpg keys in bed file range, return set of keys
    :param cpgKeys:
    :param bedFileName:
    :return:
    """
    cpgBed = calldict2bed(cpgKeys)
    coordBed = get_region_bed(bedFileName)
    intersectBed = intersect_bed_regions(cpgBed, coordBed, bedFileName)

    ret = set(bedtxt2dict(intersectBed).keys())
    return ret


def find_bed_filename(basedir, pattern):
    """
    Find BED file of concordant and discordant for each dataset or run.
    :param basedir:
    :param pattern:
    :return:
    """
    fnlist = glob.glob(os.path.join(basedir, '**', pattern), recursive=True)
    if len(fnlist) < 1:
        logger.debug(
            f'ERROR: Find no files: {fnlist}, please check the basedir={basedir} if it is correct, and the search pattern={pattern}')
        return None
    logger.debug(f'find bed file:{fnlist[0]}, len={len(fnlist)}')
    return fnlist[0]


def compute_and_gen_venn_data(set_dict, namelist, outdir, tagname='tagname'):
    """
    Compute and generate 7 data for three set or 31 data for five set joining Venn Diagram plotting
    Return 2^n-1 of set intersections
    :param set_dict:
    :param outdir:
    :param tagname:
    :return:
    """
    retlist = []
    for k in range(len(set_dict)):
        for combin in combinations(namelist, k + 1):
            join_set = set(set_dict[combin[0]])
            for t in range(1, len(combin)):
                join_set = join_set.intersection(set_dict[combin[t]])
            retlist.append(len(join_set))
    if (len(namelist) == 5 and len(retlist) != 31) or (len(namelist) == 3 and len(retlist) != 7):
        raise Exception(
            f'Number of set need to have combinations corresponding, but get cnt={len(retlist)} for set number={len(namelist)}, code bugs.')

    outfn = os.path.join(outdir, f'venn.data.{tagname}.dat')
    with open(outfn, 'w') as outf:
        for num in retlist:
            outf.write(f'{num}\n')
    logger.debug(f'Note the venn data set, tool-name order must be: {namelist}')
    logger.debug(f'save {len(retlist)} points venn data for {tagname} to {outfn}')


def ontcalls_to_setsfile_for_venn_analysis(call_keys, outfn):
    """
    Ouput keys into a sets file for venn analysis later.
    key is (chr, pos, strand), output as chr1_1234_+
    :param call_keys:
    :param outfn:
    :return:
    """
    outf = gzip.open(outfn, 'wt')
    for key in call_keys:
        outf.write(f"{key[0]}_{key[1]}_{key[2]}\n")
    outf.close()
    logger.debug(f"save to {outfn}")


def bedtool_convert_0_to_1(bed):
    """
    Convert 0-based start BED file into 1-based start, assume end is always 1-based
    :param bed: BedTool object
    :return:
    """
    df = bed.to_dataframe()
    df['start'] = df['start'] + 1
    ret = BedTool.from_dataframe(df)
    return ret


def get_region_tagname(infn):
    """
    Get the infn's region name
    Args:
        infn: None for genome-wide

    Returns: None for not recognized in config file

    """
    if infn is None or infn == 'x.x.Genome-wide':
        return genome_wide_tagname

    basefn = os.path.basename(infn)
    if basefn.endswith('.concordant.bed.gz'):
        return 'Concordant'
    elif basefn.endswith('.discordant.bed.gz'):
        return 'Discordant'

    if basefn in region_filename_dict:
        return region_filename_dict[basefn][0]

    logger.debug(f"ERROR: Not correct infn={infn}")
    return None


def is_0base_region_files(infn):
    """
    Return True if region file is 0-based start
    :param infn:
    :return:
    """
    if os.path.basename(infn) in region_filename_dict:  # get from config file
        return region_filename_dict[os.path.basename(infn)][1] == 0
    # default is 1-based BED region
    return False


def get_region_bed(infn, enable_base_detection_bedfile=enable_base_detection_bedfile):
    """
    Get the BED file for infn, return 1-based start BED object
    :param infn:
    :return:
    """
    if not os.path.isfile(infn):
        # Can not find this file
        logger.debug(f"ERROR: Can not find region file={infn}")
        return None
    coordBed = BedTool(infn).sort()
    if enable_base_detection_bedfile and is_0base_region_files(infn):
        coordBed = bedtool_convert_0_to_1(coordBed)
    return coordBed


def get_region_bed_tuple(infn, enable_base_detection_bedfile=enable_base_detection_bedfile,
                         enable_cache=False, using_cache=False, cache_dir=None):
    """
    return (infn, tagname, bedobject) of regions file
    :param infn:
    :return:
    """
    if infn is None:  # means genome wide
        return (None, genome_wide_tagname, None)

    if enable_base_detection_bedfile:
        baseFormat = 1
    else:
        baseFormat = -1
    if enable_cache and using_cache:
        ret = check_cache_available(infn=infn, file_type='genome-annotation',
                                    baseFormat=baseFormat, cache_dir=cache_dir)
        ## If pkl in cache is avalable, and the BED file's tempdir/filename.tmp is available, import it
        ## or else, need to reload. Note: temp file may be released, even the cache pkl is there.
        if ret is not None and ret[2] is not None and os.path.exists(ret[2].fn):
            # logger.debug(f'BED import {os.path.basename(infn)} finished from cache file, BED tuple={ret}')
            return ret
        logger.debug(f'BED file not cached yet, we load from raw file={os.path.basename(infn)}')

    tagname = get_region_tagname(infn)
    region_bed = get_region_bed(infn, enable_base_detection_bedfile)

    ret = (infn, tagname, region_bed)
    if region_bed is None:
        logger.debug(
            f"WARN: region {tagname} file is not loaded, check if exists for filename={infn}")
        return ret
    if enable_cache:
        save_to_cache(infn, ret, file_type='genome-annotation', baseFormat=baseFormat, cache_dir=cache_dir)
    return ret


def update_progress_bar(*a):
    """
    Update progress for multiprocessing
    :param a:
    :return:
    """
    global progress_bar_global
    progress_bar_global.update()


def get_region_bed_pairs_list_mp(infn_list, processors=1, enable_base_detection_bedfile=enable_base_detection_bedfile,
                                 enable_cache=False, using_cache=False, cache_dir=None):
    """
    Get list of pairs [(infn, tagname, bedobject), ...] of regions file
    :param infn_list:
    :return:
    """
    logger.debug(f"get_region_bed_pairs_list_mp start, processors={processors}")
    logger.debug(f"regions={infn_list}, len={len(infn_list)}")
    ret_list = []
    with Pool(processors) as pool:
        # region_bed_list = pool.map(get_region_bed_pairs, infn_list)
        global progress_bar_global
        progress_bar_global = tqdm(total=len(infn_list))
        progress_bar_global.set_description("Load all BED regions")
        for infn in infn_list:
            ret_list.append(pool.apply_async(get_region_bed_tuple, args=(
                infn, enable_base_detection_bedfile, enable_cache, using_cache, cache_dir,),
                                             callback=update_progress_bar))
        pool.close()
        pool.join()
        progress_bar_global.close()
    region_bed_list = [ret.get() for ret in ret_list]
    logger.debug("get_region_bed_pairs_list_mp finished")
    return region_bed_list


# Deprecated
def get_region_bed_pairs_list_mt(infn_list, max_threads=1):
    """
    Get list of pairs [(basefn, tagname, bedobject), ...] of regions file
    :param infn_list:
    :return:
    """
    logger.debug(f"get_region_bed_pairs_list_mt start, max_threads={max_threads}")
    with ThreadPoolExecutor(max_workers=max_threads) as t:
        taskList = []
        for infn in infn_list:
            task = t.submit(get_region_bed_tuple, infn)
            taskList.append(task)
        wait(taskList)
        task_ret_list = [task.result() for task in taskList]
    logger.debug("get_region_bed_pairs_list_mt finished")
    return task_ret_list


def intersect_bed_regions(bed_a, bed_region, bedfn=""):
    """
    Intersect A with region, consider if bedfn need strandness.
    :param bed_a:
    :param bed_region:
    :param bedfn: default is strand insensitive
    :return:
    """
    ## Check if need strand for intersections, repetitive regions
    strand_sensitive_bed = False  # default is not strand sensitive
    if os.path.basename(bedfn) in region_filename_dict:
        strand_sensitive_bed = region_filename_dict[os.path.basename(bedfn)][2]

    if strand_sensitive_bed:
        intersectBed = bed_a.intersect(bed_region, u=True, wa=True, s=True)
    else:
        intersectBed = bed_a.intersect(bed_region, u=True, wa=True)
    return intersectBed


def filter_corrdata_df_by_bedfile(df, coord_bed, coord_fn):
    """
    Filter lines in correlation data, within coordinate BED file
    :param df:
    :param coord_fn:
    :return:
    """
    if coord_bed is None:  # No need to filter for genome wide
        return df
    # suppress warnings /home/liuya/anaconda3/envs/nanocompare/lib/python3.6/site-packages/pybedtools/bedtool.py:3287: UserWarning: Default names for filetype bed are:
    # ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
    # but file has 20 fields; you can supply custom names with the `names` kwarg
    #   % (self.file_type, _names, self.field_count()))
    warnings.filterwarnings('ignore', category=UserWarning)

    ## In order to support NaN in dataframe
    ## Need replace NaN to ., then for sort, or else will encounter error: Differing number of BED fields encountered at line: 1339
    ## Ref: NA values https://github.com/daler/pybedtools/issues/257
    bed_of_df = BedTool.from_dataframe(df).sort()
    bed_of_intersect = intersect_bed_regions(bed_of_df, coord_bed, coord_fn)

    if len(bed_of_intersect) > 0:
        ## replace NAN will be used to deal with non-joined corr meth data csv
        retdf = bed_of_intersect.to_dataframe(names=df.columns).replace('.', np.NaN)
    else:
        retdf = None
    return retdf


def get_meteore_format_set(infn, read_level=True):
    df = pd.read_csv(infn, sep='\t')
    ret = set()
    for index, row in df.iterrows():
        id = row['ID']
        chr = row['Chr']
        pos = int(row['Pos'])
        strand = row['Strand']
        if read_level:
            key = (id, chr, pos, strand)
        else:
            key = (chr, pos, strand)

        if key not in ret:
            ret.add(key)
    return ret


def sanity_check_meteore_combine_sites():
    bdir = "/projects/li-lab/Nanopore_compare/suppdata/METEORE_results/METEORE_raw_input"

    apl_deepsignal = "APL_DeepSignal-METEORE-perRead-score.tsv.gz"
    apl_megalodon = "APL_Megalodon-METEORE-perRead-score.tsv.gz"

    apl_meteore_comb = "/projects/li-lab/Nanopore_compare/suppdata/METEORE_results/APL.METEORE.megalodon_deepsignal-optimized-model-perRead.combine.tsv.gz"

    apl_deepsignal_set = get_meteore_format_set(os.path.join(bdir, apl_deepsignal))
    apl_megalodon_set = get_meteore_format_set(os.path.join(bdir, apl_megalodon))
    apl_meteore_set = get_meteore_format_set(apl_meteore_comb)

    logger.debug(
        f"Read level: deepsignal={len(apl_deepsignal_set):,}, megalodon={len(apl_megalodon_set):,}, meteore={len(apl_meteore_set):,}")

    apl_deepsignal_set = get_meteore_format_set(os.path.join(bdir, apl_deepsignal), read_level=False)
    apl_megalodon_set = get_meteore_format_set(os.path.join(bdir, apl_megalodon), read_level=False)
    apl_meteore_set = get_meteore_format_set(apl_meteore_comb, read_level=False)
    logger.debug(
        f"Site level: deepsignal={len(apl_deepsignal_set):,}, megalodon={len(apl_megalodon_set):,}, meteore={len(apl_meteore_set):,}")


def sanity_check_merge_bedtools():
    infn = "/projects/li-lab/yang/results/2021-06-30/bed-bk/hg38.gc5Base.bin100.bed.gz"
    bin100_bed = BedTool(infn).sort().merge()
    logger.debug(f"bin100_bed={len(bin100_bed)}")

    infn = "/projects/li-lab/yang/results/2021-06-30/bed-bk/hg38.gc5Base.bin20.bed.gz"
    bin20_bed = BedTool(infn).sort().merge()
    logger.debug(f"bin20_bed={len(bin20_bed)}")

    merge_bin_100_20 = bin100_bed.cat(bin20_bed)
    logger.debug(f"merge_bin_100_20={len(merge_bin_100_20)}")


def sanity_check_get_dna_seq():
    # sanity_check_dna_sequence('chr1', 202108456)
    # sanity_check_dna_sequence('chr1', 202108476)
    #
    # sanity_check_dna_sequence('chr19', 40957177)
    # sanity_check_dna_sequence('chr19', 40958725)

    # sanity_check_dna_sequence('NC_000913.3', 3503572)
    # sanity_check_dna_sequence('NC_000913.3', 3503591)
    #
    # sanity_check_dna_sequence('NC_000913.3', 3507107)
    # sanity_check_dna_sequence('NC_000913.3', 3507083)
    return


def uniqueOption(deduplicate=False):
    """
    return unique option -u if True
    Args:
        deduplicate:

    Returns:

    """
    if deduplicate:
        uniqueOption = '-u'
    else:
        uniqueOption = '  '
    return uniqueOption


def sort_bed_file(infn, outfn, has_header=False, deduplicate=False):
    """
    Sort and save bed files into outfn
    Args:
        infn:
        outfn:

    Returns:
    """
    command = f"zcat {infn} | sort -V {uniqueOption(deduplicate)} -k1,1 -k2,2n | gzip -f > {outfn}"
    if has_header:
        command = f"""
        zcat {infn}  |  awk 'NR<2{{print $0;next}}{{print $0| "sort -V -k1,1 -k2,2n"}}' | gzip -f  > {outfn}
        """
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE) \
        .stdout.read().decode("utf-8")
    logger.debug(f'Sort {infn} and save into {outfn}, has_header={has_header}')


def sort_set_txt_file(infn, outfn, deduplicate=False):
    """
    Sort and save bed files into outfn
    Args:
        infn:
        outfn:

    Returns:
    """
    command = f"zcat {infn} | sort -t '_' -V {uniqueOption(deduplicate)} -k1,1 -k2,2n -k3,3 | gzip -f > {outfn}"
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE) \
        .stdout.read().decode("utf-8")
    logger.debug(f'Sort {infn} and save into {outfn}')
    return True


def sort_per_read_tsv_file(infn, outfn, delimiter="", deduplicate=False):
    """
    Sort and save bed files into outfn
    Args:
        infn:
        outfn:

    Returns:
    """
    command = f"zcat {infn} | sort {delimiter}  -V {uniqueOption(deduplicate)} -k2,2 -k3,3n -k4,4 -k1,1 | gzip -f > {outfn}"
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE) \
        .stdout.read().decode("utf-8")
    logger.debug(f'Sort {infn} and save into {outfn}')
    return True


def sort_per_read_csv_file(infn, outfn, delimiter="-t ,", deduplicate=False):
    """
    Sort and save csv files into outfn, need add -t option for sort
    Args:
        infn:
        outfn:

    Returns:
    """
    return sort_per_read_tsv_file(infn, outfn, delimiter=delimiter, deduplicate=deduplicate)


def convert_size(size_bytes):
    """
    Convert byte number to string
    Args:
        size_bytes:

    Returns:

    """
    if size_bytes == 0:
        return "0B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return "%s %s" % (s, size_name[i])


def get_current_memory_usage():
    """
    Get the memory usage info, ref: https://github.com/dask/distributed/issues/1409
    Returns:

    """
    process = psutil.Process(os.getpid())
    ret = f"VMS:{convert_size(process.memory_info().vms)}, RSS:{convert_size(process.memory_info().rss)}"
    return ret


def freq_to_label(freq, fully_cutoff=1.0, eps=EPSLONG):
    """
    Convert methylation frequency value into 5mC class label
    Args:
        freq:
        fully_cutoff:
        eps:

    Returns:

    """
    if freq <= eps:
        return 0
    if freq >= fully_cutoff - eps:
        return 1
    raise Exception(f"Encounter not fully meth or unmeth value, freq={freq}")


def tool_pred_class_label(log_likelyhood, cutoff=0):
    """
    Infer class label based on log-likelyhood
    Args:
        log_likelyhood:

    Returns:

    """
    if log_likelyhood > cutoff + EPSLONG:
        return 1
    return 0


def load_tool_read_level_unified_as_df(data_file_path, toolname, filterChrSet=None, chunksize=CHUNKSIZE):
    """
    Load read-level unified input
    Args:
        data_file_path:
        toolname:
        filterChrSet:

    Returns:

    """
    logger.debug(f"Load {toolname}:{data_file_path}")

    if filterChrSet is not None:
        iter_df = pd.read_csv(data_file_path, header=0, index_col=False, sep="\t", iterator=True,
                              chunksize=chunksize)
        data_file = pd.concat([chunk[chunk['Chr'].isin(filterChrSet)] for chunk in iter_df])
    else:
        data_file = pd.read_csv(data_file_path, header=0, index_col=False, sep="\t")

    data_file['Pos'] = data_file['Pos'].astype(np.int64)
    data_file.drop_duplicates(subset=['Chr', "ID", "Pos", "Strand"], inplace=True)
    data_file.rename(columns={"Score": toolname}, inplace=True)
    data_file.dropna(inplace=True)
    data_file.reset_index(inplace=True, drop=True)
    return data_file


def prob_to_log(prob):
    """
    convert probability of 5mC to log-likelyhood
    Args:
        prob:

    Returns:

    """
    return math.log2((prob + EPSLONG) / (1 - prob + EPSLONG))


if __name__ == '__main__':
    set_log_debug_level()
    # refGenome = get_ref_fasta('/projects/li-lab/Nanopore_compare/nf_input/reference_genome/ecoli/Ecoli_k12_mg1655.fasta')
    # sanity_check_get_dna_seq()
    # refGenome = get_ref_fasta()
    # sanity_check_get_dna_seq()
    logger.debug("DONE")
