"""
Common functions used by read level and site level evaluations in Nanocompare paper.

Such as import_DeepSignal, import_BGTruth, etc.
"""

import csv
import glob
import gzip
import pickle
import re
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from pybedtools import BedTool
from scipy import stats
from sklearn.metrics import roc_curve, auc, average_precision_score, f1_score, precision_score, recall_score
from tqdm import tqdm

from nanocompare.global_config import *
from nanocompare.global_settings import humanChrSet, ToolEncodeList, BGTruthEncodeList, narrowCoordFileList, narrowCoordFileTag


def importPredictions_Nanopolish(infileName, chr_col=0, start_col=2, strand_col=1, log_lik_ratio_col=5, sequence_col=-1, num_motifs_col=-2, baseFormat=1, llr_cutoff=2.0, output_first=False, include_score=False, filterChr=humanChrSet):
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
    cpgDict = defaultdict(list)
    call_cnt = 0
    meth_cnt = 0
    unmeth_cnt = 0

    infile = open_file_gz_or_txt(infileName)

    for row in infile:
        if isinstance(row, str):
            row1 = row
        else:  # byte array need decode firstly
            row1 = row.decode()

        tmp = row1.strip().split("\t")

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
                if abs(llr) < llr_cutoff * num_sites:  # Consider all sites as a group when there are multiple sites
                    continue

                meth_score = llr
                if llr > 0:
                    meth_indicator = 1

                elif llr < 0:
                    meth_indicator = 0

                strand_info = tmp[strand_col]
                if strand_info == '-':  # - strand, point to positive sequence C, so need +1 to point to G
                    start = start + 1

                if strand_info not in ['-', '+']:
                    raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')
            except:
                logger.error(f'###\tError when parsing row=[{row}] in {infileName}')
                continue

            if num_sites == 1:  # we have singleton, i.e. only one CpG within the area
                if baseFormat == 0:
                    key = (tmp[chr_col], start, strand_info)
                elif baseFormat == 1:
                    key = (tmp[chr_col], start + 1, strand_info)
                else:
                    logger.error("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
                    sys.exit(-1)

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
                        logger.error("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
                        sys.exit(-1)

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
    logger.info(f"###\timportPredictions_Nanopolish SUCCESS: {call_cnt:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file with LLR-cutoff={llr_cutoff:.2f}")
    return cpgDict


def importPredictions_DeepSignal(infileName, chr_col=0, start_col=1, strand_col=2, meth_prob_col=7, meth_col=8, baseFormat=1, include_score=False, filterChr=humanChrSet):
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
    infile = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict(list)
    row_count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        if isinstance(row, str):
            row1 = row
        else:  # byte array need decode firstly
            row1 = row.decode()

        tmp = row1.strip().split("\t")

        if tmp[chr_col] not in filterChr:
            continue

        if baseFormat == 1:
            start = int(tmp[start_col]) + 1
            strand = tmp[strand_col]
        elif baseFormat == 0:
            start = int(tmp[start_col])
            strand = tmp[strand_col]
        else:
            logger.error("###\timportPredictions_DeepSignal InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')
        key = (tmp[chr_col], start, strand)

        if include_score:
            cpgDict[key].append((int(tmp[meth_col]), float(tmp[meth_prob_col])))
        else:
            cpgDict[key].append(int(tmp[meth_col]))

        row_count += 1
        if int(tmp[meth_col]) == 1:
            meth_cnt += 1
        else:
            unmeth_cnt += 1

    infile.close()

    logger.info(f"###\timportPredictions_DeepSignal SUCCESS: {row_count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")

    return cpgDict


def importPredictions_Tombo(infileName, chr_col=0, start_col=1, strand_col=5, meth_col=4, baseFormat=1, cutoff=(-1.5, 2.5), output_first=False, include_score=False, filterChr=humanChrSet):
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
    infile = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict(list)
    row_count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        if isinstance(row, str):
            row1 = row
        else:  # byte array need decode firstly
            row1 = row.decode()

        tmp = row1.strip().split("\t")

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
            logger.error(f"###\timportPredictions_Tombo InputValueError: baseCount value set to '{baseFormat}'. It should be equal to 0 or 1")
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')

        key = (tmp[chr_col], start, strand)

        try:
            methCall = float(tmp[meth_col])
        except Exception as e:
            logger.error(f" ####Tombo parse error at row={row}, exception={e}")
            continue

        meth_score = -methCall
        if methCall < cutoff[0]:  # below -1.5 is methylated by default
            meth_indicator = 1
            meth_cnt += 1
        elif methCall > cutoff[1]:  # above 2.5 is methylated by default
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

    logger.info(f"###\timportPredictions_Tombo SUCCESS: {row_count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-call={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs with meth-cutoff={cutoff} from {infileName} file")
    return cpgDict


def importPredictions_DeepMod(infileName, chr_col=0, start_col=1, strand_col=5, coverage_col=-3, meth_freq_col=-2, meth_cov_col=-1, baseFormat=1, sep=' ', output_first=False, include_score=False, filterChr=humanChrSet, total_cols=13):
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
    infile = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict(list)
    count_calls = 0
    meth_cnt = 0
    unmeth_cnt = 0

    first_row_checked = False

    for row in infile:
        if isinstance(row, str):
            row1 = row
        else:  # byte array need decode firstly
            row1 = row.decode()
        tmp = row1.strip().split(sep)

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
            logger.debug("###\timportPredictions_DeepMod InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')

        key = (tmp[chr_col], start, strand)

        methReads = int(tmp[meth_cov_col])
        coverage = int(tmp[coverage_col])

        meth_freq = methReads / coverage

        if include_score:
            methCallsList = [(1, 1.0)] * methReads + [(0, 0.0)] * (coverage - methReads)
        else:
            methCallsList = [1] * methReads + [0] * (coverage - methReads)

        if key in cpgDict:
            raise Exception(f'In DeepMod.C results, we found duplicate key={key}, this is not correct')

        cpgDict[key] = methCallsList

        count_calls += len(methCallsList)
        meth_cnt += methReads
        unmeth_cnt += coverage - methReads

    infile.close()

    logger.info(f"###\timportPredictions_DeepMod SUCCESS: {count_calls:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")
    return cpgDict


def importPredictions_DeepMod_clustered(infileName, chr_col=0, start_col=1, strand_col=5, coverage_col=-4, meth_cov_col=-2, clustered_meth_freq_col=-1, baseFormat=1, sep=' ', output_first=False, siteLevel=True, include_score=False, filterChr=humanChrSet, total_cols=14):
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

    cpgDict = defaultdict()
    count_calls = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        if isinstance(row, str):
            row1 = row
        else:  # byte array need decode firstly
            row1 = row.decode()
        tmp = row1.strip().split(sep)

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
            logger.error("###\timportPredictions_DeepMod InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
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

    logger.info(f"###\timportDeepMod_clustered Parsing SUCCESS: {count_calls:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")

    return cpgDict


def importPredictions_Megalodon(infileName, chr_col=1, start_col=3, strand_col=2, mod_log_prob_col=4, can_log_prob_col=5, baseFormat=1, cutoff=0.8, sep='\t', output_first=False, include_score=False, filterChr=humanChrSet):
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
    infile = open_file_gz_or_txt(infileName)

    cpgDict = defaultdict(list)
    call_cnt = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        if isinstance(row, str):
            row1 = row
        else:  # byte array need decode firstly
            row1 = row.decode()

        tmp = row1.strip().split(sep)

        if tmp[chr_col] not in filterChr:
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
            logger.error(f"###\timportPredictions_Megalodon_Read_Level InputValueError: baseFormat value set to '{baseFormat}'. It should be equal to 0 or 1")
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')

        key = (tmp[chr_col], start, strand)

        try:
            meth_prob = float(np.e ** float(tmp[mod_log_prob_col]))  # Calculate mod_prob
        except Exception as e:
            logger.error(f" ####Megalodon parse error at row={row}, exception={e}")
            continue

        if meth_prob > cutoff:  ##Keep methylated reads
            meth_indicator = 1
            meth_cnt += 1
        elif meth_prob < 1 - cutoff:  ##Count unmethylated reads
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

    logger.info(f"###\timportPredictions_Megalodon SUCCESS: {call_cnt:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-call={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs with meth-cutoff={cutoff:.2f} from {infileName} file")
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
    for cpg in ontDict:
        if type(ontDict[cpg]) == list:  # value is [0 0 1 1 0 ...]
            if len(ontDict[cpg]) >= minCov:
                result[cpg] = [sum(ontDict[cpg]) / float(len(ontDict[cpg])), len(ontDict[cpg])]
        elif type(ontDict[cpg]) == tuple:  # Used by DeepMod_cluster results, value is (freq, cov) {'freq':0.72, 'cov':18}
            if ontDict[cpg][1] >= minCov:  # no change for site level format, e.g., DeepModCluster
                result[cpg] = ontDict[cpg]
        else:
            raise Exception('Not support type of value of dict')

    logger.info(f"###\treadLevelToSiteLevelWithCov: completed filtering with minCov={minCov}, {len(result):,} CpG sites left for {toolname}")
    return result


def open_file_gz_or_txt(infn):
    """
    Open a txt or gz file, based on suffix
    :param infn:
    :return:
    """
    logger.info(f"open file: {infn}")
    if infn.endswith('.gz'):
        infile = gzip.open(infn, 'rb')
    else:
        infile = open(infn, 'r')
    return infile


# encode format
def importGroundTruth_from_Encode(infileName, chr_col=0, start_col=1, methfreq_col=-1, cov_col=4, strand_col=5, covCutoff=1, baseFormat=1, includeCov=True):
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
    """

    cpgDict = {}

    infile = open_file_gz_or_txt(infileName)

    nrow = 0
    for row in infile:
        nrow += 1

        if isinstance(row, str):
            row1 = row
        else:  # byte array need decode firstly
            row1 = row.decode()

        tmp = row1.strip().split("\t")

        if tmp[chr_col] not in humanChrSet:  # Filter out non-human chrs
            continue

        strand = tmp[strand_col]

        if strand not in ['+', '-']:
            raise Exception(f'input file format error when parsing: {row}, tmp={tmp}')

        try:
            if baseFormat == 1:
                start = int(tmp[start_col]) + 1
            elif baseFormat == 0:
                start = int(tmp[start_col])
            else:
                logger.error("###\timportGroundTruth_BedMethyl_from_Encode InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
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
            cpgDict[key] = [meth_freq, cov_cnt]
        else:
            cpgDict[key] = meth_freq

    infile.close()
    logger.debug(f"###\timportGroundTruth_BedMethyl_from_Encode: loaded information for {len(cpgDict):,} CpGs, with cutoff={covCutoff} ({nrow:,} rows)")

    return cpgDict


# bismark format, Deal with file name like 'bismark_bt2.CpG_report.txt.gz'
def importGroundTruth_from_Bismark(infn, chr_col=0, start_col=1, strand_col=2, meth_col=3, unmeth_col=4, covCutoff=1, baseFormat=1, includeCov=True):
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

    df = df.iloc[:, [chr_col, start_col, df.columns.get_loc("end"), df.columns.get_loc("meth-freq"), df.columns.get_loc("cov"), strand_col]]

    df.columns = ['chr', 'start', 'end', 'meth-freq', 'cov', 'strand']

    cpgDict = defaultdict(list)
    for index, row in df.iterrows():
        chr = row['chr']

        if chr not in humanChrSet:  # Filter out non-human chrs
            continue

        start = int(row['start'])
        strand = row['strand']

        if strand not in ['+', '-']:
            raise Exception(f'strand={strand}  for row={row} is not acceptable, please check use correct function to parse bgtruth file {infn}')

        key = (chr, int(start), strand)
        if key not in cpgDict:
            if includeCov:
                cpgDict[key] = [row['meth-freq'] / 100.0, row['cov']]
            else:
                cpgDict[key] = row['meth-freq'] / 100.0
        else:
            raise Exception(f'In genome-wide, we found duplicate sites: for key={key} in row={row}, please check input file {infn}')

    logger.info(f"###\timportGroundTruth_genome_wide_from_Bismark: loaded {len(cpgDict):,} CpGs with cutoff={covCutoff} from file {infn}")

    return cpgDict


# Will deprecated, before is deal with K562 bed with strand info
def importGroundTruth_bed_file_format(infileName, chr_col=0, start_col=1, meth_col=3, meth_reads_col=4, unmeth_reads_col=5, strand_col=6, covCutt=10, baseFormat=0, chrFilter=False, gzippedInput=True, includeCov=False):
    '''
    We modified this function due to the histogram shows it is 0-based format, NOT 1-based format.

    includeCov  is True if return key->[meth-freq, meth-cov]
                is False if return key-> meth-freq

    ### Description of the columns in this format:

    1. Reference chromosome or scaffold
    2. Start position in chromosome (1-based)  error here
    3. End position in chromosome
    4. methylation percentage (0-100)
    5. methylated reads number
    6. unmethylated reads number

    I think that this output is by default 1-based. This is based on (https://www.bioinformatics.babraham.ac.uk/projects/bismark/):
        bismark2bedGraph: This module does now produce these two output files:
        (1) A bedGraph file, which now contains a header line: 'track type=bedGraph'. The genomic start coords are 0-based, the end coords are 1-based.
        (2) A coverage file ending in .cov. This file replaces the former 'bedGraph --counts' file and is required to proceed with the subsequent step to generate a genome-wide cytosine report (the module doing this has been renamed to coverage2cytosine to reflect this file name change)
    comparison of bedGraph and coverage files suggests that in the latter we deal with 1-based. Also, to get 0-based one have to activate appropriate flag in Bismark, which is not done in MINE, pipeline. I emphasize "mine" as files from somebody else might use different setting, so be careful!

    e.g. structure of the coverage file:
    chr8    205945  205945  100     2       0
    chr8    206310  206310  100     10      0
    chr8    206317  206317  100     10      0
    chr8    206319  206319  90      9       1
    chr8    206322  206322  80      8       2

    ### Output files:
    cpgDict = {"chr\tstart\tend\n" : methylation level (format: float (0-1))} # have all cpgs with specified cutoff (not only 100% 5mC or 100% 5C)
    output coordinates are 1-based genome coordinates.

    Sample gbtruth file is followings:

    gzip -cd /projects/li-lab/yang/workspace/nano-compare/data/bgtruth-data/K562_joined.bed.gz | head
    chr17    60576    60576    0.0    0    1    +
    chr17    61640    61640    0.0    0    1    -
    chr17    61649    61649    100.0    1    0    -
    chr17    61772    61772    100.0    1    0    -
    chr17    61785    61785    100.0    1    0    -
    chr17    61801    61801    0.0    0    1    -
    chr17    61804    61804    100.0    1    0    -
    chr17    62005    62005    100.0    1    0    -
    chr17    62025    62025    100.0    1    0    -
    chr17    62147    62147    100.0    1    0    -


    gzip -cd /pod/2/li-lab/Nanopore_methyl_compare/result/BS_seq_result/HL60_RRBS_ENCFF000MDA.Read_R1.Rep_1_trimmed_bismark_bt2.bismark.cov.gz | head
    chr4	10351	10351	100	2	0
    chr4	10435	10435	100	2	0
    chr4	10443	10443	100	2	0
    chr4	10632	10632	100	6	0
    chr4	10634	10634	100	6	0
    chr4	10643	10643	100	6	0
    chr4	10651	10651	100	6	0
    chr4	11150	11150	100	8	0
    chr4	11155	11155	100	8	0
    chr4	11161	11161	100	8	0

    '''

    cpgDict = {}

    row_cnt = 0

    if gzippedInput:
        infile = gzip.open(infileName, 'rb')
    else:
        infile = open(infileName, 'r')

    for row in infile:
        tmp = row.decode('ascii').strip().split("\t")

        if tmp[chr_col] not in humanChrSet:  # Filter out non-human chrs
            continue

        if baseFormat == 1:
            try:
                start = int(tmp[start_col]) + 1
                end = start
                strand = tmp[strand_col]
            except:
                logger.error(f" ### error when parse ground_truth row={row}")
                continue
        elif baseFormat == 0:
            try:
                start = int(tmp[start_col])
                end = start + 1
                strand = tmp[strand_col]
            except:
                logger.error(f" ### error when parse ground_truth row={row}")
                continue
        else:
            logger.error("###\timportGroundTruth_coverage_output_from_Bismark InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
            sys.exit(-1)

        if chrFilter == False or chrFilter == tmp[chr_col]:
            try:
                temp_meth_and_unmeth = int(tmp[meth_reads_col]) + int(tmp[unmeth_reads_col])
            except:
                logger.error(f" ### Error parse gbTruth row = {row}")
                continue

            row_cnt += 1

            if temp_meth_and_unmeth >= covCutt:
                try:
                    # key = "{}\t{}\t{}\t{}\n".format(tmp[chr_col], start, end, strand)
                    key = (tmp[chr_col], start, strand)
                    if key not in cpgDict:
                        # TODO: add coverage to values also
                        if includeCov:
                            cpgDict[key] = [float(tmp[meth_col]) / 100.0, temp_meth_and_unmeth]
                        else:
                            cpgDict[key] = float(tmp[meth_col]) / 100.0
                    else:
                        logger.error("###\timportGroundTruth_coverage_output_from_Bismark SanityCheckError: One CpG should not have more than 1 entry")
                        sys.exit(-1)
                except:
                    logger.error(f" ### Error parse gbTruth row = {row}")
                    continue

    infile.close()
    logger.info(f"###\timportGroundTruth_coverage_output_from_Bismark: loaded information for {len(cpgDict):,} CpGs with cutoff={covCutt}, before cutoff={row_cnt:,}")
    return cpgDict


# Need proof check before use
def importGroundTruth_oxBS_deprecated(infileName, chr_col='#chromosome', start_col='start', meth_col="pmC", covCutt=4, baseCount=1, chrFilter=False):
    '''
    Note that this function was optimized to parse the data from my oxBS results. More specifically, the (sudo)format which I have created to be able to have the information from both BS and corresponding oxBS-seq.

    ### Description of the columns in this format:
    1. chromosome
    2. start
    3. end
    4. pmC - percentage of methylated
    5. artifact = 0
    6. pC - percentage of unmethylated
    7. artifact = 0
    8. qA - No. of 5mC + 5hmC reads (oxBS-seq based // quadrant A)
    9. qB - No. of C reads (oxBS-seq based // quadrant B)
    10. artifact = 0
    11. artifact = 0
    12. artifact = 0

    e.g.:
    #chromosome	start	end	pmC	phmC	pC	err	qA	qB	qC	qD	N
    chr1	10662	10663	1.0	0.0	0.0	0	4	0	4	0	1.0
    chr1	10665	10666	1.0	0.0	0.0	0	6	0	4	0	1.5
    chr1	10667	10668	1.0	0.0	0.0	0	6	0	4	0	1.5

    ### Output files:
    cpgDict = {"chr\tstart\tend\n" : methylation level (format: float (0-1))} # have all cpgs with specified cutoff (not only 100% 5mC or 100% 5C)
    output coordinates are 1-based genome coordinates.

    '''

    cpgDict = {}

    infile = open(infileName, 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')
    for row in csvfile:
        if baseCount == 1:
            start = int(row[start_col])
        elif baseCount == 0:
            start = int(row[start_col]) + 1
        else:
            logger.error("###\timportGroundTruth_oxBS InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
            sys.exit(-1)

        if chrFilter == False or chrFilter == row[chr_col]:
            cov = int(row["qA"]) + int(row["qB"])
            if cov >= covCutt:
                #                 key = (row[chr_col], start)
                key = "{}\t{}\t{}\n".format(row[chr_col], start, start)
                if key not in cpgDict:
                    cpgDict[key] = float(row['pmC'])
                else:
                    logger.error("###\timportGroundTruth_oxBS SanityCheckError: One CpG should not have more than 1 entry")

    infile.close()
    return cpgDict


# Not used now, if use need proof check
def importGroundTruth_coverage_output_from_Bismark_BedGraph(infileName, chr_col=0, start_col=1, meth_col=3, baseCount=1, gzippedInput=True):
    '''

    We are going to parse followsings:

    2020-01-17 23:26:28,242 - [Methylation_correlation_plotting.py:709] - ERROR:  ### Error parse gbTruth row = b'chr19\t344210\t344211\t3.44827586206897\n'
    2020-01-17 23:26:28,243 - [Methylation_correlation_plotting.py:709] - ERROR:  ### Error parse gbTruth row = b'chr19\t344212\t344213\t0\n'
    2020-01-17 23:26:28,243 - [Methylation_correlation_plotting.py:709] - ERROR:  ### Error parse gbTruth row = b'chr19\t344220\t344221\t4.76190476190476\n'
    2020-01-17 23:26:28,243 - [Methylation_correlation_plotting.py:709] - ERROR:  ### Error parse gbTruth row = b'chr19\t344226\t344227\t0\n'
    2020-01-17 23:26:28,243 - [Methylation_correlation_plotting.py:709] - ERROR:  ### Error parse gbTruth row = b'chr19\t344228\t344229\t9.52380952380952\n'
    2020-01-17 23:26:28,243 - [Methylation_correlation_plotting.py:709] - ERROR:  ### Error parse gbTruth row = b'chr19\t344231\t344232\t0\n'

    ### Description of the columns in this format:

    1. Reference chromosome or scaffold
    2. Start position in chromosome (1-based)
    3. End position in chromosome
    4. methylation percentage (0-100)
    5. methylated reads number
    6. unmethylated reads number

    I think that this output is by default 1-based. This is based on (https://www.bioinformatics.babraham.ac.uk/projects/bismark/):
        bismark2bedGraph: This module does now produce these two output files:
        (1) A bedGraph file, which now contains a header line: 'track type=bedGraph'. The genomic start coords are 0-based, the end coords are 1-based.
        (2) A coverage file ending in .cov. This file replaces the former 'bedGraph --counts' file and is required to proceed with the subsequent step to generate a genome-wide cytosine report (the module doing this has been renamed to coverage2cytosine to reflect this file name change)
    comparison of bedGraph and coverage files suggests that in the latter we deal with 1-based. Also, to get 0-based one have to activate appropriate flag in Bismark, which is not done in MINE, pipeline. I emphasize "mine" as files from somebody else might use different setting, so be careful!

    e.g. structure of the coverage file:
    chr8    205945  205945  100     2       0
    chr8    206310  206310  100     10      0
    chr8    206317  206317  100     10      0
    chr8    206319  206319  90      9       1
    chr8    206322  206322  80      8       2

    ### Output files:
    cpgDict = {"chr\tstart\tend\n" : methylation level (format: float (0-1))} # have all cpgs with specified cutoff (not only 100% 5mC or 100% 5C)
    output coordinates are 1-based genome coordinates.

    '''

    logger.info(f"in importGroundTruth_coverage_output_from_Bismark_BedGraph, infileName={infileName}")

    cpgDict = {}

    if gzippedInput:
        infile = gzip.open(infileName, 'rb')
    else:
        infile = open(infileName, 'r')

    for row in infile:
        tmp = row.decode('ascii').strip().split("\t")
        if baseCount == 1:
            try:
                start = int(tmp[start_col])
            except:
                logger.error(f" ### error when parse ground_truth row={row}")
                continue
        elif baseCount == 0:
            try:
                start = int(tmp[start_col]) + 1
            except:
                logger.error(f" ### error when parse ground_truth row={row}")
                continue
        else:
            logger.error("###\timportGroundTruth_coverage_output_from_Bismark InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
            raise Exception(f'Check error here.')

        try:
            key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)
            if key not in cpgDict:
                cpgDict[key] = float(tmp[meth_col]) / 100.0
            else:
                logger.error("###\timportGroundTruth_coverage_output_from_Bismark SanityCheckError: One CpG should not have more than 1 entry")
        except:
            logger.error(f" ### Error parse gbTruth row = {row}")
            continue

    infile.close()
    logger.info("###\timportGroundTruth_coverage_output_from_Bismark: loaded information for {} CpGs".format(len(cpgDict)))
    return cpgDict


def import_call(infn, encode, baseFormat=1, include_score=False, deepmod_cluster_freq_cov_format=True, enable_cache=False, using_cache=False, filterChr=humanChrSet):
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
    logger.debug(f"Start load {encode}")

    if enable_cache and using_cache:
        ret = check_cache_available(infn=infn, encode=encode, baseFormat=baseFormat, include_score=include_score, deepmod_cluster_freq_cov_format=deepmod_cluster_freq_cov_format)
        if ret:
            logger.debug(f'Import {encode} finished!\n')
            return ret
        logger.debug(f'Not cached yet, we load from raw file')
    call0 = None
    if encode == 'DeepSignal':
        calls0 = importPredictions_DeepSignal(infn, baseFormat=baseFormat, include_score=include_score, filterChr=filterChr)
    elif encode == 'Tombo':
        calls0 = importPredictions_Tombo(infn, baseFormat=baseFormat, include_score=include_score, filterChr=filterChr)
    elif encode == 'Nanopolish':
        calls0 = importPredictions_Nanopolish(infn, baseFormat=baseFormat, llr_cutoff=2.0, include_score=include_score, filterChr=filterChr)
    # elif encode == 'DeepMod':
    #     calls0 = importPredictions_DeepMod_Read_Level(infn, baseFormat=baseFormat, include_score=include_score)
    elif encode == 'Megalodon':  # Original format
        calls0 = importPredictions_Megalodon(infn, baseFormat=baseFormat, include_score=include_score, filterChr=filterChr)
    elif encode == 'Megalodon.ZW':  # Here, ziwei preprocess the raw file to this format: chr_col=0, start_col=1, strand_col=3
        calls0 = importPredictions_Megalodon(infn, baseFormat=baseFormat, include_score=include_score, chr_col=0, start_col=1, strand_col=3, filterChr=filterChr)
    elif encode == 'DeepMod.Cluster':  # import DeepMod Clustered output for site level, itself tool reports by cluster, key->value={'freq':68, 'cov':10}
        calls0 = importPredictions_DeepMod_clustered(infn, baseFormat=baseFormat, siteLevel=deepmod_cluster_freq_cov_format, include_score=include_score, filterChr=filterChr)
    elif encode == 'DeepMod.C':  # import DeepMod itself tool for read level
        calls0 = importPredictions_DeepMod(infn, baseFormat=baseFormat, include_score=include_score, filterChr=filterChr)
    else:
        raise Exception(f'Not support {encode} for file {infn} now')

    if enable_cache:
        save_to_cache(infn, calls0, encode=encode, baseFormat=baseFormat, include_score=include_score, deepmod_cluster_freq_cov_format=deepmod_cluster_freq_cov_format)

    logger.debug(f'Import {encode} finished!\n')
    return calls0


def import_bgtruth(infn, encode, covCutoff=5, baseFormat=1, includeCov=True, enable_cache=False, using_cache=False):
    """
    General purpose for import any BG-Truth input files.

    Import bgtruth from file fn using encode, when use new dataset, MUST check input file start baseFormat and import functions are consistent!!!
    :param infn:
    :param encode:
    :param includeCov:  if true return [freq, cov] as value of each key=(chr, (int)start, strand), or just value=freq. Note: freq is in range of [0,1]
    :return:
    """
    if enable_cache and using_cache:
        ret = check_cache_available(infn, encode=encode, cov=covCutoff, baseFormat=baseFormat, includeCov=includeCov)
        if ret:
            logger.debug(f'Import BG-Truth using encode={encode} finished!\n')
            return ret
        logger.debug(f'Not cached yet, we load from raw file')

    if encode == "encode":  # This is official K562 BS seq data, 0-based start
        bgTruth = importGroundTruth_from_Encode(infn, covCutoff=covCutoff, baseFormat=baseFormat, includeCov=includeCov)
    elif encode == "bismark":  # for genome-wide Bismark Report ouput format results, such as HL60, etc. 1-based start input
        bgTruth = importGroundTruth_from_Bismark(infn, covCutoff=covCutoff, baseFormat=baseFormat, includeCov=includeCov)
    elif encode == "oxBS_sudo":  # not used now, function need to check if use
        raise Exception(f'Please check this function is ok, encode={encode}')
        bgTruth = importGroundTruth_oxBS_deprecated(infn, covCutt=covCutoff, baseCount=baseFormat)
    elif encode == "bed":  # like K562 joined bg-truth, 0-based start input
        raise Exception(f'Please check this function is ok, encode={encode}')
        bgTruth = importGroundTruth_bed_file_format(infn, covCutt=covCutoff, baseFormat=baseFormat, includeCov=includeCov)
    elif encode == "bismark_bedgraph":  # not used now, function need to check if use
        raise Exception(f'Please check this function is ok, encode={encode}')
        bgTruth = importGroundTruth_coverage_output_from_Bismark_BedGraph(infn, baseCount=baseFormat)
    else:
        raise Exception(f"encode={encode} is not supported yet, for inputfile={infn}")

    if enable_cache:
        save_to_cache(infn, bgTruth, encode=encode, cov=covCutoff, baseFormat=baseFormat, includeCov=includeCov)

    logger.debug(f'Import BG-Truth using encode={encode} finished!\n')
    return bgTruth


def calldict2txt(inputDict):
    """
    Convert all keys in dict to a string txt, txt file is:
    chr1  123  123  +

    :param inputDict:
    :return:
    """
    text = ""
    for key in inputDict:
        text += f'{key[0]}\t{key[1]}\t{key[1]}\t{key[2]}\n'
    return text


def txt2dict(pybed, strand_col=3):
    """
    convert bed txt to a dict with keys in bed file

    From
    chr123  123  123 +
    To
    ('chr123', 123, '+') -> 1

    :param pybed:
    :return:
    """
    d = {}
    for t in pybed:
        ret = str(t)[:-1].split('\t')
        key = (ret[0], int(ret[1]), ret[strand_col])
        d[key] = 1
    return d


def load_single_sites_bed_as_set(infn):
    """
    Load the single sites BED files, and return a Set object represents sites.
    :param infn:
    :return:
    """
    infile = open(infn, 'r')
    ret = set()
    for row in infile:  # each row: chr123  123   123  .  .  +
        rowsplit = row[:-1].split('\t')
        key = (rowsplit[0], int(rowsplit[1]), rowsplit[5])
        ret.add(key)
    infile.close()
    return ret


def computePerReadPerfStats(ontCalls, bgTruth, title, coordBedFileName=None, secondFilterBedFileName=None, cutoff_fully_meth=1.0, outdir=None, tagname=None, is_save=True):
    """
    Compute ontCalls with bgTruth performance results by per-read count.
    coordBedFileName        -   full file name of coordinate used to eval
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

    # Firstly reduce ontCalls with in bgTruth keys
    ontCallsKeySet = set(ontCalls.keys()).intersection(set(bgTruth.keys()))
    # newOntCalls = {}
    # for key in bgTruth:
    #     if key in ontCalls:
    #         newOntCalls[key] = ontCalls[key]

    # switch = 0  # 1 if genome-wide
    ontCalls_narrow_set = None  # Intersection of ontCall with coord, or None if genome-wide
    if coordBedFileName:
        # Try ontCall intersect with coord (Genomewide, Singletons, etc.)
        ontCalls_bed = BedTool(calldict2txt(ontCallsKeySet), from_string=True)
        ontCalls_bed = ontCalls_bed.sort()

        coordBed = BedTool(coordBedFileName)
        coordBed = coordBed.sort()
        ontCalls_intersect = ontCalls_bed.intersect(coordBed, u=True, wa=True)
        ontCalls_narrow_set = set(txt2dict(ontCalls_intersect).keys())
    # else:
    #     switch = 1

    ## Second optional filter, organized in the same fashion as the first one. Designed to accomodate for example list of CpGs covered by second program
    # secondSwitch = 0
    ontCalls_narrow_second_set = None  # if using joined sites of all tools, or None for not using joined sites
    if secondFilterBedFileName:
        joined_set = load_single_sites_bed_as_set(secondFilterBedFileName)
        infile = open(secondFilterBedFileName, 'r')
        secondFilterDict = {}
        for row in infile:  # each row: chr123  123   123  .  .  +
            rowsplit = row[:-1].split('\t')
            key = (rowsplit[0], int(rowsplit[1]), rowsplit[5])
            secondFilterDict[key] = 0
        infile.close()
        ontCalls_narrow_second_set = set(ontCalls.keys()).intersection(joined_set)
    # else:
    #     secondSwitch = 1

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
    if ontCalls_narrow_set:
        targetedSet = targetedSet.intersection(ontCalls_narrow_set)

    if ontCalls_narrow_second_set:
        targetedSet = targetedSet.intersection(ontCalls_narrow_second_set)

    for cpgKey in targetedSet:  # key = (chr, start, strand)
        ##### for each sites, we perform per read stats:
        if satisfy_fully_meth_or_unmeth(bgTruth[cpgKey][0]):
            # if bgTruth[cpgKey][0] >= (cutoff_fully_meth - 1e-6) or bgTruth[cpgKey][0] <= 1e-5:  # we only consider absolute states here, 0, or 1 in bgtruth
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
                # totalCalls += 1

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
                yscore_of_ont_tool.append(perCall[1])

                if is_fully_meth(bgTruth[cpgKey][0]):  # BG Truth label
                    y_of_bgtruth.append(1)
                else:
                    y_of_bgtruth.append(0)
        else:
            raise Exception(f'We must see all certain sites here, but see meth_freq={bgTruth[cpgKey][0]}')

    ### compute all per read stats:
    #     Accuracy:
    try:
        accuracy = (TP_5mC + TN_5mC) / float(TP_5mC + FP_5mC + FN_5mC + TN_5mC)
    except ZeroDivisionError:
        accuracy = 0

    #     Positive predictive value (PPV), Precision = (TP) / E(Predicted condition positive)
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

    #     True positive rate (TPR), Recall, Sensitivity, probability of detection = (TP) / (TP+FN)
    try:
        recall_5mC = TP_5mC / float(TP_5mC + FN_5mC)
    except ZeroDivisionError:
        recall_5mC = 0

    try:
        recall_5C = TP_5C / float(TP_5C + FN_5C)
    except ZeroDivisionError:
        recall_5C = 0

    #     F1 score:

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

    ## plot AUC curve:
    # fig = plt.figure(figsize=(5, 5), dpi=300)

    fprSwitch = 1
    try:
        fpr, tpr, _ = roc_curve(y_of_bgtruth, yscore_of_ont_tool)
        average_precision = average_precision_score(y_of_bgtruth, yscore_of_ont_tool)
    except ValueError:
        logger.error(f"###\tERROR for roc_curve: y(Truth):{y_of_bgtruth}, scores(Call pred):{ypred_of_ont_tool}, \nother settings: {title}, {coordBedFileName}, {secondFilterBedFileName}")
        fprSwitch = 0
        roc_auc = 0.0
        average_precision = 0.0

    if fprSwitch == 1:
        roc_auc = auc(fpr, tpr)

    ########################
    basefn = 'x.x.Genome-wide' if not coordBedFileName else os.path.basename(coordBedFileName)

    if is_save:
        # save y and y-pred and y-score for later plot:
        curve_data = {'yTrue': y_of_bgtruth, 'yPred': ypred_of_ont_tool, 'yScore': yscore_of_ont_tool}

        os.makedirs(os.path.join(outdir, 'curve_data'), exist_ok=True)
        outfn = os.path.join(outdir, 'curve_data', f'{tagname}.{basefn}.curve_data.pkl')
        with open(outfn, 'wb') as handle:
            pickle.dump(curve_data, handle)

    return accuracy, roc_auc, average_precision, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, Csites_BGTruth, mCsites_BGTruth


def save_keys_to_single_site_bed(keys, outfn, callBaseFormat=1, outBaseFormat=1, nonstr='.'):
    """
    Save all keys in set of ('chr  123  123  .  .  +\n', etc.) to outfn.
    We use non-string like . in 3rd, 4th columns by BED file format.
    :param keys:
    :param outfn:
    :return:
    """
    outfile = open(outfn, 'w')
    for key in keys:
        if outBaseFormat == 0:
            outfile.write(f'{key[0]}\t{key[1] - callBaseFormat + outBaseFormat}\t{key[1] - callBaseFormat + outBaseFormat + 1}\t{nonstr}\t{nonstr}\t{key[2]}\n')
        else:
            outfile.write(f'{key[0]}\t{key[1] - callBaseFormat + outBaseFormat}\t{key[1] - callBaseFormat + outBaseFormat}\t{nonstr}\t{nonstr}\t{key[2]}\n')
    outfile.close()


def SingletonsAndNonSingletonsScanner(referenceGenomeFile, outfileName_s, outfileName_ns, kbp=10):
    """
    Generate singleton and non-singletons BED file, based on Reference Genome and KBP up and down streams.

    The output file is in 1-based at start coordinate system.
    kbp is up and down k-bp regions and evaluated on positive strand.
    Singletons: only one CpG in the region
    Nonsingletons: more than one CpG in the region
    """
    reference = SeqIO.to_dict(SeqIO.parse(referenceGenomeFile, "fasta"))
    logger.debug(f"###\tSingletonsAndNonSingletonsScanner: {referenceGenomeFile} reference genome file is parsed, up and down bp={kbp}")

    outfile_s = open(outfileName_s, "w")  # "s" stands for Singletons
    outfile_ns = open(outfileName_ns, "w")  # "s" stands for Non-Singletons

    for chromosome in list(reference.keys()):
        if chromosome not in humanChrSet:
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

    logger.debug(f"###\tSingletonsAndNonSingletonsScanner: {referenceGenomeFile} file processed, kbp={kbp}, save to Singletons:{outfile_s}, and Nonsingletons:{outfile_ns}")


# 10-bp regions, current work on it
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
        if abs(evalcpg[0] - cpg[0]) - 2 >= kbp:  # not within kbp regions, skip (Note: kbp is number of X between CG, ex. CGXXCGXXCG)
            continue
        if evalcpg[1] != cpg[1]:  # found there is a CPG not same state, return Discordant
            return False
    return True  # Concordant


def nonSingletonsPostprocessing(absoluteBGTruth, nsRegionsBedFileName, nsConcordantFileName, nsDisCordantFileName, kbp=10, print_first=False):
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
    bedBGTruth = BedTool(calldict2txt(absoluteBGTruth), from_string=True)
    bedBGTruth = bedBGTruth.sort()

    infn = os.path.join(data_base_dir, 'genome-annotation', nsRegionsBedFileName)
    regionNonsingletons = BedTool(infn)
    regionNonsingletons = regionNonsingletons.sort()

    regionWithBGTruth = regionNonsingletons.intersect(bedBGTruth, wa=True, wb=True)  # chr start end   chr start end strand

    regionDict = defaultdict(list)  # key->value, key=region of (chr, start, end), value=list of [f1,f2,etc.] , suche as {regionCoords : [methylation percentage list]}

    is_print_first = print_first
    cntBedLines = 0

    for ovr in regionWithBGTruth:
        cntBedLines += 1
        if is_print_first:
            logger.debug(f'ovr={ovr}')

        regionKey = (ovr[0], int(ovr[1]), int(ovr[2]))  # chr  start  end
        methKey = (ovr[3], int(ovr[4]), ovr[6])  # chr, start, strand

        tmpStart = int(ovr[4])
        tmpStrand = ovr[6]

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
            logger.info(f'regionDict[regionKey]={regionDict[regionKey]}, regionKey={regionKey}')
        is_print_first = False

    logger.info(f'cntBedLines={cntBedLines}')

    is_print_first = print_first
    concordantList = {}  # list of (chr, start)
    discordantList = {}
    meth_cnt_dict = defaultdict(int)  # count meth states in concordant and discordant
    unmeth_cnt_dict = defaultdict(int)
    for region in tqdm(regionDict):  # region is (chr, start, end)
        cpgList = regionDict[region]  # get values as the list of  (start, 0/1, cov)
        chrOut = region[0]

        if is_print_first:
            logger.info(f'cpgList={cpgList}')

        for k in range(len(cpgList)):  # for each CpG
            # Since before we group + - to +, now we assure use start
            startOut = cpgList[k][0]
            methIndicator = cpgList[k][1]
            methCov = cpgList[k][2]

            if eval_concordant_within_kbp_region(cpgList, cpgList[k], kbp=kbp):
                # add this cpg to concordant
                concordantList.update({(chrOut, startOut): (chrOut, startOut, methIndicator, methCov)})  # list of (chr, start)
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
            logger.info(f'concordantList={concordantList}, discordantList={discordantList}')
        is_print_first = False

    concordantSet = list(set(concordantList) - set(discordantList))  # remove duplicate and same in discordant
    outfile_concordant = open(nsConcordantFileName, "w")
    for cpgKey in concordantSet:  # (chrOut, startOut) -> (chrOut, startOut, methIndicator, methCov)
        cpg = concordantList[cpgKey]
        region_txt = '\t'.join([cpg[0], str(cpg[1]), str(cpg[1] + 1), str(cpg[2]), str(cpg[3])]) + '\n'
        outfile_concordant.write(region_txt)
    outfile_concordant.close()

    discordantSet = list(set(discordantList))  # remove duplicate
    outfile_discordant = open(nsDisCordantFileName, "w")
    for cpgKey in discordantSet:
        cpg = discordantList[cpgKey]
        region_txt = '\t'.join([cpg[0], str(cpg[1]), str(cpg[1] + 1), str(cpg[2]), str(cpg[3])]) + '\n'
        outfile_discordant.write(region_txt)
    outfile_discordant.close()

    logger.info(f'save to {[nsConcordantFileName, nsDisCordantFileName]}')
    logger.debug(f'meth_cnt={meth_cnt_dict}, unmeth_cnt={unmeth_cnt_dict}')
    ret = {'Concordant.5mC': meth_cnt_dict['Concordant'], 'Concordant.5C': unmeth_cnt_dict['Concordant'], 'Discordant.5mC': meth_cnt_dict['Discordant'], 'Discordant.5C': unmeth_cnt_dict['Discordant']}
    return ret


def load_nanopolish_df(infn='/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/APL/APL.nanopolish_methylation_calls.tsv'):
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


def load_tombo_df(infn='/projects/li-lab/yang/workspace/nano-compare/data/tools-call-data/K562/K562.tombo_perReadsStats.bed'):
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


# def get_dna_sequence_from_samfile(chr, start, end, bamfile):
#     """
#     Get the specific location DNA sequence at SAM files
#
#     start is 0-based format
#
#     Note: all functions of pysam need to be checked, see https://readthedocs.org/projects/pysam/downloads/pdf/latest/
#
#     :param chr:
#     :param start:
#     :param end:
#     :param bamfile:
#     :return:
#     """
#     for read in bamfile.fetch(chr, start=start, end=end):
#
#         alignedRefPositions = read.query_alignment_start
#         # refStart = alignedRefPositions
#
#         refStart = read.get_reference_positions()[0]
#
#         # refSequence = read.get_reference_sequence()
#         readSequence = read.query_alignment_sequence  # current use
#
#         readSequence = read.query_sequence
#
#         if readSequence is None:  # some read has no sequence, may return None, we only report the first read has sequence
#             continue
#         logger.debug(read.query_name)
#         logger.debug(readSequence)
#
#         # logger.debug(f'ref-start={refStart} len(align) = {len(readSequence)}, len(seq) = {len(readSequence1)} compare two:\n{readSequence}\n{readSequence1}')
#
#         # logger.info(refStart)
#         # logger.info(refSequence)
#         # logger.info(readSequence)
#
#         # logger.debug(readSequence[start - refStart:end - refStart])
#         return readSequence[start - refStart:end - refStart]
#     # if all reads return None, we report None
#     return None
#

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


def get_ref_fasta(ref_fn='/projects/li-lab/Ziwei/Nanopore/data/reference/hg38.fa'):
    ref_fasta = SeqIO.to_dict(SeqIO.parse(open(ref_fn), 'fasta'))
    logger.debug(f'load ref file from {ref_fn}')
    return ref_fasta


def sanity_check_sequence(tlist):
    ret = get_ref_fasta()
    for site in tlist:
        ret_seq = get_dna_base_from_reference('chr8', site, ref_fasta=ret)
        logger.info(f'{site}:{ret_seq}')


#
# def scatter_plot_x_y(xx, yy, tagname='', outdir=None):
#     """
#     :param xx:
#     :param yy:
#     :return:
#     """
#     if outdir is None:
#         outdir = pic_base_dir
#
#     import seaborn as sns
#     import matplotlib.pyplot as plt
#
#     plt.figure(figsize=(4, 4))
#
#     ax = sns.scatterplot(x=xx, y=yy)
#     ax.set_title(f"Scatter of {len(xx):,} CpG sites")
#
#     plt.tight_layout()
#     outfn = os.path.join(outdir, f'scatter-plot-{tagname}.jpg')
#     plt.savefig(outfn, dpi=600, bbox_inches='tight')
#     plt.show()
#     plt.close()
#     logger.info(f"save to {outfn}")
#
#     pass
#
#
# def scatter_plot_cov_compare_df(infn=None, df=None, outdir=None):
#     if df is None:
#         if infn is None:
#             raise Exception("No data source specified.")
#
#         if os.path.splitext(infn)[1] == '.csv':
#             df = pd.read_csv(infn)
#         elif os.path.splitext(infn)[1] == '.pkl':
#             df = pd.read_pickle(infn)
#
#     if outdir is None:
#         outdir = pic_base_dir
#
#     xx = df.iloc[:, 0]
#     yy = df.iloc[:, 1]
#
#     tagname = f'cov-of-{xx.name}-vs-{yy.name}'
#     scatter_plot_x_y(xx, yy, tagname=tagname, outdir=outdir)
#
#
# def scatter_analysis_cov(Tombo_calls1, Nanopolish_calls1, outdir, RunPrefix, tool1_name='Tombo', tool2_name='nanopolish'):
#     tomboCpGs = set(Tombo_calls1.keys())
#     intersectCpGs = tomboCpGs.intersection(Nanopolish_calls1.keys())
#     key_list = []
#     tombo_list = []
#     tool_list = []
#     for cpg in intersectCpGs:
#         key_list.append(cpg[:-1])
#         tombo_list.append(len(Tombo_calls1[cpg]))
#         tool_list.append(len(Nanopolish_calls1[cpg]))
#
#     scatterDf = pd.DataFrame(data={f'{tool1_name}-cov': tombo_list, f'{tool2_name}-cov': tool_list}, index=key_list)
#     outfn = os.path.join(outdir, f'{RunPrefix}-{tool1_name}-{tool2_name}-scatter.csv')
#     scatterDf.to_csv(outfn)
#
#     outfn = os.path.join(outdir, f'{RunPrefix}-{tool1_name}-{tool2_name}-scatter.pkl')
#     scatterDf.to_pickle(outfn)
#
#     scatter_plot_cov_compare_df(df=scatterDf, outdir=outdir)


def get_cache_filename(infn, params):
    basefn = os.path.basename(infn)
    cachefn = f'cachefile.{basefn}.encode.{params["encode"]}.base.{params["baseFormat"]}'

    if params["encode"] in ToolEncodeList:
        cachefn += f'.inscore.{params["include_score"]}'
        if params["encode"] == 'DeepMod.Cluster':
            cachefn += f'.deepmod_cluster_cov.{params["deepmod_cluster_freq_cov_format"]}'
    elif params["encode"] in BGTruthEncodeList:
        cachefn += f'.cov.{params["cov"]}.incov.{params["includeCov"]}'
    else:
        raise Exception(f'Encode {params["encode"]} is not support now')
    cachefn = os.path.join(cache_dir, cachefn + '.pkl')
    return cachefn


def save_to_cache(infn, data, **params):
    if not data:
        return
    # logger.debug(f'infn={infn}, data={len(data)}, params={params}')
    cache_fn = get_cache_filename(infn, params)

    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir, exist_ok=True)
    with open(cache_fn, 'wb') as outf:
        pickle.dump(data, outf)
    logger.debug(f'Cached to file:{cache_fn}')


def check_cache_available(infn, **params):
    cachefn = get_cache_filename(infn, params)

    if os.path.exists(cachefn):
        logger.debug(f'Start to get from cache:{cachefn}')
        try:
            with open(cachefn, 'rb') as inf:
                ret = pickle.load(inf)
                logger.info(f'Get from cache:{cachefn}')
            return ret
        except:
            return None
    return None


def filter_cpg_dict(cpgDict, filterDict):
    """
    Filter and keep only filterDict keys
    :param cpgDict:
    :param filterDict:
    :return:
    """
    retDict = defaultdict(list)
    joinedKeys = set(filterDict.keys()).intersection(set(cpgDict.keys()))
    for k in joinedKeys:
        retDict[k] = cpgDict[k]
    return retDict


# def compare_cpg_key(item1, item2):
#     """
#     First compare chr, then start number
#     usage:
#     k1 = ('chr1', 123, '+')
#     k2 = ('chr1', 124, '+')
#     k3 = ('chr3', 1, '-')
#     k4 = ('chr1', 124, '-')
#     k5 = ('chr3', 1, '-')
#     l = [k3, k1, k2, k4, k5]
#     import functools
#     l_sorted = sorted(l, key=functools.cmp_to_key(compare_cpg_key))
#     logger.info(l_sorted)
#
#     output:
#     [('chr1', 123, '+'), ('chr1', 124, '+'), ('chr1', 124, '-'), ('chr3', 1, '-'), ('chr3', 1, '-')]
#
#     :param item1:
#     :param item2:
#     :return:
#     """
#     chr1 = int(item1[0].replace('chr', '').replace('X', '23').replace('Y', '24'))
#     chr2 = int(item2[0].replace('chr', '').replace('X', '23').replace('Y', '24'))
#
#     # return chr1 < chr2 or (chr1 == chr2 and item1[1] < item2[1])
#
#     if chr1 != chr2:
#         return chr1 - chr2
#     elif item1[1] != item2[1]:
#         return item1[1] - item2[1]
#     elif item1[2][0] != item2[2][0]:
#         return ord(item1[2][0]) - ord(item2[2][0])
#     return 0


def is_fully_meth(methfreq, eps=1e-5, cutoff_fully_meth=1.0):
    if methfreq > 1.0 + eps or methfreq < 0:
        raise Exception(f'detect non-freq value for freq={methfreq}')

    if methfreq > cutoff_fully_meth - eps:  # near 1
        return True
    return False


def is_fully_unmeth(methfreq, eps=1e-5):
    if methfreq > 1.0 + eps or methfreq < 0:
        raise Exception(f'detect non-freq value for freq={methfreq}')

    if methfreq < eps:  # near 0
        return True
    return False


def satisfy_fully_meth_or_unmeth(methfreq, eps=1e-5, cutoff_fully_meth=1.0):
    """
    Return true if fully meth or unmeth, eps is a near number of 0 and 1
    :param methfreq:
    :return:
    """
    if is_fully_meth(methfreq, eps=eps, cutoff_fully_meth=cutoff_fully_meth) or is_fully_unmeth(methfreq, eps=eps):
        return True
    return False


def combineBGTruthList(bgTruthList, covCutoff=1):
    """
    Combine two replicates together, we joined two replicates together as one bgtruth, and retain only cov >= covCutoff sites
    :param bgTruthList:
    :return:
    """
    unionBGTruth = {}  # used for singleton and non-singletons detect, used for performance eval
    jointBGTruth = {}  # intersection of CpG sites
    if len(bgTruthList) == 2:  # sites must in both replicates, and
        logger.info(f'Start study union of multiple BG-Truth, with coverage={covCutoff}')
        unionSet = set(bgTruthList[0].keys()).union(set(bgTruthList[1].keys()))
        jointSet = set(bgTruthList[0].keys()).intersection(set(bgTruthList[1].keys()))

        for key in unionSet:
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
                raise Exception(f"Compute joined meth_freq={meth_freq} error, using meth1={meth1}, cov1={cov1}, meth2={meth2}, cov2={cov2}, in sites key={key}, and bgTruthList[0][key]={bgTruthList[0][key]}, bgTruthList[1][key]={bgTruthList[1][key]}")

            cov = int(cov1 + cov2)

            if cov < covCutoff:
                continue

            unionBGTruth[key] = [meth_freq, cov]
            if key in jointSet:
                jointBGTruth[key] = [meth_freq, cov]
        logger.info(f'unionBGTruth = {len(unionBGTruth):,}, jointBGTruth={len(jointBGTruth):,}, with cov-cutoff={covCutoff}')
    elif len(bgTruthList) == 1:
        for key in bgTruthList[0]:
            if bgTruthList[0][key][1] >= covCutoff:
                unionBGTruth[key] = bgTruthList[0][key]
        logger.info(f'Only 1 replicates, BGTruth = {len(unionBGTruth):,}, with cov-cutoff={covCutoff}')
    else:
        raise Exception(f'len={len(bgTruthList)}, is not support now.')

    return unionBGTruth


def filter_cpgkeys_using_bedfile(cpgKeys, bedFileName):
    """
    Keep only cpg keys in bed file range, return set of keys
    :param cpgKeys:
    :param bedFileName:
    :return:
    """
    cpgBed = BedTool(calldict2txt(cpgKeys), from_string=True)
    cpgBed = cpgBed.sort()

    coordBed = BedTool(bedFileName)
    coordBed = coordBed.sort()

    intersectBed = cpgBed.intersect(coordBed, u=True, wa=True)
    ret = set(txt2dict(intersectBed).keys())
    return ret


def find_bed_filename(basedir, pattern):
    """
    Find BED file of concordant and discordant for each dataset or run.
    :param basedir:
    :param pattern:
    :return:
    """
    fnlist = glob.glob(os.path.join(basedir, '**', pattern), recursive=True)
    # logger.info(fnlist)
    if len(fnlist) != 1:
        raise Exception(f'Find more files: {fnlist}, please check the basedir is correct')
    logger.debug(f'find bed file:{fnlist[0]}')
    return fnlist[0]


def gen_venn_data(set_dict, namelist, outdir, tagname='tagname'):
    """
    Generate 7 data for three set or 31 data for five set joining Venn Diagram plotting
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
        raise Exception(f'Number of set need to have combinations corresponding, but get cnt={len(retlist)} for set number={len(namelist)}, code bugs.')

    outfn = os.path.join(outdir, f'venn.data.{tagname}.dat')
    with open(outfn, 'w') as outf:
        for num in retlist:
            outf.write(f'{num}\n')
    logger.info(f'Note the venn data set, tool-name order must be: {namelist}')
    logger.info(f'save {len(retlist)} points venn data for {tagname} to {outfn}')


def filter_corrdata_df_by_bedfile(df, df_bed, coord_fn):
    """
    Filter lines in correlation data, within coordinate BED file
    :param df:
    :param coord_fn:
    :return:
    """
    if not coord_fn:  # None means genome wide
        return df

    coordBed = BedTool(coord_fn)
    coordBed = coordBed.sort()
    df_bed_intersect = df_bed.intersect(coordBed, u=True, wa=True)

    select_lines = []
    for line in df_bed_intersect:
        select_lines.append(str(line)[:-1])
        # logger.info(select_lines)
    df = df[df['bedline'].isin(select_lines)]
    # logger.info(df)
    return df

    pass


def correlation_report_on_regions(corr_infn, beddir=None, dsname=None, outdir=pic_base_dir):
    df = pd.read_csv(corr_infn)
    logger.info(df)

    location_flist = list(narrowCoordFileList)
    location_ftag = list(narrowCoordFileTag)

    if beddir:
        concordantFileName = find_bed_filename(basedir=beddir, pattern=f'{dsname}*hg38_nonsingletons.concordant.bed')

        discordantFileName = find_bed_filename(basedir=beddir, pattern=f'{dsname}*hg38_nonsingletons.discordant.bed')
        location_flist.extend([concordantFileName, discordantFileName])
        location_ftag.extend(['Concordant', 'Discordant'])
    # logger.info(location_ftag)
    # logger.info(location_flist)

    df['bedline'] = df["chr"] + '\t' + df["start"].astype(str) + '\t' + df["end"].astype(str) + '\t' + df["strand"]  # df[['chr', 'start', 'end']].agg('\t'.join, axis=1)
    bedline_str = '\n'.join(df['bedline'].tolist())
    df_bed = BedTool(bedline_str, from_string=True)
    df_bed = df_bed.sort()

    dataset = defaultdict(list)
    for tagname, coord_fn in zip(location_ftag[:], location_flist[:]):
        logger.info(f'tagname={tagname}, coord_fn={coord_fn}')
        newdf = filter_corrdata_df_by_bedfile(df, df_bed, coord_fn)

        # Computer COE and pvalue
        newdf = newdf.filter(regex='_freq$', axis=1)
        for i in range(1, len(newdf.columns)):
            toolname = str(newdf.columns[i]).replace('_freq', '')
            try:  # too few samples will fail
                coe, pval = stats.pearsonr(newdf.iloc[:, 0], newdf.iloc[:, i])
            except:
                coe, pval = None, None

            # report to dataset
            dataset['dsname'].append(dsname)
            dataset['Tool'].append(toolname)
            dataset['Location'].append(tagname)
            dataset['#Bases'].append(len(newdf))
            dataset['COE'].append(coe)
            dataset['p-value'].append(pval)

    # logger.info(dataset)
    outdf = pd.DataFrame.from_dict(dataset)
    logger.info(outdf)

    outfn = os.path.join(outdir, f'{dsname}.corrdata.coe.pvalue.each.regions.xlsx')
    outdf.to_excel(outfn)
    logger.info(f'save to {outfn}')

    return outdf


if __name__ == '__main__':
    set_log_debug_level()
    ## check branch only

    # outs = os.path.join(pic_base_dir, 'hg38_singletons_5bp.bed')
    # outns = os.path.join(pic_base_dir, 'hg38_nonsingletons_5bp.bed')
    # SingletonsAndNonSingletonsScanner(referenceGenomeFile=referenceGenomeFile, outfileName_s=outs, outfileName_ns=outns, kbp=5)
    #
    # outs = os.path.join(pic_base_dir, 'hg38_singletons_10bp.bed')
    # outns = os.path.join(pic_base_dir, 'hg38_nonsingletons_10bp.bed')
    # SingletonsAndNonSingletonsScanner(referenceGenomeFile=referenceGenomeFile, outfileName_s=outs, outfileName_ns=outns, kbp=10)
    import os

    infn = os.path.join("/projects/li-lab/yang/workspace/nano-compare/outputs/TestData-methylation-callings", "TestData.DeepSignal.combine.tsv.gz")

    ontDeepSignal = importPredictions_DeepSignal(infn)

    logger.info("DONE")
