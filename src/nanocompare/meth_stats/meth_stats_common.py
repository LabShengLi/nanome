"""
Common functions used by Methylation_correlation_plotting.py and Universal_meth_stats_evaluation.py

Such as import_DeepSignal, import_BGTruth, etc.
"""

import csv
import gzip
import itertools
import pickle
import re
from collections import defaultdict
from inspect import signature

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from pybedtools import BedTool
from scipy.stats import pearsonr
from sklearn.metrics import roc_curve, precision_recall_curve, auc, confusion_matrix, average_precision_score, f1_score, precision_score, recall_score
from tqdm import tqdm

from nanocompare.global_config import *
from nanocompare.global_settings import humanChrs


def report2dict(cr):
    # solution of  "jolespin commented on Jan 19, 2017", + the updated by "HyungSeokPark commented on Mar 20, 2018"
    # https://github.com/scikit-learn/scikit-learn/issues/7845
    # i needed that because classification_report function does not recognize the "output_dict" parameter

    # Parse rows
    tmp = list()
    for row in cr.split("\n"):
        parsed_row = [x for x in row.split("  ") if len(x) > 0]
        if len(parsed_row) > 0:
            tmp.append(parsed_row)

    # Store in dictionary
    measures = tmp[0]
    logger.debug(measures)

    D_class_data = {}  # defaultdict(dict)
    for row in tmp[1:]:
        class_label = row[0].strip()
        for j, m in enumerate(measures):
            D_class_data[class_label][m.strip()] = float(row[j + 1].strip())
    return D_class_data


def importPredictions_NanoXGBoost(infileName, chr_col=0, start_col=1, meth_col=4, baseFormat=1):
    '''
    Note that the function requires per read stats, not frequencies of methylation.
    !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!

    ### Parameters of the function:
    chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
    start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
    meth_col - name (as header) of the column with methylation call (integer expected). If "header" variable == False, give integer number of the column.
    baseCount - 0 or 1, standing for 0-based or 1-based, respectively
    #header - True or False.

    ### Example input format from my NanoXGBoost model:
    chr20   15932019        15932019        /home/rosikw/fastscratch/APL_newSept/0/GXB01186_20180508_FAH83098_GA10000_sequencing_run_180508_18_li_001_GXB01186_001_12562_read_35498_ch_259_strand.fast5     1
    chr20   15932898        15932898        /home/rosikw/fastscratch/APL_newSept/0/GXB01186_20180508_FAH83098_GA10000_sequencing_run_180508_18_li_001_GXB01186_001_12562_read_35498_ch_259_strand.fast5     0
    chr20   15932997        15932997        /home/rosikw/fastscratch/APL_newSept/0/GXB01186_20180508_FAH83098_GA10000_sequencing_run_180508_18_li_001_GXB01186_001_12562_read_35498_ch_259_strand.fast5     1
    chr20   15933019        15933019        /home/rosikw/fastscratch/APL_newSept/0/GXB01186_20180508_FAH83098_GA10000_sequencing_run_180508_18_li_001_GXB01186_001_12562_read_35498_ch_259_strand.fast5     1
    Note: here there are no headers (probably this will change in the future)

    ### Example input format from Nanopolish:
    chromosome      start   end     read_name       log_lik_ratio   log_lik_methylated      log_lik_unmethylated    num_calling_strands     num_cpgs        sequence
    chr20   106142  106142  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    2.32    -93.28  -95.61  1       1       CTCAACGTTTG
    chr20   106226  106226  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    1.46    -193.91 -195.37 1       1       TGGCACGTGGA
    chr20   104859  104859  107a0850-7500-443c-911f-4857424c889c    4.36    -163.26 -167.62 1       1       ATTCCCGAGAG
    Note: Nanopolish output have header, yet the function needs to be universal enough to handle both cases (including whatever awaits in case of NanoMod, DeepSignal etc.)

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls]}
    output coordinates are 1-based genome coordinates.

    ============

    Future changes:
    This function cannot be for everything. Output from nanopolish is too different from mine (e.g. singletons and non-singletons), so it will need an independent parser. Maybe the same will go for other programs?
    Nevertheless, all functions will have the same output, which is the most important part.
    '''

    infile = open(infileName, "r")
    cpgDict = {}
    count = 0

    for row in infile:
        tmp = row.strip().split("\t")
        if baseFormat == 1:
            start = int(tmp[start_col])
        elif baseFormat == 0:
            start = int(tmp[start_col]) + 1
        else:
            logger.error("###\timportPredictions_NanoXGBoost InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
            sys.exit(-1)
        #         key = (tmp[chr_col], start)
        key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)
        if key not in cpgDict:
            cpgDict[key] = []
        cpgDict[key].append(int(tmp[meth_col]))
        count += 1

    infile.close()

    logger.debug("###\timportPredictions_NanoXGBoost SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
    return cpgDict


def importPredictions_Nanopolish(infileName, chr_col=0, start_col=2, strand_col=1, log_lik_ratio_col=5, sequence_col=-1, num_motifs_col=-2, baseFormat=0, llr_cutoff=2.0, output_first=False):
    """
    We assume the input is 0-based for the start col, such as chr10 122 122

    Return dict of key='chr1\t123\t123\t+', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.

    if strand-info=='-', the start still point to positive sequence CG's C position, should be +1 to point to G. We clear the bugs in followings for correct CG, and for reverse strand GC, we need add 1 more to start. Because the sequence is still positive strand sequence, so report position should be at G of CG.

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

    Their current script for handling with conversion of calls to frequencies:
    https://github.com/jts/nanopolish/blob/master/scripts/calculate_methylation_frequency.py

    It looks that we both do the same, or not?
    """
    cpgDict = defaultdict(list)
    count = 0
    meth_cnt = 0
    unmeth_cnt = 0
    infile = open(infileName, 'r')

    for row in infile:
        tmp = row.strip().split("\t")
        if tmp[chr_col] != "chromosome":
            if output_first:
                logger.debug(list(enumerate(tmp)))
                output_first = False

            if tmp[chr_col] not in humanChrs:
                continue

            try:  # try to find if these columns are interpretable
                start = int(tmp[start_col])
                num_sites = int(tmp[num_motifs_col])
                llr = float(tmp[log_lik_ratio_col])
                if abs(llr) < llr_cutoff:
                    continue

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

                cpgDict[key].append(meth_indicator)
                count += 1
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

                    cpgDict[key].append(meth_indicator)
                    count += 1
                    if meth_indicator == 1:
                        meth_cnt += 1
                    else:
                        unmeth_cnt += 1

    infile.close()
    logger.info(f"###\timportPredictions_Nanopolish SUCCESS: {count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file with LLR-cutoff={llr_cutoff:.2f}")
    return cpgDict


def importPredictions_Nanopolish_v2(infileName, baseCount=0, logLikehoodCutt=2.5, IncludeNonSingletons=True):
    """
    Nanopolish Parser function based on parsing script from Nanopolish:
    https://github.com/jts/nanopolish/blob/master/scripts/calculate_methylation_frequency.py
    Code was downloaded from their Github and modified for the purpose of this project on April 30th 2019.

    Generally it gives exactly the same results as my own, but i think its better to use their code, so that nobody would be able to say that we did something differently


    !!! This function will be needed for NanoCompare project !!!

    ### Example input format from Nanopolish:
    chromosome      start   end     read_name       log_lik_ratio   log_lik_methylated      log_lik_unmethylated    num_calling_strands     num_cpgs        sequence
    chr20   106142  106142  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    2.32    -93.28  -95.61  1       1       CTCAACGTTTG
    chr20   106226  106226  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    1.46    -193.91 -195.37 1       1       TGGCACGTGGA
    chr20   104859  104859  107a0850-7500-443c-911f-4857424c889c    4.36    -163.26 -167.62 1       1       ATTCCCGAGAG

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls]}
    output coordinates are 1-based genome coordinates.
    """
    count = 0
    cpgDict = {}

    infile = open(infileName, 'r')
    csv_reader = csv.DictReader(infile, delimiter='\t')

    for record in csv_reader:
        # print(record)
        # num_sites = int(record['num_cpgs'])
        # llr = float(record['log_lik_ratio'])

        try:
            num_sites = int(record['num_cpgs'])
            llr = float(record['log_lik_ratio'])
        except:  # skip not parsed results
            logger.error(f"Can not parse record: record={record}")
            continue

        # Skip ambiguous call
        if abs(llr) < logLikehoodCutt:
            continue
        # sequence = record['sequence']

        is_methylated = int(llr > 0)

        # if this is a multi-cpg group and split_groups is set, break up these sites
        if IncludeNonSingletons and num_sites > 1:
            c = str(record['chromosome'])
            s = int(record['start'])
            e = int(record['end'])

            # find the position of the first CG dinucleotide
            sequence = record['sequence']
            cg_pos = sequence.find("CG")
            first_cg_pos = cg_pos
            while cg_pos != -1:
                #                 key = (c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                if baseCount == 0:
                    key = "{}\t{}\t{}\n".format(c, s + cg_pos - first_cg_pos + 1, s + cg_pos - first_cg_pos + 1)
                elif baseCount == 1:
                    key = "{}\t{}\t{}\n".format(c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                else:
                    logger.error("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
                    sys.exit(-1)

                if key not in cpgDict:
                    cpgDict[key] = []
                cpgDict[key].append(is_methylated)
                count += 1

                cg_pos = sequence.find("CG", cg_pos + 1)
        else:
            if baseCount == 0:
                key = "{}\t{}\t{}\n".format(str(record['chromosome']), int(record['start']) + 1, int(record['end']) + 1)
            elif baseCount == 1:
                key = "{}\t{}\t{}\n".format(str(record['chromosome']), int(record['start']), int(record['end']))
            else:
                logger.error("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
                sys.exit(-1)

            if key not in cpgDict:
                cpgDict[key] = []
            cpgDict[key].append(is_methylated)
            count += 1

    logger.info("###\timportPredictions_Nanopolish SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
    return cpgDict


def importPredictions_Nanopolish_2_nofilter(infileName, baseCount=0, logLikehoodCutt=2.5, IncludeNonSingletons=True):
    '''
    Not filter any CpG sites
    Nanopolish Parser function based on parsing script from Nanopolish:
    https://github.com/jts/nanopolish/blob/master/scripts/calculate_methylation_frequency.py
    Code was downloaded from their Github and modified for the purpose of this project on April 30th 2019.

    Generally it gives exactly the same results as my own, but i think its better to use their code, so that nobody would be able to say that we did something differently


    !!! This function will be needed for NanoCompare project !!!

    ### Example input format from Nanopolish:
    chromosome      start   end     read_name       log_lik_ratio   log_lik_methylated      log_lik_unmethylated    num_calling_strands     num_cpgs        sequence
    chr20   106142  106142  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    2.32    -93.28  -95.61  1       1       CTCAACGTTTG
    chr20   106226  106226  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    1.46    -193.91 -195.37 1       1       TGGCACGTGGA
    chr20   104859  104859  107a0850-7500-443c-911f-4857424c889c    4.36    -163.26 -167.62 1       1       ATTCCCGAGAG

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls]}
    output coordinates are 1-based genome coordinates.



    '''
    count = 0
    cpgDict = {}

    infile = open(infileName, 'r')
    csv_reader = csv.DictReader(infile, delimiter='\t')

    for record in csv_reader:
        # print(record)
        # num_sites = int(record['num_cpgs'])
        # llr = float(record['log_lik_ratio'])

        try:
            num_sites = int(record['num_cpgs'])
            llr = float(record['log_lik_ratio'])
        except:  # skip not parsed results
            logger.error(f"Can not parse record: record={record}")
            continue

        # Skip ambiguous call
        if abs(llr) < logLikehoodCutt:
            is_methylated = -1
        else:
            is_methylated = int(llr > 0)
        # sequence = record['sequence']

        # if this is a multi-cpg group and split_groups is set, break up these sites
        if IncludeNonSingletons and num_sites > 1:
            c = str(record['chromosome'])
            s = int(record['start'])
            e = int(record['end'])

            # find the position of the first CG dinucleotide
            sequence = record['sequence']
            cg_pos = sequence.find("CG")
            first_cg_pos = cg_pos
            while cg_pos != -1:
                #                 key = (c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                if baseCount == 0:
                    key = "{}\t{}\t{}\n".format(c, s + cg_pos - first_cg_pos + 1, s + cg_pos - first_cg_pos + 1)
                elif baseCount == 1:
                    key = "{}\t{}\t{}\n".format(c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                else:
                    logger.error("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
                    sys.exit(-1)

                if key not in cpgDict:
                    cpgDict[key] = []
                cpgDict[key].append(is_methylated)
                count += 1

                cg_pos = sequence.find("CG", cg_pos + 1)
        else:
            if baseCount == 0:
                key = "{}\t{}\t{}\n".format(str(record['chromosome']), int(record['start']) + 1, int(record['end']) + 1)
            elif baseCount == 1:
                key = "{}\t{}\t{}\n".format(str(record['chromosome']), int(record['start']), int(record['end']))
            else:
                logger.error("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
                sys.exit(-1)

            if key not in cpgDict:
                cpgDict[key] = []
            cpgDict[key].append(is_methylated)
            count += 1

    logger.info("###\timportPredictions_Nanopolish nofilter SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
    return cpgDict


def importPredictions_Nanopolish_3(infileName, baseCount=0, logLikehoodCutt=2.5, IncludeNonSingletons=True):
    '''
    Score import for AUC

    Nanopolish Parser function based on parsing script from Nanopolish:
    https://github.com/jts/nanopolish/blob/master/scripts/calculate_methylation_frequency.py
    Code was downloaded from their Github and modified for the purpose of this project on April 30th 2019.

    Generally it gives exactly the same results as my own, but i think its better to use their code, so that nobody would be able to say that we did something differently


    !!! This function will be needed for NanoCompare project !!!

    ### Example input format from Nanopolish:
    chromosome      start   end     read_name       log_lik_ratio   log_lik_methylated      log_lik_unmethylated    num_calling_strands     num_cpgs        sequence
    chr20   106142  106142  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    2.32    -93.28  -95.61  1       1       CTCAACGTTTG
    chr20   106226  106226  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    1.46    -193.91 -195.37 1       1       TGGCACGTGGA
    chr20   104859  104859  107a0850-7500-443c-911f-4857424c889c    4.36    -163.26 -167.62 1       1       ATTCCCGAGAG

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls]}
    output coordinates are 1-based genome coordinates.



    '''
    count = 0
    cpgDict = {}

    infile = open(infileName, 'r')
    csv_reader = csv.DictReader(infile, delimiter='\t')

    for record in csv_reader:
        # print(record)
        # num_sites = int(record['num_cpgs'])
        # llr = float(record['log_lik_ratio'])

        try:
            num_sites = int(record['num_cpgs'])
            llr = float(record['log_lik_ratio'])
        except:  # skip not parsed results
            logger.error(f"Can not parse record: record={record}")
            continue

        # Skip ambiguous call
        if abs(llr) < logLikehoodCutt:
            continue
        # sequence = record['sequence']

        # is_methylated = int(llr > 0)
        is_methylated = llr

        # if this is a multi-cpg group and split_groups is set, break up these sites
        if IncludeNonSingletons and num_sites > 1:
            c = str(record['chromosome'])
            s = int(record['start'])
            e = int(record['end'])

            # find the position of the first CG dinucleotide
            sequence = record['sequence']
            cg_pos = sequence.find("CG")
            first_cg_pos = cg_pos
            while cg_pos != -1:
                #                 key = (c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                if baseCount == 0:
                    key = "{}\t{}\t{}\n".format(c, s + cg_pos - first_cg_pos + 1, s + cg_pos - first_cg_pos + 1)
                elif baseCount == 1:
                    key = "{}\t{}\t{}\n".format(c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                else:
                    logger.error("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
                    sys.exit(-1)

                if key not in cpgDict:
                    cpgDict[key] = []
                cpgDict[key].append(is_methylated)
                count += 1

                cg_pos = sequence.find("CG", cg_pos + 1)
        else:
            if baseCount == 0:
                key = "{}\t{}\t{}\n".format(str(record['chromosome']), int(record['start']) + 1, int(record['end']) + 1)
            elif baseCount == 1:
                key = "{}\t{}\t{}\n".format(str(record['chromosome']), int(record['start']), int(record['end']))
            else:
                logger.error("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
                sys.exit(-1)

            if key not in cpgDict:
                cpgDict[key] = []
            cpgDict[key].append(is_methylated)
            count += 1

    logger.info("###\timportPredictions_Nanopolish SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
    return cpgDict


def importPredictions_DeepSignal(infileName, chr_col=0, start_col=1, strand_col=2, meth_col=8, baseFormat=0):
    '''
    We treat input as 0-based format for start col.

    Return dict of key='chr1\t123\t123\t+', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.

    Note that the function requires per read stats, not frequencies of methylation.

    ### Parameters of the function:
    chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
    start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
    meth_col - name (as header) of the column with methylation call (integer expected). If "header" variable == False, give integer number of the column. *7 is probability, 8 is binary call.
    baseCount - 0 or 1, standing for 0-based or 1-based, respectively
    #header - True or False.

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

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
    output coordinates are 1-based genome coordinates.

    ** by default if this probability will be higher than 0.5, DeppSignal will tell that this is methylated site, lower, unmethylated

    ============

    '''

    infile = open(infileName, "r")
    cpgDict = defaultdict(list)
    row_count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        tmp = row.strip().split("\t")

        if tmp[chr_col] not in humanChrs:
            continue

        if baseFormat == 1:
            start = int(tmp[start_col]) + 1
            end = start
            strand = tmp[strand_col]
        elif baseFormat == 0:
            start = int(tmp[start_col])
            end = start + 1
            strand = tmp[strand_col]
        else:
            logger.error("###\timportPredictions_DeepSignal InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')
        #         key = (tmp[chr_col], start)
        # key = "{}\t{}\t{}\t{}\n".format(tmp[chr_col], start, end, strand)
        key = (tmp[chr_col], start, strand)
        # if key not in cpgDict:
        #     cpgDict[key] = []
        #         cpgDict[key].append(float(tmp[meth_col])) ##### uncomment this line to get probabilities instead of final, binary calls
        cpgDict[key].append(int(tmp[meth_col]))
        row_count += 1
        if int(tmp[meth_col]) == 1:
            meth_cnt += 1
        else:
            unmeth_cnt += 1

    infile.close()

    logger.info(f"###\timportPredictions_DeepSignal SUCCESS: {row_count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")
    return cpgDict


def importPredictions_DeepSignal3(infileName, chr_col=0, start_col=1, meth_col=7, baseCount=0):
    '''

    Score import for AUC

    Note that the function requires per read stats, not frequencies of methylation.
    !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!

    ### Parameters of the function:
    chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
    start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
    meth_col - name (as header) of the column with methylation call (integer expected). If "header" variable == False, give integer number of the column. *7 is probability, 8 is binary call.
    baseCount - 0 or 1, standing for 0-based or 1-based, respectively
    #header - True or False.

    ### Example input format from DeepSignal:
    chr11   127715633       -       7370988 791c33e4-5b63-4cce-989c-186aff79db9b    t       0.16227172      0.83772826      1       AGGAAAATCGCTTGAAC
    chr11   127715585       -       7371036 791c33e4-5b63-4cce-989c-186aff79db9b    t       0.69724774      0.3027523       0       TGCCACTGCGCTCTAGC
    chr11   127715554       -       7371067 791c33e4-5b63-4cce-989c-186aff79db9b    t       0.786389        0.21361093      0       CAGAACTCCGTCTCAAA
    chr11   127715423       -       7371198 791c33e4-5b63-4cce-989c-186aff79db9b    t       0.5939311       0.40606892      0       TTTTGTAGCGTTGTACA

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

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
    output coordinates are 1-based genome coordinates.

    ** by default if this probability will be higher than 0.5, DeppSignal will tell that this is methylated site, lower, unmethylated

    ============

    '''

    infile = open(infileName, "r")
    cpgDict = {}
    count = 0

    for row in infile:
        tmp = row.strip().split("\t")
        if baseCount == 1:
            start = int(tmp[start_col])
        elif baseCount == 0:
            start = int(tmp[start_col]) + 1
        else:
            logger.error("###\timportPredictions_DeepSignal InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
            sys.exit(-1)
        #         key = (tmp[chr_col], start)
        key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)
        if key not in cpgDict:
            cpgDict[key] = []
        #         cpgDict[key].append(float(tmp[meth_col])) ##### uncomment this line to get probabilities instead of final, binary calls
        cpgDict[key].append(float(tmp[meth_col]))
        count += 1

    infile.close()

    logger.info("###\timportPredictions_DeepSignal SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
    return cpgDict


def importPredictions_Tombo(infileName, chr_col=0, start_col=1, strand_col=5, meth_col=4, baseFormat=0, cutoff=[-1.5, 2.5], output_first=False):
    """
    We treate input as 0-based format.

    Return dict of key='chr1\t123\t123\t+', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.

    Note that the function requires per read stats, not frequencies of methylation.

    if strand-info=='-', the start still point to positive sequence CG's C position, should be more +1 to point to G. We clear the bugs in followings for correct CG, and for reverse strand GC, we need add 1 more to start. Because the sequence is still positive strand sequence, so report position should be at G of CG.

    ### Parameters of the function:
    chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
    start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
    meth_col - name (as header) of the column with methylation call (integer expected). If "header" variable == False, give integer number of the column.
    baseCount - 0 or 1, standing for 0-based or 1-based, respectively
    cutoff - sumilarly as in case of Nanopolish, here we have cutoff for the value from the statistical test. From this conversations (https://github.com/nanoporetech/tombo/issues/151), I know this value is by default 2.5.

    ### Example input format from Tombo (v1.5), we need to pre-process the Tombo results by filtering out non-CG patterns firstly

    more /projects/li-lab/yang/results/2020-12-21/K562.tombo_perReadsStats-with-seq-info-n350-t001-chr1.tsv

    chr1    48020    48020    3526811b-6958-49f8-b78c-a205c1b5fc6e    1.185219591257949    +    TATTACACCCG
    chr1    48022    48022    3526811b-6958-49f8-b78c-a205c1b5fc6e    1.6267354150537658    +    TTACACCCGTT
    chr1    48023    48023    3526811b-6958-49f8-b78c-a205c1b5fc6e    2.6122662196889728    +    TACACCCGTTA
    chr1    48024    48024    3526811b-6958-49f8-b78c-a205c1b5fc6e    2.771131774766473    +    ACACCCGTTAA
    chr1    48041    48041    3526811b-6958-49f8-b78c-a205c1b5fc6e    6.524775544143312    +    GATTTCTAAAT
    chr1    48048    48048    3526811b-6958-49f8-b78c-a205c1b5fc6e    1.9142728191641216    +    AAATGCATTGA
    chr1    48054    48054    3526811b-6958-49f8-b78c-a205c1b5fc6e    1.8675210090110548    +    ATTGACATTTG

    chr1    8447736    8447736    c9339e26-1898-4483-a312-b78c3fafc6a9    8.073560995614967    -    CTGTGCTGTGT
    chr1    8447745    8447745    c9339e26-1898-4483-a312-b78c3fafc6a9    2.4467964154940858    -    GTTGACCGTGT
    chr1    8447746    8447746    c9339e26-1898-4483-a312-b78c3fafc6a9    1.966921521322515    -    TTGACCGTGTA
    chr1    8447754    8447754    c9339e26-1898-4483-a312-b78c3fafc6a9    5.387457000225035    -    GTATGCAATGG
    chr1    8447761    8447761    c9339e26-1898-4483-a312-b78c3fafc6a9    -0.8580941645036908    -    ATGGACACAGA


    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
    output coordinates are 1-based genome coordinates.

    ** by default if this probability will be higher than 0.5, DeppSignal will tell that this is methylated site, lower, unmethylated

    ============

    """
    infile = open(infileName, "r")
    cpgDict = defaultdict(list)
    row_count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        tmp = row.strip().split("\t")

        if tmp[chr_col] not in humanChrs:
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

        # TODO: check how to treate tombo prediction value

        if methCall < cutoff[0]:
            meth_indicator = 1
            meth_cnt += 1
        elif methCall > cutoff[1]:
            meth_indicator = 0
            unmeth_cnt += 1
        else:
            continue

        # if abs(methCall) > cutoff:
        #     meth_indicator = 1
        #     meth_cnt += 1
        # else:
        #     meth_indicator = 0
        #     unmeth_cnt += 1

        cpgDict[key].append(meth_indicator)
        row_count += 1

    infile.close()

    logger.info(f"###\timportPredictions_Tombo SUCCESS: {row_count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-call={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs with meth-cutoff={cutoff} from {infileName} file")
    return cpgDict


def importPredictions_Tombo_nofilter(infileName, chr_col=0, start_col=1, meth_col=4, baseCount=0, cutoff=2.5):
    '''
    Note that the function requires per read stats, not frequencies of methylation.
    !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!

    ### Parameters of the function:
    chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
    start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
    meth_col - name (as header) of the column with methylation call (integer expected). If "header" variable == False, give integer number of the column.
    baseCount - 0 or 1, standing for 0-based or 1-based, respectively
    cutoff - sumilarly as in case of Nanopolish, here we have cutoff for the value from the statistical test. From this conversations (https://github.com/nanoporetech/tombo/issues/151), I know this value is by default 2.5.

    ### Example input format from Tombo (v1.5):
    chr1    66047   66047   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    6.057825558813564       +
    chr1    66053   66053   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    -0.3359579051241508     +
    chr1    66054   66054   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    0.1202407639936725      +
    chr1    66055   66055   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    2.1077369345267907      +
    chr1    66076   66076   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    0.8979673996582611      +

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
    output coordinates are 1-based genome coordinates.

    ** by default if this probability will be higher than 0.5, DeppSignal will tell that this is methylated site, lower, unmethylated

    ============

    '''

    logger.debug(f"importPredictions_Tombo_nofilter infileName={infileName}")
    infile = open(infileName, "r")
    cpgDict = {}
    count = 0

    for row in infile:
        tmp = row.strip().split("\t")
        if baseCount == 1:
            try:
                start = int(tmp[start_col])
            except:
                logger.error(f" ####Tombo parse error at row={row}")
                continue
        elif baseCount == 0:
            try:
                start = int(tmp[start_col]) + 1
            except:
                logger.error(f" ####Tombo parse error at row={row}")
                continue
        else:
            logger.error("###\timportPredictions_Tombo InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
            sys.exit(-1)
        #         key = (tmp[chr_col], start)
        key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)

        try:
            methCall = float(tmp[meth_col])
        except:
            logger.error(f" ####Tombo parse error at row={row}")
            continue

        switch = 0
        #         if methCall > cutoff:
        if methCall < -cutoff:
            methCall = 1
            switch = 1
        #         elif methCall < -cutoff:
        elif methCall > cutoff:
            methCall = 0
            switch = 1
        else:
            methCall = -1
            switch = 1

        if switch == 1:
            if key not in cpgDict:
                cpgDict[key] = []
            cpgDict[key].append(methCall)
        count += 1

    infile.close()

    logger.info("###\timportPredictions_Tombo SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
    return cpgDict


def importPredictions_Tombo3(infileName, chr_col=0, start_col=1, meth_col=4, baseCount=0, cutoff=2.5):
    '''

    Score import for AUC

    Note that the function requires per read stats, not frequencies of methylation.
    !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!

    ### Parameters of the function:
    chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
    start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
    meth_col - name (as header) of the column with methylation call (integer expected). If "header" variable == False, give integer number of the column.
    baseCount - 0 or 1, standing for 0-based or 1-based, respectively
    cutoff - sumilarly as in case of Nanopolish, here we have cutoff for the value from the statistical test. From this conversations (https://github.com/nanoporetech/tombo/issues/151), I know this value is by default 2.5.

    ### Example input format from Tombo (v1.5):
    chr1    66047   66047   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    6.057825558813564       +
    chr1    66053   66053   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    -0.3359579051241508     +
    chr1    66054   66054   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    0.1202407639936725      +
    chr1    66055   66055   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    2.1077369345267907      +
    chr1    66076   66076   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    0.8979673996582611      +

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
    output coordinates are 1-based genome coordinates.

    ** by default if this probability will be higher than 0.5, DeppSignal will tell that this is methylated site, lower, unmethylated

    ============

    '''

    logger.debug(f"importPredictions_Tombo infileName={infileName}")
    infile = open(infileName, "r")
    cpgDict = {}
    count = 0

    for row in infile:
        tmp = row.strip().split("\t")
        if baseCount == 1:
            try:
                start = int(tmp[start_col])
            except:
                logger.error(f" ####Tombo parse error at row={row}")
                continue
        elif baseCount == 0:
            try:
                start = int(tmp[start_col]) + 1
            except:
                logger.error(f" ####Tombo parse error at row={row}")
                continue
        else:
            logger.error("###\timportPredictions_Tombo InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
            sys.exit(-1)
        #         key = (tmp[chr_col], start)
        key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)

        try:
            methCall = float(tmp[meth_col])
        except:
            logger.error(f" ####Tombo parse error at row={row}")
            continue

        switch = 0

        # TODO check if this is inversed
        #         if methCall > cutoff:
        if methCall < -cutoff:
            methCall = -1 * methCall
            switch = 1
        #         elif methCall < -cutoff:
        elif methCall > cutoff:
            methCall = -1 * methCall
            switch = 1

        if switch == 1:
            if key not in cpgDict:
                cpgDict[key] = []
            cpgDict[key].append(methCall)
        count += 1

    infile.close()

    logger.info("###\timportPredictions_Tombo SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
    return cpgDict


def importPredictions_DeepMod(infileName, chr_col=0, start_col=1, strand_col=5, meth_reads_col=-2, coverage_col=-4, baseFormat=0, sep=' ', output_first=False):
    '''
    We treate input as 0-based format for start col. Due to we pre-processed original DeepMod results by filter out non-CG sites, the input of this funciton is sep=TAB instead!!!

    Return dict of key='chr1  123  123  +', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.

    Note that the function requires genome-level stats.

    DeepMod BED format ref: https://github.com/WGLab/DeepMod/blob/master/docs/Results_explanation.md#2-format-of-output

    The output is in a BED format like below. The first six columns are Chr, Start pos, End pos, Base, Capped coverage, and Strand, and the last three columns are Real coverage, Mehylation percentage and Methylation coverage.

    Example input format from DeepMod (standard), we also preprocess the DeepMod initial results by filetering non-CG sites, which is out of our interest.

    head /projects/li-lab/yang/results/2020-12-21/K562.deepmod_combined-with-seq-info-n100-t006-chr1.tsv

    chr1    75694844    75694845    C    1    +        75694844    75694845    0,0,0    1  100    1    TAAGTCGTTCA
    chr1    75696163    75696164    C    1    +        75696163    75696164    0,0,0    1  100    1    CACTCCGGGAC
    chr1    75696217    75696218    C    1    -        75696217    75696218    0,0,0    1  100    1    CATACGGGATA
    chr1    75696583    75696584    C    1    +        75696583    75696584    0,0,0    1  100    1    TATGTCGACTC

    Description (https://github.com/WGLab/DeepMod/blob/master/docs/Results_explanation.md):
    The output is in a BED format like below. The first six columns are Chr,
    Start pos, End pos, Base, Capped coverage, and Strand, and the last three
    columns are Real coverage, Mehylation percentage and Methylation coverage.

    ### Example input format from DeepMod (clustered - following Step 4 from "Example 3: Detect 5mC on Na12878" section; https://github.com/WGLab/DeepMod/blob/master/docs/Reproducibility.md):
    chr2 241991445 241991446 C 3 -  241991445 241991446 0,0,0 3 100 3 69
    chr2 241991475 241991476 C 3 -  241991475 241991476 0,0,0 3 33 1 75
    chr2 241991481 241991482 C 2 -  241991481 241991482 0,0,0 2 50 1 76

    Note: it is space-separated in original result file, not tab-separated file

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
    ============

    '''

    infile = open(infileName, "r")
    cpgDict = defaultdict(list)
    count_calls = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        tmp = row.strip().split(sep)

        if tmp[chr_col] not in humanChrs:
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
            logger.debug("###\timportPredictions_DeepMod InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')

        key = (tmp[chr_col], start, strand)

        methReads = int(tmp[meth_reads_col])
        coverage = int(tmp[coverage_col])

        methCallsList = [1] * methReads + [0] * (coverage - methReads)
        cpgDict[key] = methCallsList

        count_calls += len(methCallsList)
        meth_cnt += methReads
        unmeth_cnt += coverage - methReads

    infile.close()

    logger.info(f"###\timportPredictions_DeepMod SUCCESS: {count_calls:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")
    return cpgDict


def importPredictions_DeepMod_Read_Level(infileName, chr_col=0, start_col=1, strand_col=5, meth_col=-1, baseFormat=0, sep='\t', output_first=False):
    infile = open(infileName, "r")
    # cpgDict = {}
    cpgDict = defaultdict(list)
    count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        tmp = row.strip().split(sep)

        if tmp[chr_col] not in humanChrs:
            continue

        if output_first:
            logger.debug(f'row = {list(enumerate(tmp))}')
            output_first = False
        if baseFormat == 1:
            start = int(tmp[start_col]) + 1
            end = start
            strand = tmp[strand_col]
        elif baseFormat == 0:
            start = int(tmp[start_col])
            end = start + 1
            strand = tmp[strand_col]
        else:
            logger.debug("###\timportPredictions_DeepMod_Read_Level InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseFormat))
            sys.exit(-1)

        if strand not in ['-', '+']:
            raise Exception(f'The file [{infileName}] can not recognized strand-info from row={row}, please check it')
        #         key = (tmp[chr_col], start)
        key = (tmp[chr_col], start, strand)

        meth_indicator = int(tmp[meth_col])

        if meth_indicator == 0:
            unmeth_cnt += 1
        else:
            meth_cnt += 1

        # methCallsList = [1] * methReads + [0] * (coverage - methReads)
        cpgDict[key].append(meth_indicator)

        count += 1

    infile.close()

    logger.info(f"###\timportPredictions_DeepMod_Read_Level SUCCESS: {count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-calls={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs from {infileName} file")
    return cpgDict


def importPredictions_DeepMod3(infileName, chr_col=0, start_col=1, meth_percentage_col=11, coverage_col=10, clusteredResult=False, clustered_meth_freq_col=13, baseCount=0):
    '''

    each CpGs sites a methylation percentage
    Note that the function requires per read stats, not frequencies of methylation.
    !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!

    ### Parameters of the function:
    chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
    start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
    meth_reads_col - name (as header) of the column with number of methylated reads mapped.
    coverage_col - name (as header) of the column with coverage of the site.
    [[TO DO]] clusteredResult - True / False. Input file is in the "clustered" format (additional post-processing step). False (default option) - standard output with calls.
    [[TO DO]] clustered_meth_freq_col - column with the methylation frequency after additional postprocessing step.
    baseCount - 0 or 1, standing for 0-based or 1-based, respectively

    ### Example input format from DeepMod (standard):
    chr2 110795922 110795923 C 4 -  110795922 110795923 0,0,0 4 75 3
    chr2 110795929 110795930 C 3 -  110795929 110795930 0,0,0 3 66 2
    chr2 110796453 110796454 C 4 -  110796453 110796454 0,0,0 4 25 1

    Description (https://github.com/WGLab/DeepMod/blob/master/docs/Results_explanation.md):
    The output is in a BED format like below. The first six columns are Chr,
    Start pos, End pos, Base, Capped coverage, and Strand, and the last three
    columns are Real coverage, Mehylation percentage and Methylation coverage.


    ### Example input format from DeepMod (clustered - following Step 4 from "Example 3: Detect 5mC on Na12878" section; https://github.com/WGLab/DeepMod/blob/master/docs/Reproducibility.md):
    chr2 241991445 241991446 C 3 -  241991445 241991446 0,0,0 3 100 3 69
    chr2 241991475 241991476 C 3 -  241991475 241991476 0,0,0 3 33 1 75
    chr2 241991481 241991482 C 2 -  241991481 241991482 0,0,0 2 50 1 76

    Note: it is space-separated, not tab-separated file

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
    output coordinates are 1-based genome coordinates.

    ============

    '''

    infile = open(infileName, "r")
    cpgDict = {}
    count = 0

    for row in infile:
        tmp = row.strip().split(" ")
        if baseCount == 1:
            start = int(tmp[start_col])
        elif baseCount == 0:
            start = int(tmp[start_col]) + 1
        else:
            logger.error("###\timportPredictions_DeepMod InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
            sys.exit(-1)
        #         key = (tmp[chr_col], start)
        key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)

        methval = int(tmp[meth_percentage_col])
        coverage = int(tmp[coverage_col])

        # methCallsList = [1] * methCalls + [0] * (coverage - methCalls)
        if key not in cpgDict:
            cpgDict[key] = methval / 100.0
        else:
            raise Exception(f"Duplicate CpG sites in DeepMod, key={key}")
        count += 1

    infile.close()

    logger.info("###\timportPredictions_DeepMod SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
    return cpgDict


def importPredictions_DeepMod_clustered(infileName, chr_col=0, start_col=1, strand_col=5, coverage_col=-4, clustered_meth_freq_col=-1, baseFormat=0, sep=' '):
    '''
    Note: results of cluster is differ from other tools like:
    [methFrequency, coverage] such as key -> values [50 (freq 0-100), 10 (cov)]

    Note that the function requires per read stats, not frequencies of methylation.
    !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!

    ### Parameters of the function:
    chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
    start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
    meth_reads_col - name (as header) of the column with number of methylated reads mapped.
    coverage_col - name (as header) of the column with coverage of the site.
    [[TO DO]] clusteredResult - True / False. Input file is in the "clustered" format (additional post-processing step). False (default option) - standard output with calls.
    [[TO DO]] clustered_meth_freq_col - column with the methylation frequency after additional postprocessing step.
    baseCount - 0 or 1, standing for 0-based or 1-based, respectively

    ### Example input format from DeepMod (standard):
    chr2 110795922 110795923 C 4 -  110795922 110795923 0,0,0 4 75 3
    chr2 110795929 110795930 C 3 -  110795929 110795930 0,0,0 3 66 2
    chr2 110796453 110796454 C 4 -  110796453 110796454 0,0,0 4 25 1

    Description (https://github.com/WGLab/DeepMod/blob/master/docs/Results_explanation.md):
    The output is in a BED format like below. The first six columns are Chr,
    Start pos, End pos, Base, Capped coverage, and Strand, and the last three
    columns are Real coverage, Mehylation percentage and Methylation coverage.


    ### Example input format from DeepMod (clustered - following Step 4 from "Example 3: Detect 5mC on Na12878" section; https://github.com/WGLab/DeepMod/blob/master/docs/Reproducibility.md):
    chr2 241991445 241991446 C 3 -  241991445 241991446 0,0,0 3 100 3 69
    chr2 241991475 241991476 C 3 -  241991475 241991476 0,0,0 3 33 1 75
    chr2 241991481 241991482 C 2 -  241991481 241991482 0,0,0 2 50 1 76

    Note: it is space-separated, not tab-separated file

    ### Output format:
    result = {"chr\tstart\tend\n" : [methFrequency, coverage]}
    output coordinates are 1-based genome coordinates.

    ============

    '''

    infile = open(infileName, "r")
    cpgDict = defaultdict(dict)
    count = 0
    output_first = True

    for row in infile:
        tmp = row.strip().split(sep)

        if tmp[chr_col] not in humanChrs:
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

        methFreq = int(tmp[clustered_meth_freq_col])
        coverage = int(tmp[coverage_col])

        cpgDict[key] = {'freq': methFreq / 100.0, 'cov': coverage}

        count += 1

    infile.close()

    logger.info(f"###\tDeepMod_clusteredResult Parsing SUCCESS: imported for {len(cpgDict):,} CpGs from {infileName} file")

    return cpgDict


def importPredictions_Megalodon_Read_Level(infileName, chr_col=0, start_col=1, strand_col=3, mod_log_prob_col=4, can_log_prob_col=5, baseFormat=0, cutoff=0.8, sep='\t', output_first=False):
    '''

    0-based start for Magelodon：
        1. baseFormat=0， start=Megalondon start；
        2.  baseFormat=1， start=Megalondon start +1

    We assume the input is 0-based for the start col, such as 'chr10\t122'
    Return dict of key='chr1\t123\t124\t+', and values=list of [1 1 0 0 1 1], in which 0-unmehylated, 1-methylated.

    Note that the function requires per-read file, not per-site methlyation level.

    ### Parameters of the function:
    chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
    start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
    mod_log_prob_col - name (as header) of the column with natural log probability that the base is modified (float expected).
    can_log_prob_col - name (as header) of the column with natural log probability that the base is unmodified (float expected).
    baseCount - 0 or 1, standing for 0-based or 1-based, respectively.

    ### Example default input format from Megalondon (pre-processed to containing strand-info):
    head /fastscratch/c-panz/K562/megalodon/per_read_modified_base_calls.merged.sorted.bed
    chr1	10468	a0b5e934-3784-4163-b46c-f575ac1015bf	+	-1.6103846163370363	-0.2229070153110119	    m
    chr1	10470	a0b5e934-3784-4163-b46c-f575ac1015bf	+	-2.035311286540776	-0.13999775139928522	m
    chr1	10483	a0b5e934-3784-4163-b46c-f575ac1015bf	+	-1.5477196338381982	-0.2391872270542014 	m

    ### Output format:
    result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
    output coordinates are 1-based genome coordinates.

    ** By default setting in Megalodon, for a modified base, mod_prob > 0.8; for a canonical base, can_prob > 0.8 (In other words, mod_prob < 0.2)
    For any position where no probability is greater than 0.8 neither canonical not modified bases get a count.
    (https://github.com/nanoporetech/megalodon/issues/47#issuecomment-673742805)

    ============
    '''
    infile = open(infileName, "r")
    cpgDict = defaultdict(list)
    row_count = 0
    meth_cnt = 0
    unmeth_cnt = 0

    for row in infile:
        tmp = row.strip().split(sep)

        if tmp[chr_col] not in humanChrs:
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
            methCall = float(np.e ** float(tmp[mod_log_prob_col]))  # Calculate mod_prob
        except Exception as e:
            logger.error(f" ####Megalodon parse error at row={row}, exception={e}")
            continue

        if methCall > cutoff:  ##Keep methylated reads
            meth_indicator = 1
            meth_cnt += 1
        elif methCall < 1 - cutoff:  ##Count unmethylated reads
            meth_indicator = 0
            unmeth_cnt += 1
        else:  ## Neglect other cases 0.2<= prob <=0.8
            continue

        cpgDict[key].append(meth_indicator)
        row_count += 1

    infile.close()

    logger.info(f"###\timportPredictions_Megalodon SUCCESS: {row_count:,} methylation calls (meth-calls={meth_cnt:,}, unmeth-call={unmeth_cnt:,}) mapped to {len(cpgDict):,} CpGs with meth-cutoff={cutoff:.2f} from {infileName} file")
    return cpgDict


def coverageFiltering(calls_dict, minCov=4, byLength=True, toolname="Tool"):
    """
    Convert orignial call object from dict[cpg] = {cpg: [meth_freq_read1, ..., meth_freq_readn]} to dict[cpg] = {cpg:[meth_freq, coverage_number]}

    meth_freq   in [0.0,1.0]
    cov_num     in int
    Read-level -> Genome-level

    :param calls_dict:
    :param minCov:
    :param byLength: if False, will deal with DeepMod_cluster results
    :return:
    """
    result = defaultdict()
    for cpg in calls_dict:
        if type(calls_dict[cpg]) == list:  # value is [0 0 1 1 0 ...]
            if len(calls_dict[cpg]) >= minCov:
                result[cpg] = [sum(calls_dict[cpg]) / float(len(calls_dict[cpg])), len(calls_dict[cpg])]
        else:  # Used by DeepMod_cluster results, value is {'freq':0.72, 'cov':18}
            if calls_dict[cpg]['cov'] >= minCov:
                result[cpg] = [calls_dict[cpg]['freq'], calls_dict[cpg]['cov']]

    logger.info(f"###\tcoverageFiltering: completed filtering with minCov={minCov}, {len(result):,} CpG sites left for {toolname}")
    return result


def importGroundTruth_oxBS(infileName, chr_col='#chromosome', start_col='start', meth_col="pmC", covCutt=4, baseCount=1, chrFilter=False):
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


def importGroundTruth_BedMethyl_from_Encode(infileName, chr_col=0, start_col=1, meth_col=10, cov_col=9, covCutt=10, baseCount=1, chrFilter=False):
    '''

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

    I think that this output is 0-based.

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

    ### Output files:
    cpgDict = {"chr\tstart\tend\n" : methylation level (format: float (0-1))} # have all cpgs with specified cutoff (not only 100% 5mC or 100% 5C)
    output coordinates are 1-based genome coordinates.

    '''

    cpgDict = {}

    infile = open(infileName, 'r')

    nrow = 0
    for row in infile:
        nrow += 1
        tmp = row.strip().split("\t")
        if baseCount == 1:
            start = int(tmp[start_col]) + 1
            end = start
        elif baseCount == 0:
            start = int(tmp[start_col])
            end = start + 1
        else:
            logger.error("###\timportGroundTruth_BedMethyl_from_Encode InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
            sys.exit(-1)

        if chrFilter == False or chrFilter == tmp[chr_col]:
            if int(tmp[cov_col]) >= covCutt:
                key = "{}\t{}\t{}\n".format(tmp[chr_col], start, end)
                if key not in cpgDict:
                    cpgDict[key] = float(tmp[meth_col]) / 100.0
                else:
                    logger.error("###\timportGroundTruth_BedMethyl_from_Encode SanityCheckError: One CpG should not have more than 1 entry")

    infile.close()
    logger.debug("###\timportGroundTruth_BedMethyl_from_Encode: loaded information for {} CpGs, ({} rows)".format(len(cpgDict), nrow))
    return cpgDict


# Deal with K562 bed with strand info
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

        if tmp[chr_col] not in humanChrs:
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
            sys.exit(-1)

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


def plot_AUC_curve(scores, y, ax, title="", outfile=None):
    try:
        fpr, tpr, _ = roc_curve(y, scores)
        roc_auc = auc(fpr, tpr)
        lw = 2

        plt.plot(fpr, tpr,
                 lw=lw, label='{0} - ROC curve (area = {1:.4f})'.format(title, roc_auc))
        plt.plot([0, 1], [0, 1], color='lightgrey', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title)
        plt.legend(loc="lower right")
    except ValueError:
        logger.error(f"###\tERROR for plot_AUC_curve: y: {y}, scores: {scores}")


def importGroundTruth_BS():
    '''
    !!! This function will be needed for NanoCompare project !!!
    !!!!!!! its not a final version of the function, but with only some small changes it will be

    note that this function was written to parse "bedMethyl" format from Encode. (description here: https://www.encodeproject.org/documents/964e2676-d0be-4b5d-aeec-f4f02310b221/@@download/attachment/WGBS%20pipeline%20overview.pdf)
    additionally, the output file was preprocessed, such that only sites with coverage >= 10 reads and either 100% or 0% methylated, are included.
    these will be taken to compute true positive, true negative etc.

    #####
    example usage here: http://helix122:9912/edit/li-lab/NanoporeData/WR_ONT_analyses/ai/APL_nanopolishStats/automatedSingleReadPrecission_3.py

    '''

    #     infile = open(infileName, 'r')

    #     bsseqDict = {} # {"chr:start:end" : 100 or 0} - 100 in case if methylated, 0 if unmethylated
    #     for row in infile:
    #         tmp = row.strip().split("\t")
    #         bsseqDict["{}:{}:{}".format(tmp[0], tmp[1], tmp[2])] = int(tmp[-1])
    #     infile.close()
    #     return bsseqDict
    pass


def importNarrowCpGsList():
    pass


def dict2txt(inputDict):
    """
    Convert all keys in dict to a string txt, txt file is:
    chr1  123  123  +

    :param inputDict:
    :return:
    """
    text = ""
    for key in inputDict:
        #         print(key)
        # text += key
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


# Deprecated
def computePerReadStats_deprecated(ontCalls, bsReference, title, bedFile=False, ontCutt_perRead=1, ontCutt_4corr=4, secondFilterBed=False, secondFilterBed_4Corr=False):
    '''
    ontCalls - dictionary of CpG coordinates with their per read methylation call (1 or 0) // Format: {"chr\tstart\tend\n" : [list of methylation calls]}
    bsReference - dictionary of CpG coordinates with their methylation frequencies (range 0 - 1). This list is already prefiltered to meet minimal coverage (e.g. 4x) at this point. // Format: {"chr\tstart\tend\n" : methylation level (format: float (0-1))}
    title - prefix of the analysis, output plots etc. - should be as short as possible, but unique in context of other analyses
    bedFile - BED file which will be used to narrow down the list of CpGs for example to those inside CGIs or promoters etc.. By default "False" - which means no restrictions are done (i.e. genome wide)
    secondFilterBed - these should be CpGs covered in some reference list. Format: BED


    ============================================

    Basically i want to fill in the table below:
               Positive	Negative	Total
     Presence	a	        b	        a+b
     Absence	c	        d	        c+d
     Total	    a+c	        b+d	        a+b+c+d

    Nice summary also at wiki: https://en.wikipedia.org/wiki/F1_score
    , where "Positive" and "Negative" corresponds with ONT based observations, while "Presence" and "Absence" is based on BS-Seq

    '''

    switch = 0
    ontCalls_narrow = []

    coord_base_dir = os.path.join(data_base_dir, 'genome-annotation')

    if bedFile != False:
        ontCalls_bed = BedTool(dict2txt(ontCalls), from_string=True)
        ontCalls_bed = ontCalls_bed.sort()

        infn = os.path.join(coord_base_dir, bedFile)
        # narrowBed = BedTool(bedFile)
        narrowBed = BedTool(infn)
        narrowBed = narrowBed.sort()
        ontCalls_intersect = ontCalls_bed.intersect(narrowBed, u=True, wa=True)
        ontCalls_narrow = txt2dict(ontCalls_intersect)
        suffix = bedFile
    else:
        switch = 1
        suffix = "GenomeWide"

    ## Second optional filter, organized in the same fashion as the first one. Designed to accomodate for example list of CpGs covered by second program
    secondSwitch = 0
    ontCalls_narrow_second = {}
    if secondFilterBed != False:
        infile = open(secondFilterBed, 'r')
        secondFilterDict = {}
        for row in infile:
            secondFilterDict[row] = 0
        infile.close()
        ontCalls_narrow_second = dict.fromkeys(set(ontCalls.keys()).intersection(set(secondFilterDict.keys())), 0)
    else:
        secondSwitch = 1
        suffix = "GenomeWide"

    # Second optional filter, shoudl be used in combination with second optional filter above
    secondSwitch_4corr = 0
    ontCalls_narrow_second_4corr = {}
    if secondFilterBed_4Corr != False:
        infile = open(secondFilterBed_4Corr, 'r')
        secondFilterDict = {}
        for row in infile:
            secondFilterDict[row] = 0
        infile.close()
        ontCalls_narrow_second_4corr = dict.fromkeys(set(ontCalls.keys()).intersection(set(secondFilterDict.keys())), 0)
    else:
        secondSwitch_4corr = 1
        suffix = "GenomeWide"

    suffixTMP = suffix.split("/")
    if len(suffixTMP) > 1:
        suffix = suffixTMP[-1]

    TP_5mC = FP_5mC = FN_5mC = TN_5mC = TP_5C = FP_5C = FN_5C = TN_5C = 0
    y = []
    scores = []

    ontSites = 0
    mCsites = 0
    Csites = 0
    referenceCpGs = 0

    ## four tuples for correlation:
    ontFrequencies_4corr_mix = []  # mix(ed) are those, which in reference have methylation level >0 and <1
    refFrequencies_4corr_mix = []

    ontFrequencies_4corr_all = []  # all are all:) i.e. all CpGs with methyaltion level in refence in range 0-1
    refFrequencies_4corr_all = []
    leftovers = {}
    leftovers1 = {}

    for cpg_ont in ontCalls:
        ##### for per read stats:
        if cpg_ont in bsReference and len(ontCalls[cpg_ont]) >= ontCutt_perRead and (switch == 1 or cpg_ont in ontCalls_narrow) and (
                secondSwitch == 1 or cpg_ont in ontCalls_narrow_second):  # we should not take onCuttoffs for per read stats - shouldn't we? actually, we need to have the option to use this parameter, because at some point we may want to narrow down the per read stats to cover only the sites which were also covered by correlation with BS-Seq. Using this cutoff here is the easiest way to do just that
            #         if cpg_ont in bsReference and (switch == 1 or cpg_ont in ontCalls_narrow):
            if bsReference[cpg_ont] == 1 or bsReference[cpg_ont] == 0:  # we only consider absolute states here
                referenceCpGs += 1

                for ontCall in ontCalls[cpg_ont]:  # TODO  mCsites and Csites will be count more due to more reads
                    if bsReference[cpg_ont] == 1:
                        mCsites += 1
                    if bsReference[cpg_ont] == 0:
                        Csites += 1
                    ontSites += 1

                    ### variables needed to compute precission, recall etc.:
                    if ontCall == 1 and bsReference[cpg_ont] == 1:  # true positive
                        TP_5mC += 1
                    elif ontCall == 1 and bsReference[cpg_ont] == 0:  # false positive
                        FP_5mC += 1
                    elif ontCall == 0 and bsReference[cpg_ont] == 1:  # false negative
                        FN_5mC += 1
                    elif ontCall == 0 and bsReference[cpg_ont] == 0:  # true negative
                        TN_5mC += 1

                    if ontCall == 0 and bsReference[cpg_ont] == 0:  # true positive
                        TP_5C += 1
                    elif ontCall == 0 and bsReference[cpg_ont] == 1:  # false positive
                        FP_5C += 1
                    elif ontCall == 1 and bsReference[cpg_ont] == 0:  # false negative
                        FN_5C += 1
                    elif ontCall == 1 and bsReference[cpg_ont] == 1:  # true negative
                        TN_5C += 1

                    ### AUC related:
                    scores.append(ontCall)
                    y.append(bsReference[cpg_ont])

        ##### for correlation stats:
        if cpg_ont in bsReference and len(ontCalls[cpg_ont]) >= ontCutt_4corr and (switch == 1 or cpg_ont in ontCalls_narrow) and (secondSwitch_4corr == 1 or cpg_ont in ontCalls_narrow_second_4corr):
            ontMethFreq = np.mean(ontCalls[cpg_ont])
            ontFrequencies_4corr_all.append(ontMethFreq)
            refFrequencies_4corr_all.append(bsReference[cpg_ont])
            if bsReference[cpg_ont] > 0 and bsReference[cpg_ont] < 1:
                ontFrequencies_4corr_mix.append(ontMethFreq)
                refFrequencies_4corr_mix.append(bsReference[cpg_ont])

    ### compute all per read stats:

    #     Accuracy:
    try:
        accuracy = (TP_5mC + TN_5mC) / float(TP_5mC + FP_5mC + FN_5mC + TN_5mC)
    except ZeroDivisionError:
        accuracy = 0
    #     print("Accuracy: {0:1.4f}".format(accuracy))

    #     Positive predictive value (PPV), Precision = (TP) / E(Predicted condition positive)
    try:
        predicted_condition_positive_5mC = float(TP_5mC + FP_5mC)
        precision_5mC = TP_5mC / predicted_condition_positive_5mC
    except ZeroDivisionError:
        precision_5mC = 0
    #     print("Precision_5mC: {0:1.4f}".format(precision_5mC))
    try:
        predicted_condition_positive_5C = float(TP_5C + FP_5C)
        precision_5C = TP_5C / predicted_condition_positive_5C
    except ZeroDivisionError:
        precision_5C = 0

    #     print("Precision_5C: {0:1.4f}".format(precision_5C))

    #     True positive rate (TPR), Recall, Sensitivity, probability of detection = (TP) / (TP+FN)
    try:
        recall_5mC = TP_5mC / float(TP_5mC + FN_5mC)
    except ZeroDivisionError:
        recall_5mC = 0
    #     print("Recall_5mC: {0:1.4f}".format(recall_5mC))

    try:
        recall_5C = TP_5C / float(TP_5C + FN_5C)
    except ZeroDivisionError:
        recall_5C = 0
    #     print("Recall_5C: {0:1.4f}".format(recall_5C))

    #     F1 score:
    try:
        F1_5mC = 2 * ((precision_5mC * recall_5mC) / (precision_5mC + recall_5mC))
    except ZeroDivisionError:
        F1_5mC = 0
    #     print("F1 score_5mC: {0:1.4f}".format(F1_5mC))

    try:
        F1_5C = 2 * ((precision_5C * recall_5C) / (precision_5C + recall_5C))
    except ZeroDivisionError:
        F1_5C = 0
    #     print("F1 score_5C: {0:1.4f}".format(F1_5C))

    #     print("ontSites:", ontSites)

    ## plot AUC curve:
    fig = plt.figure(figsize=(5, 5), dpi=300)

    fprSwitch = 1
    try:
        fpr, tpr, _ = roc_curve(y, scores)
    except ValueError:
        logger.error(f"###\tERROR for roc_curve: y:{y}, scores: {scores}, \nother settings: {title}, {bedFile}, {secondFilterBed}, {secondFilterBed_4Corr}")
        fprSwitch = 0
        roc_auc = 0

    if fprSwitch == 1:
        roc_auc = auc(fpr, tpr)
        #     print("AUC: {0:1.4f}".format(roc_auc))
        #     print(title)

        lw = 2

        plt.plot(fpr, tpr, lw=lw, label='ROC curve (area = {0:.4f})'.format(roc_auc))
        plt.plot([0, 1], [0, 1], color='lightgrey', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(suffix)
        plt.legend(loc="lower right")

        fig.savefig("{}.{}.AUC.pdf".format(title, suffix), bbox_inches='tight', dpi=300)
        plt.close()

        ## Plot confusion matrix:
        cnf_matrix = confusion_matrix(y, scores)
        np.set_printoptions(precision=2)
        plt.figure(figsize=(5, 5))
        plot_confusion_matrix(cnf_matrix, classes=[0, 1], normalize=True, title=suffix)
        plt.savefig("{}.{}.ConfusionMatrix.pdf".format(title, suffix), bbox_inches='tight', dpi=300)
        plt.close()

        ## plot Precission-recall:
        average_precision = average_precision_score(y, scores)
        #     print('Average precision-recall score: {0:0.2f}'.format(average_precision))
        plt.figure(figsize=(5, 5))
        precision, recall, _ = precision_recall_curve(y, scores)
        # In matplotlib < 1.5, plt.fill_between does not have a 'step' argument
        step_kwargs = ({'step': 'post'}
                       if 'step' in signature(plt.fill_between).parameters
                       else {})
        plt.step(recall, precision, color='b', alpha=0.2,
                 where='post')
        plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)

        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('2-class Precision-Recall curve: AP={0:0.2f}\n{1}'.format(average_precision, suffix))
        plt.savefig("{}.{}.PrecissionRecall.pdf".format(title, suffix), bbox_inches='tight', dpi=300)
        plt.close()

    ########################
    # correlation based stats:
    try:
        corrMix, pvalMix = pearsonr(ontFrequencies_4corr_mix, refFrequencies_4corr_mix)
    except:
        corrMix, pvalMix = (0, 0)

    try:
        corrAll, pvalAll = pearsonr(ontFrequencies_4corr_all, refFrequencies_4corr_all)
    except:
        corrAll, pvalAll = (0, 0)

    return accuracy, roc_auc, precision_5C, recall_5C, F1_5C, Csites, precision_5mC, recall_5mC, F1_5mC, mCsites, referenceCpGs, corrMix, len(ontFrequencies_4corr_mix), corrAll, len(ontFrequencies_4corr_all)  # , leftovers, leftovers1


def load_single_sites_bed_as_set(infn):
    infile = open(infn, 'r')
    ret = set()
    for row in infile:  # each row: chr123  123   123  .  .  +
        rowsplit = row[:-1].split('\t')
        key = (rowsplit[0], int(rowsplit[1]), rowsplit[5])
        ret.add(key)
    infile.close()
    return ret


def computePerReadStats(ontCalls, bgTruth, title, coordBedFileName=None, ontCutt_perRead=1, ontCutt_4corr=3, secondFilterBedFile=None, secondFilterBed_4CorrFile=None, cutoff_meth=1.0, outdir=None, tagname=None):
    """
    Compute ontCalls with bgTruth performance results by per-read count.
    bedFile                 -   coordinate used to eval
    secondFilterBed         -   joined sets of four tools with bg-truth, or False with out joined sets
    secondFilterBed_4Corr   -   for report corr purpose, currently not fixed bugs.

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


    ============================================

    Basically i want to fill in the table below:
               Positive	Negative	Total
     Presence	a	        b	        a+b
     Absence	c	        d	        c+d
     Total	    a+c	        b+d	        a+b+c+d

    Nice summary also at wiki: https://en.wikipedia.org/wiki/F1_score
    , where "Positive" and "Negative" corresponds with ONT based observations, while "Presence" and "Absence" is based on BS-Seq

    """

    switch = 0
    ontCalls_narrow_set = None
    if coordBedFileName:
        # logger.debug(bedFile)
        ontCalls_bed = BedTool(dict2txt(ontCalls), from_string=True)
        ontCalls_bed = ontCalls_bed.sort()

        narrowBed = BedTool(coordBedFileName)
        narrowBed = narrowBed.sort()
        ontCalls_intersect = ontCalls_bed.intersect(narrowBed, u=True, wa=True)
        ontCalls_narrow_set = set(txt2dict(ontCalls_intersect).keys())
    else:
        switch = 1

    ## Second optional filter, organized in the same fashion as the first one. Designed to accomodate for example list of CpGs covered by second program
    secondSwitch = 0
    ontCalls_narrow_second_set = None
    if secondFilterBedFile:
        joined_set = load_single_sites_bed_as_set(secondFilterBedFile)
        infile = open(secondFilterBedFile, 'r')
        secondFilterDict = {}
        for row in infile:  # each row: chr123  123   123  .  .  +
            rowsplit = row[:-1].split('\t')
            key = (rowsplit[0], int(rowsplit[1]), rowsplit[5])
            secondFilterDict[key] = 0
        infile.close()
        ontCalls_narrow_second_set = set(ontCalls.keys()).intersection(joined_set)
    else:
        secondSwitch = 1

    # Second optional filter, shoudl be used in combination with second optional filter above
    secondSwitch_4corr = 0
    ontCalls_narrow_second_4corr_set = None
    if secondFilterBed_4CorrFile != False:
        narrow_second_4corr_set = load_single_sites_bed_as_set(secondFilterBed_4CorrFile)
        infile = open(secondFilterBed_4CorrFile, 'r')
        secondFilterDict = {}
        for row in infile:
            rowsplit = row[:-1].split('\t')
            key = (rowsplit[0], int(rowsplit[1]), rowsplit[5])
            secondFilterDict[key] = 0
        infile.close()
        ontCalls_narrow_second_4corr_set = set(ontCalls.keys()).intersection(narrow_second_4corr_set)
    else:
        secondSwitch_4corr = 1

    # Initial evaluation vars
    TP_5mC = FP_5mC = FN_5mC = TN_5mC = TP_5C = FP_5C = FN_5C = TN_5C = 0
    y_of_bgtruth = []
    ypred_of_ont_tool = []

    ontSites = 0  # count how many reads is methy or unmethy
    mCalls = 0  # count how many read call is methy
    cCalls = 0  # count how many read call is unmethy

    referenceCpGs = 0

    # really CpGs for fully methylation and unmethylation,
    mCsites_BGTruth = 0
    Csites_BGTruth = 0

    ## four tuples for correlation:
    ontFrequencies_4corr_mix = []  # mix(ed) are those, which in reference have methylation level >0 and <1
    refFrequencies_4corr_mix = []

    ontFrequencies_4corr_all = []  # all are all:) i.e. all CpGs with methyaltion level in refence in range 0-1
    refFrequencies_4corr_all = []

    # We find the narrowed CpG set to evaluate, try to reduce running time
    targetedSet = set(ontCalls.keys())
    if ontCalls_narrow_set:
        targetedSet = targetedSet.intersection(ontCalls_narrow_set)

    if ontCalls_narrow_second_set:
        targetedSet = targetedSet.intersection(ontCalls_narrow_second_set)

    if ontCalls_narrow_second_4corr_set:
        targetedSet = targetedSet.union(ontCalls_narrow_second_4corr_set)

    for cpgKey in targetedSet:  # key = (chr, start, strand)
        ##### for per read stats:
        if cpgKey in bgTruth and len(ontCalls[cpgKey]) >= ontCutt_perRead and (switch == 1 or cpgKey in ontCalls_narrow_set) and (
                secondSwitch == 1 or cpgKey in ontCalls_narrow_second_set):  # we should not take onCuttoffs for per read stats - shouldn't we? actually, we need to have the option to use this parameter, because at some point we may want to narrow down the per read stats to cover only the sites which were also covered by correlation with BS-Seq. Using this
            # cutoff here is the easiest way to do just that
            #         if cpg_ont in bsReference and (switch == 1 or cpg_ont in ontCalls_narrow):
            if bgTruth[cpgKey][0] >= (cutoff_meth - 1e-6) or bgTruth[cpgKey][0] <= 1e-5:  # we only consider absolute states here, 0, or 1 in bgtruth
                referenceCpGs += 1

                if bgTruth[cpgKey][0] >= (cutoff_meth - 1e-6):
                    mCsites_BGTruth += 1
                elif bgTruth[cpgKey][0] <= 1e-5:
                    Csites_BGTruth += 1

                for perCall in ontCalls[cpgKey]:  # perCall only 0 or 1
                    if perCall >= (cutoff_meth - 1e-6):
                        mCalls += 1
                    elif perCall <= 1e-5:
                        cCalls += 1
                    ontSites += 1

                    ### variables needed to compute precission, recall etc.:
                    if perCall == 1 and bgTruth[cpgKey][0] >= (cutoff_meth - 1e-6):  # true positive
                        TP_5mC += 1
                    elif perCall == 1 and bgTruth[cpgKey][0] <= 1e-5:  # false positive
                        FP_5mC += 1
                    elif perCall == 0 and bgTruth[cpgKey][0] >= (cutoff_meth - 1e-6):  # false negative
                        FN_5mC += 1
                    elif perCall == 0 and bgTruth[cpgKey][0] <= 1e-5:  # true negative
                        TN_5mC += 1

                    if perCall == 0 and bgTruth[cpgKey][0] <= 1e-5:  # true positive
                        TP_5C += 1
                    elif perCall == 0 and bgTruth[cpgKey][0] >= (cutoff_meth - 1e-6):  # false positive
                        FP_5C += 1
                    elif perCall == 1 and bgTruth[cpgKey][0] <= 1e-5:  # false negative
                        FN_5C += 1
                    elif perCall == 1 and bgTruth[cpgKey][0] >= (cutoff_meth - 1e-6):  # true negative
                        TN_5C += 1

                    ### AUC related:
                    ypred_of_ont_tool.append(perCall)

                    if bgTruth[cpgKey][0] >= (cutoff_meth - 1e-6):
                        y_of_bgtruth.append(1)
                    else:
                        y_of_bgtruth.append(0)

        ##### for correlation stats: # TODO Currently not used
        if cpgKey in bgTruth and len(ontCalls[cpgKey]) >= ontCutt_4corr and (switch == 1 or cpgKey in ontCalls_narrow_set) and (secondSwitch_4corr == 1 or cpgKey in ontCalls_narrow_second_4corr_set):
            ontMethFreq = np.mean(ontCalls[cpgKey])
            ontFrequencies_4corr_all.append(ontMethFreq)
            refFrequencies_4corr_all.append(bgTruth[cpgKey][0])
            if bgTruth[cpgKey][0] > 0 and bgTruth[cpgKey][0] < (cutoff_meth - 1e-6):
                ontFrequencies_4corr_mix.append(ontMethFreq)
                refFrequencies_4corr_mix.append(bgTruth[cpgKey][0])

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
        fpr, tpr, _ = roc_curve(y_of_bgtruth, ypred_of_ont_tool)
    except ValueError:
        logger.error(f"###\tERROR for roc_curve: y(Truth):{y_of_bgtruth}, scores(Call pred):{ypred_of_ont_tool}, \nother settings: {title}, {coordBedFileName}, {secondFilterBedFile}, {secondFilterBed_4CorrFile}")
        fprSwitch = 0
        roc_auc = 0

    if fprSwitch == 1:
        roc_auc = auc(fpr, tpr)

        ## Precission-recall:
        average_precision = average_precision_score(y_of_bgtruth, ypred_of_ont_tool)
        precision, recall, _ = precision_recall_curve(y_of_bgtruth, ypred_of_ont_tool)

    ########################
    # save y and y-pred for later plot:
    curve_data = {'yTrue': y_of_bgtruth, 'yPred': ypred_of_ont_tool}

    basefn = 'x.x.GenomeWide' if not coordBedFileName else os.path.basename(coordBedFileName)

    os.makedirs(os.path.join(outdir, 'curve_data'), exist_ok=True)
    outfn = os.path.join(outdir, 'curve_data', f'{tagname}.{basefn}.curve_data.pkl')
    with open(outfn, 'wb') as handle:
        pickle.dump(curve_data, handle)

    ########################
    # correlation based stats:
    try:
        corrMix, pvalMix = pearsonr(ontFrequencies_4corr_mix, refFrequencies_4corr_mix)
    except:
        corrMix, pvalMix = (0, 0)

    try:
        corrAll, pvalAll = pearsonr(ontFrequencies_4corr_all, refFrequencies_4corr_all)
    except:
        corrAll, pvalAll = (0, 0)

    return accuracy, roc_auc, average_precision, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, corrMix, len(ontFrequencies_4corr_mix), corrAll, len(ontFrequencies_4corr_all), Csites_BGTruth, mCsites_BGTruth  # , leftovers,
    # leftovers1


# only care about AUC
def computePerReadStats_v2_for_roc_auc(ontCalls, bsReference, title, bedFile=False, ontCutt_perRead=1, ontCutt_4corr=4, secondFilterBed=False, secondFilterBed_4Corr=False):
    '''

    Only care about AUC data

    ontCalls - dictionary of CpG coordinates with their per read methylation call (1 or 0) // Format: {"chr\tstart\tend\n" : [list of methylation calls]}
    bsReference - dictionary of CpG coordinates with their methylation frequencies (range 0 - 1). This list is already prefiltered to meet minimal coverage (e.g. 4x) at this point. // Format: {"chr\tstart\tend\n" : methylation level (format: float (0-1))}
    title - prefix of the analysis, output plots etc. - should be as short as possible, but unique in context of other analyses
    bedFile - BED file which will be used to narrow down the list of CpGs for example to those inside CGIs or promoters etc.. By default "False" - which means no restrictions are done (i.e. genome wide)
    secondFilterBed - these should be CpGs covered in some reference list. Format: BED


    ============================================

    Basically i want to fill in the table below:
               Positive	Negative	Total
     Presence	a	        b	        a+b
     Absence	c	        d	        c+d
     Total	    a+c	        b+d	        a+b+c+d

    Nice summary also at wiki: https://en.wikipedia.org/wiki/F1_score
    , where "Positive" and "Negative" corresponds with ONT based observations, while "Presence" and "Absence" is based on BS-Seq

    '''

    switch = 0
    ontCalls_narrow = []
    logger.debug("in computePerReadStats2")
    if bedFile != False:
        ontCalls_bed = BedTool(dict2txt(ontCalls), from_string=True)
        ontCalls_bed = ontCalls_bed.sort()
        logger.debug("ontCalls_bed.sort()")

        infn = os.path.join(nanocompare_basedir, "reports", bedFile)
        # narrowBed = BedTool(bedFile)
        narrowBed = BedTool(infn)
        narrowBed = narrowBed.sort()
        logger.debug("narrowBed.sort()")

        ontCalls_intersect = ontCalls_bed.intersect(narrowBed, u=True, wa=True)
        ontCalls_narrow = txt2dict(ontCalls_intersect)
        suffix = bedFile
    else:
        switch = 1
        suffix = "GenomeWide"

    ## Second optional filter, organized in the same fashion as the first one. Designed to accomodate for example list of CpGs covered by second program
    secondSwitch = 0
    ontCalls_narrow_second = {}
    if secondFilterBed != False:
        infile = open(secondFilterBed, 'r')
        secondFilterDict = {}
        for row in infile:
            secondFilterDict[row] = 0
        infile.close()
        ontCalls_narrow_second = dict.fromkeys(set(ontCalls.keys()).intersection(set(secondFilterDict.keys())), 0)
    else:
        secondSwitch = 1
        suffix = "GenomeWide"

    # Second optional filter, shoudl be used in combination with second optional filter above
    secondSwitch_4corr = 0
    ontCalls_narrow_second_4corr = {}
    if secondFilterBed_4Corr != False:
        infile = open(secondFilterBed_4Corr, 'r')
        secondFilterDict = {}
        for row in infile:
            secondFilterDict[row] = 0
        infile.close()
        ontCalls_narrow_second_4corr = dict.fromkeys(set(ontCalls.keys()).intersection(set(secondFilterDict.keys())), 0)
    else:
        secondSwitch_4corr = 1
        suffix = "GenomeWide"

    suffixTMP = suffix.split("/")
    if len(suffixTMP) > 1:
        suffix = suffixTMP[-1]

    TP_5mC = FP_5mC = FN_5mC = TN_5mC = TP_5C = FP_5C = FN_5C = TN_5C = 0
    y = []
    scores = []

    ontSites = 0
    mCsites = 0
    Csites = 0
    referenceCpGs = 0

    for cpg_ont in ontCalls:
        ##### for per read stats:
        if cpg_ont in bsReference and len(ontCalls[cpg_ont]) >= ontCutt_perRead and (switch == 1 or cpg_ont in ontCalls_narrow) and (
                secondSwitch == 1 or cpg_ont in ontCalls_narrow_second):  # we should not take onCuttoffs for per read stats - shouldn't we? actually, we need to have the option to use this parameter, because at some point we may want to narrow down the per read stats to cover only the sites which were also covered by correlation with BS-Seq. Using this cutoff here is the easiest way to do just that
            #         if cpg_ont in bsReference and (switch == 1 or cpg_ont in ontCalls_narrow):
            if bsReference[cpg_ont] == 1 or bsReference[cpg_ont] == 0:  # we only consider absolute states here
                referenceCpGs += 1

                for ontCall in ontCalls[cpg_ont]:
                    if bsReference[cpg_ont] == 1:
                        mCsites += 1
                    if bsReference[cpg_ont] == 0:
                        Csites += 1
                    ontSites += 1

                    ### variables needed to compute precission, recall etc.:
                    if ontCall == 1 and bsReference[cpg_ont] == 1:  # true positive
                        TP_5mC += 1
                    elif ontCall == 1 and bsReference[cpg_ont] == 0:  # false positive
                        FP_5mC += 1
                    elif ontCall == 0 and bsReference[cpg_ont] == 1:  # false negative
                        FN_5mC += 1
                    elif ontCall == 0 and bsReference[cpg_ont] == 0:  # true negative
                        TN_5mC += 1

                    if ontCall == 0 and bsReference[cpg_ont] == 0:  # true positive
                        TP_5C += 1
                    elif ontCall == 0 and bsReference[cpg_ont] == 1:  # false positive
                        FP_5C += 1
                    elif ontCall == 1 and bsReference[cpg_ont] == 0:  # false negative
                        FN_5C += 1
                    elif ontCall == 1 and bsReference[cpg_ont] == 1:  # true negative
                        TN_5C += 1

                    ### AUC related:
                    scores.append(ontCall)
                    y.append(bsReference[cpg_ont])

    return scores, y
    ### compute all per read stats:

    #     Accuracy:
    try:
        accuracy = (TP_5mC + TN_5mC) / float(TP_5mC + FP_5mC + FN_5mC + TN_5mC)
    except ZeroDivisionError:
        accuracy = 0
    #     print("Accuracy: {0:1.4f}".format(accuracy))

    #     Positive predictive value (PPV), Precision = (TP) / E(Predicted condition positive)
    try:
        predicted_condition_positive_5mC = float(TP_5mC + FP_5mC)
        precision_5mC = TP_5mC / predicted_condition_positive_5mC
    except ZeroDivisionError:
        precision_5mC = 0
    #     print("Precision_5mC: {0:1.4f}".format(precision_5mC))
    try:
        predicted_condition_positive_5C = float(TP_5C + FP_5C)
        precision_5C = TP_5C / predicted_condition_positive_5C
    except ZeroDivisionError:
        precision_5C = 0

    #     print("Precision_5C: {0:1.4f}".format(precision_5C))

    #     True positive rate (TPR), Recall, Sensitivity, probability of detection = (TP) / (TP+FN)
    try:
        recall_5mC = TP_5mC / float(TP_5mC + FN_5mC)
    except ZeroDivisionError:
        recall_5mC = 0
    #     print("Recall_5mC: {0:1.4f}".format(recall_5mC))

    try:
        recall_5C = TP_5C / float(TP_5C + FN_5C)
    except ZeroDivisionError:
        recall_5C = 0
    #     print("Recall_5C: {0:1.4f}".format(recall_5C))

    #     F1 score:
    try:
        F1_5mC = 2 * ((precision_5mC * recall_5mC) / (precision_5mC + recall_5mC))
    except ZeroDivisionError:
        F1_5mC = 0
    #     print("F1 score_5mC: {0:1.4f}".format(F1_5mC))

    try:
        F1_5C = 2 * ((precision_5C * recall_5C) / (precision_5C + recall_5C))
    except ZeroDivisionError:
        F1_5C = 0
    #     print("F1 score_5C: {0:1.4f}".format(F1_5C))

    #     print("ontSites:", ontSites)

    ## plot AUC curve:
    fig = plt.figure(figsize=(5, 5), dpi=300)

    fprSwitch = 1
    try:
        fpr, tpr, _ = roc_curve(y, scores)
    except ValueError:
        logger.error(f"###\tERROR for roc_curve: y: {y}, scores: {scores}, \nother settings: {title}, {bedFile}, {secondFilterBed}, {secondFilterBed_4Corr}")
        fprSwitch = 0
        roc_auc = 0

    if fprSwitch == 1:
        roc_auc = auc(fpr, tpr)
        #     print("AUC: {0:1.4f}".format(roc_auc))
        #     print(title)

        lw = 2

        plt.plot(fpr, tpr, lw=lw, label='ROC curve (area = {0:.4f})'.format(roc_auc))
        plt.plot([0, 1], [0, 1], color='lightgrey', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(suffix)
        plt.legend(loc="lower right")

        fig.savefig("{}.{}.AUC.pdf".format(title, suffix), bbox_inches='tight', dpi=300)
        plt.close()

        ## Plot confusion matrix:
        cnf_matrix = confusion_matrix(y, scores)
        np.set_printoptions(precision=2)
        plt.figure(figsize=(5, 5))
        plot_confusion_matrix(cnf_matrix, classes=[0, 1], normalize=True, title=suffix)
        plt.savefig("{}.{}.ConfusionMatrix.pdf".format(title, suffix), bbox_inches='tight', dpi=300)
        plt.close()

        ## plot Precission-recall:
        average_precision = average_precision_score(y, scores)
        #     print('Average precision-recall score: {0:0.2f}'.format(average_precision))
        plt.figure(figsize=(5, 5))
        precision, recall, _ = precision_recall_curve(y, scores)
        # In matplotlib < 1.5, plt.fill_between does not have a 'step' argument
        step_kwargs = ({'step': 'post'}
                       if 'step' in signature(plt.fill_between).parameters
                       else {})
        plt.step(recall, precision, color='b', alpha=0.2,
                 where='post')
        plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)

        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('2-class Precision-Recall curve: AP={0:0.2f}\n{1}'.format(average_precision, suffix))
        plt.savefig("{}.{}.PrecissionRecall.pdf".format(title, suffix), bbox_inches='tight', dpi=300)
        plt.close()

    ########################
    # correlation based stats:
    try:
        corrMix, pvalMix = pearsonr(ontFrequencies_4corr_mix, refFrequencies_4corr_mix)
    except:
        corrMix, pvalMix = (0, 0)

    try:
        corrAll, pvalAll = pearsonr(ontFrequencies_4corr_all, refFrequencies_4corr_all)
    except:
        corrAll, pvalAll = (0, 0)

    return accuracy, roc_auc, precision_5C, recall_5C, F1_5C, Csites, precision_5mC, recall_5mC, F1_5mC, mCsites, referenceCpGs, corrMix, len(ontFrequencies_4corr_mix), corrAll, len(ontFrequencies_4corr_all)  # , leftovers, leftovers1


# def correlateResults_fromOrg(MyMethFreqInfile, RefMethFreqInfile, OutfilePrefix, MixedOnly = True, MyMethCutt = 4):
#     '''
#     Assumed format of the input files:
#     MyMethFreqInfile:
#     chromosome  start   end coverage    methylatedSites unmethylatedSites   methPercentage
#     chr1	27625914	27625915	1	1	0	1.0
#     chr1	143222602	143222603	21	20	1	0.9523809523809523

#     RefMethFreqInfile:
#     #chromosome     start   end     pmC     phmC    pC      err     qA      qB      qC      qD      N
#     chr16   14116   14116   0.25    0       0.75    0       1       3       0       0       0
#     chr16   14155   14155   0.75    0       0.25    0       3       1       0       0       0
#     chr16   16505   16505   1       0       0       0       4       0       0       0       0

#     '''

#     MyMethFreqDF = pd.read_csv(MyMethFreqInfile, sep="\t")
#     RefMethFreqDF = pd.read_csv(RefMethFreqInfile, sep="\t")

#     if MixedOnly == True:
#         RefMethFreqDF = RefMethFreqDF[(RefMethFreqDF.pmC > 0) & (RefMethFreqDF.pC > 0) & (RefMethFreqDF.pmC < 1) & (RefMethFreqDF.pC < 1)]

#     MyMethFreq = []
#     MyMethFreqDict = {}

#     RefMethFreq = []

#     for index, row in MyMethFreqDF.iterrows():
#         if int(row["coverage"]) >= MyMethCutt:
#             MyMethFreqDict[(row["chromosome"], row["start"], row["end"])] = float(row["methPercentage"])

#     for index, row in RefMethFreqDF.iterrows():
#         key = (row["#chromosome"], row["start"], row["end"])
#         if key in MyMethFreqDict:
#             MyMethFreq.append(MyMethFreqDict[key])
#             RefMethFreq.append(float(row["pmC"]))

#     corr, pval = stats.pearsonr(MyMethFreq, RefMethFreq)
#     return corr, pval, len(RefMethFreq)

def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          title="Confusion matrix",
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    #         print("Normalized confusion matrix")
    #     else:
    #         print('Confusion matrix, without normalization')

    #     print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.tight_layout()


def save_keys_to_single_site_bed(keys, outfn, callBaseFormat=0, outBaseFormat=0, nonstr='.'):
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


def save_call_or_bgtruth_to_bed(call, outfn, callBaseFormat=0, outBaseFormat=0):
    """
    Save to BED6 file format

    See also: https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    BED6: A BED file where each feature is described by chrom, start, end, name, score, and strand.
    For example: chr1 11873 14409 uc001aaa.3 0 +

    :param call:
    :param outfn:
    :return:
    """
    outfile = open(outfn, 'w')
    for key in call:
        # outfile.write(key[:-3])  # chr \t start \t end
        if outBaseFormat == 0:
            outfile.write(f'{key[0]}\t{key[1] - callBaseFormat}\t{key[1] - callBaseFormat + 1}')
        else:
            outfile.write(f'{key[0]}\t{key[1] - callBaseFormat + outBaseFormat}\t{key[1] - callBaseFormat + outBaseFormat}')

        ret = call[key]
        if type(ret) is list:
            if len(ret) != 2:
                raise Exception(f'Error when get value of key={key}, value={ret}')
            for k in ret:
                outfile.write(f'\t{k}')  # each other columns
            outfile.write(f'\t{key[2]}')  # strand info
            outfile.write(f'\n')
        else:
            outfile.write(f'\t.\t.\t{key[2]}\n')
    outfile.close()


def save_call_to_methykit_txt(call, outfn, callBaseFormat=0, is_cov=True):
    """
    For methykit analysis format

    is_cov  True if key->value value is [freq, cov]
            False if key->value value is list of [0000111, etc.]
    :param call:
    :param outfn:
    :return:
    """
    logger.debug(f'save_call_to_methykit_txt:{outfn}')
    outfile = open(outfn, 'w')
    outfile.write("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")
    print_first = True
    for key in call:
        if print_first:
            logger.debug(f'key={key}, value={call[key]}')
            print_first = False

        # strlist = key[:-1].split('\t')
        chr = key[0]
        base = key[1]
        strand = key[2]
        if strand == '+':
            strand = 'F'
        else:
            strand = 'R'

        if is_cov:
            freq = call[key][0]
            cov = call[key][1]
        else:
            try:
                freq = sum(call[key]) / float(len(call[key]))
                cov = len(call[key])
            except:
                logger.error(f'key={key}, value={call[key]}')
                exit(-1)
            pass
        outstr = f'{chr}.{base}\t{chr}\t{base}\t{strand}\t{cov}\t{freq:.2f}\t{1 - freq:.2f}\n'
        outfile.write(outstr)
    outfile.close()


def combine2programsCalls(calls1, calls2, outfileName=None):
    '''
    call1 and call2 should have the following format:
    {"chr\tstart\tend\n" : [list of methylation calls]}

    result: dictionary of shared sites in the same format as input
    '''
    tmp = dict.fromkeys(set(calls1.keys()).intersection(set(calls2.keys())), -1)
    if outfileName is not None:
        outfile = open(outfileName, 'w')
        for key in tmp:
            # outfile.write(key)
            outfile.write('\t'.join([key[0], str(key[1]), str(key[1]), '.', '.', key[2]]))
            outfile.write('\n')
        outfile.close()
    return tmp
    #     return dict.fromkeys(set(calls1.keys()).intersection(set(calls2.keys())), -1)
    # else:
    #     return dict.fromkeys(set(calls1.keys()).intersection(set(calls2.keys())), -1)


def combine2programsCalls_4Corr(calls1, calls2, cutt=3, outfileName=False, only_bgtruth=False, print_first=False):
    '''
    call1 and call2 should have the following format:
    {"chr\tstart\tend\n" : [list of methylation calls]}

    result: dictionary of shared sites in the same format as input with coverage of calls1 over calls 2 - use this with caution

    '''

    if only_bgtruth:  # value = [freq, cov]
        filteredCalls1 = {}
        for key in calls1:
            if calls1[key][1] >= cutt:
                filteredCalls1[key] = cutt
        return filteredCalls1

    # Next only allow value=[0 1 0 0 ...] or value=cov
    filteredCalls1 = {}

    for key in calls1:
        if print_first:
            logger.debug(key)
            print_first = False
        if isinstance(calls1[key], list):
            if len(calls1[key]) >= cutt:
                filteredCalls1[key] = cutt  # calls1[call]
        elif isinstance(calls1[key], int):
            if calls1[key] >= cutt:
                filteredCalls1[key] = cutt  # calls1[call]
    #         else:
    #             print("WARNING ### combine2programsCalls_4Corr ### calls1[call]: {}".format(type(calls1[call])))

    filteredCalls2 = {}

    for key in calls2:
        if print_first:
            logger.debug(key)
            print_first = False
        if isinstance(calls2[key], list):
            if len(calls2[key]) >= cutt:
                filteredCalls2[key] = cutt  # calls2[call]
        elif isinstance(calls2[key], int):
            if calls2[key] >= cutt:
                filteredCalls2[key] = cutt  # calls2[call]
    #         else:
    #             logger.error("WARNING ### combine2programsCalls_4Corr ### calls2[call]: {}".format(type(calls2[call])))
    # logger.debug(f'{len(filteredCalls1)}, {len(filteredCalls2)}')
    tmp = dict.fromkeys(set(filteredCalls1.keys()).intersection(set(filteredCalls2.keys())), cutt)
    if outfileName != False:
        outfile = open(outfileName, 'w')
        for key in tmp:
            outfile.write(key)
        outfile.close()

    # logger.debug("combine2programsCalls_4Corr DONE")
    return tmp
    # else:
    #     print("combine2programsCalls_4Corr DONE")
    #     return dict.fromkeys(set(filteredCalls1.keys()).intersection(set(filteredCalls2.keys())), cutt)


# Deprecated
def combine_ONT_and_BS(ontCalls, bsReference, analysisPrefix, narrowedCoordinates=False, ontCutt=4, secondFilterBed=False, secondFilterBed_4Corr=False):
    d = {"prefix"              : [],
            "coord"            : [],
            "accuracy"         : [],
            "roc_auc"          : [],
            "precision_5C"     : [],
            "recall_5C"        : [],
            "F1_5C"            : [],
            "Csites"           : [],
            "precision_5mC"    : [],
            "recall_5mC"       : [],
            "F1_5mC"           : [],
            "mCsites"          : [],
            "referenceCpGs"    : [],
            "corrMix"          : [],
            "Corr_mixedSupport": [],
            "corrAll"          : [],
            "Corr_allSupport"  : []}
    for coord in narrowedCoordinates:
        accuracy, roc_auc, precision_5C, recall_5C, F1_5C, Csites, precision_5mC, recall_5mC, F1_5mC, mCsites, referenceCpGs, corrMix, Corr_mixedSupport, corrAll, Corr_allSupport = computePerReadStats_deprecated(ontCalls, bsReference, analysisPrefix, bedFile=coord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)
        #         accuracy, roc_auc, precision_5C, recall_5C, F1_5C, Csites, precision_5mC, recall_5mC, F1_5mC, mCsites, referenceCpGs, corrMix, Corr_mixedSupport, corrAll, Corr_allSupport, leftovers, leftovers1 = computePerReadStats(ontCalls, bsReference, analysisPrefix, bedFile = coord, secondFilterBed = secondFilterBed, secondFilterBed_4Corr = secondFilterBed_4Corr)
        d["prefix"].append(analysisPrefix)
        d["coord"].append(coord)
        d["accuracy"].append(accuracy)
        d["roc_auc"].append(roc_auc)
        d["precision_5C"].append(precision_5C)
        d["recall_5C"].append(recall_5C)
        d["F1_5C"].append(F1_5C)
        d["Csites"].append(Csites)
        d["precision_5mC"].append(precision_5mC)
        d["recall_5mC"].append(recall_5mC)
        d["F1_5mC"].append(F1_5mC)
        d["mCsites"].append(mCsites)
        d["referenceCpGs"].append(referenceCpGs)
        d["corrMix"].append(corrMix)
        d["Corr_mixedSupport"].append(Corr_mixedSupport)
        d["corrAll"].append(corrAll)
        d["Corr_allSupport"].append(Corr_allSupport)

    #     df = pd.DataFrame.from_dict(d, orient='index')
    df = pd.DataFrame.from_dict(d)
    return df


def report_per_read_performance(ontCalls, bgTruth, analysisPrefix, narrowedCoordinatesList=None, ontCutt=4, secondFilterBed=None, secondFilterBed_4Corr=None, cutoff_meth=1.0, outdir=None, tagname=None, test=False):
    """
    New performance evaluation by Yang
    referenceCpGs is number of all CpGs that is fully-methylated (>=cutoff_meth) or unmethylated in BG-Truth

    :param ontCalls:
    :param bgTruth:
    :param analysisPrefix:
    :param narrowedCoordinatesList:
    :param ontCutt:
    :param secondFilterBed: Joined of 4 tools and bgtruth bed file
    :param secondFilterBed_4Corr: Joined of 4 tools and bgtruth bed file for cor analysis (with cutoff cov)
    :param cutoff_meth:
    :return:
    """
    d = defaultdict(list)

    for coord_fn in tqdm(narrowedCoordinatesList):
        accuracy, roc_auc, ap, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, corrMix, Corr_mixedSupport, corrAll, Corr_allSupport, cSites_BGTruth, mSites_BGTruth = \
            computePerReadStats(ontCalls, bgTruth, analysisPrefix, coordBedFileName=coord_fn, secondFilterBedFile=secondFilterBed,
                                secondFilterBed_4CorrFile=secondFilterBed_4Corr,
                                cutoff_meth=cutoff_meth, outdir=outdir, tagname=tagname)

        coord = os.path.basename(f'{coord_fn if coord_fn else "x.x.Genome-wide"}')

        d["prefix"].append(analysisPrefix)
        d["coord"].append(coord)
        d["Accuracy"].append(accuracy)
        d["Average-Precision"].append(ap)
        d["Macro-F1"].append(f1_macro)
        d["Micro-F1"].append(f1_micro)
        d["Macro-Precision"].append(precision_macro)
        d["Micro-Precision"].append(precision_micro)
        d["Macro-Recall"].append(recall_macro)
        d["Micro-Recall"].append(recall_micro)
        d["ROC-AUC"].append(roc_auc)
        d["Precision_5C"].append(precision_5C)
        d["Recall_5C"].append(recall_5C)
        d["F1_5C"].append(F1_5C)
        d["Csites_called"].append(cCalls)
        d["Csites"].append(cSites_BGTruth)
        d["Precision_5mC"].append(precision_5mC)
        d["Recall_5mC"].append(recall_5mC)
        d["F1_5mC"].append(F1_5mC)
        d["mCsites_called"].append(mCalls)
        d["mCsites"].append(mSites_BGTruth)
        d["referenceCpGs"].append(referenceCpGs)
        d["Corr_Mix"].append(corrMix)
        d["Corr_mixedSupport"].append(Corr_mixedSupport)
        d["Corr_All"].append(corrAll)
        d["Corr_allSupport"].append(Corr_allSupport)
        if test:
            break
    df = pd.DataFrame.from_dict(d)
    return df


def save_ontcalls_to_pkl(ontCalls, outfn):
    with open(outfn, 'wb') as handle:
        pickle.dump(ontCalls, handle)
    logger.info(f"save to {outfn}")


def load_ontcalls_pkl(infn):
    with open(infn, 'rb') as handle:
        ontCalls = pickle.load(handle)
    return ontCalls


def NonSingletonsScanner(referenceGenomeFile, outfileName_s, outfileName_ns):
    '''
    The output file is in 1-based coordinate system.
    '''
    reference = SeqIO.to_dict(SeqIO.parse(referenceGenomeFile, "fasta"))
    logger.debug("###\tNonSingletonsScanner: {} reference genome parsed".format(referenceGenomeFile))

    outfile_s = open(outfileName_s, "w")  # "s" stands for Singletons
    outfile_ns = open(outfileName_ns, "w")  # "s" stands for Non-Singletons

    for chromosome in list(reference.keys()):
        idxs = re.finditer('CG', str(reference[chromosome].seq).upper())

        singleton = -1  # 1 will stand for yes, 0 for no
        for idx in idxs:
            #             print(chromosome, idx, idx.start(), idx.end())
            if singleton == -1:
                s = idx.start() + 1  # here 8: mock1 <_sre.SRE_Match object; span=(8, 10), match='CG'> 8 10
                e = idx.end()  # here 10: mock1 <_sre.SRE_Match object; span=(8, 10), match='CG'> 8 10
                singleton = 1
            else:
                if (idx.start() - e) < 5:
                    # we just found a non-singleton. I.e. accordingly to the Nanopolish approach, CGs closer than 5bp, are considered as non-singletons
                    e = idx.end()
                    singleton = 0
                else:
                    # current CG is not part of non-singleton. It might mean that its not part of a big non-singleton or singleton upstream from it. We test which of these options below
                    if singleton == 1:
                        #                         print(chromosome, s, e, "SINGLETON")
                        outfile_s.write("{}\t{}\t{}\n".format(chromosome, s, e))
                    else:
                        #                         print(chromosome, s, e, "NON-SINGLETON")
                        outfile_ns.write("{}\t{}\t{}\n".format(chromosome, s, e))
                    s = idx.start() + 1
                    e = idx.end()
                    singleton = 1

        if singleton == 1:  # this code repetition takes care of the last instance in the long list of CG indexes
            #             print(chromosome, s, e, "SINGLETON")
            outfile_s.write("{}\t{}\t{}\n".format(chromosome, s, e))
        else:
            #             print(chromosome, s, e, "NON-SINGLETON")
            outfile_ns.write("{}\t{}\t{}\n".format(chromosome, s, e))

        logger.debug("###\tNonSingletonsScanner: chromosome {} processed".format(chromosome))

    outfile_s.close()
    outfile_ns.close()

    logger.debug("###\tNonSingletonsScanner: {} file processed".format(referenceGenomeFile))


def concat_dir_fn(outdir, fn):
    if outdir is not None:
        outfn = os.path.join(outdir, fn)
    else:
        outfn = fn
    return outfn


def nonSingletonsPostprocessing(referenceMeth, regionsBedFile, runPrefix, outdir, print_first=False):
    '''

    Concordant: all CpG in the region are 0, or 1. And they will be same 0 or 1.   00000000000    111111
    Discordant: all CpG in the region are 0, or 1. And they be both 0 and 1.   001000   110111
    Mix:        All CpG are valued 0.2, not 0 or 1.
    Other:      Not in case of before.

    This function will take the input *.bed file from "NonSingletonsScanner" funtion, which corresponds with non-singletons.
    Next it will separate them into concordant non-singletons (i.e. fully methylated or fully unmethylated), and disconcordant (those with at least one CpG fully methylated and at least one fully unmethylated), or fully mixed (i.e. all CpGs in non-singletons have methylation level >0 and < 100)
    This kind of preprocessing will have to be done for each studied library separately.
    '''
    logger.debug("nonSingletonsPostprocessing")
    refMeth = BedTool(dict2txt(referenceMeth), from_string=True)
    refMeth = refMeth.sort()

    infn = os.path.join(data_base_dir, 'genome-annotation', regionsBedFile)
    regions = BedTool(infn)
    regions = regions.sort()

    regions_refMeth = regions.intersect(refMeth, wa=True, wb=True)  # chr start end   chr start end strand

    regions_refMeth_dict = defaultdict(list)  # {}  # {regionCoords : [methylation percentage list]}

    for ovr in regions_refMeth:
        if print_first:
            logger.debug(f'ovr={ovr}')
            print_first = False
        # regionKey = "{}\t{}\t{}\n".format(ovr[0], ovr[1], ovr[2])
        regionKey = (ovr[0], int(ovr[1]), int(ovr[2]))  # chr  start  end
        # methKey = "{}\t{}\t{}\t{}\n".format(ovr[3], ovr[4], ovr[5], ovr[6])
        methKey = (ovr[3], int(ovr[4]), ovr[6])  # chr, start, strand
        # if regionKey not in regions_refMeth_dict:
        #     regions_refMeth_dict[regionKey] = []
        regions_refMeth_dict[regionKey].append(referenceMeth[methKey][0])

    outfile_prefix = regionsBedFile.replace(".bed", '')

    fn_concordant = f"{outdir}/{runPrefix}.{outfile_prefix}.concordant.bed"
    fn_discordant = f"{outdir}/{runPrefix}.{outfile_prefix}.discordant.bed"
    fn_fullyMixed = f"{outdir}/{runPrefix}.{outfile_prefix}.fullyMixed.bed"
    fn_other = f"{outdir}/{runPrefix}.{outfile_prefix}.other.bed"

    outfile_concordant = open(fn_concordant, "w")
    outfile_discordant = open(fn_discordant, "w")
    outfile_fullyMixed = open(fn_fullyMixed, "w")
    outfile_other = open(fn_other, "w")

    meth_cnt_dict = defaultdict(int)
    unmeth_cnt_dict = defaultdict(int)
    for region in regions_refMeth_dict:
        cntMeth = 0  # count how many methylated sites
        cntUnmeth = 0  # count how many unmeth
        fullMeth = 0  # indicate there is a fully meth case site in the region
        nullMeth = 0  # indicate there is a fully unmeth case site in the region
        mixMeth = 0  # indicate there is mixed ( 0.2 0.3) case site in the region
        for meth in regions_refMeth_dict[region]:
            if meth >= 1 - 1e-5:
                fullMeth = 1
                cntMeth += 1
            elif meth <= 1e-5:
                nullMeth = 1
                cntUnmeth += 1
            else:
                mixMeth = 1

        # region of chr start end
        region_txt = '\t'.join([region[0], str(region[1]), str(region[2])]) + '\n'

        if (fullMeth + nullMeth) == 1 and mixMeth == 0:
            #             print("Concordant")
            outfile_concordant.write(region_txt)
            meth_cnt_dict['Concordant'] += cntMeth
            unmeth_cnt_dict['Concordant'] += cntUnmeth
        elif (fullMeth + nullMeth) == 2 and mixMeth == 0:
            #             print("Discordant")
            outfile_discordant.write(region_txt)
            meth_cnt_dict['Discordant'] += cntMeth
            unmeth_cnt_dict['Discordant'] += cntUnmeth
        elif (fullMeth + nullMeth) == 0 and mixMeth == 1:
            #             print("mixed")
            outfile_fullyMixed.write(region_txt)
        else:
            #             print("What do we have here? ", fullMeth, nullMeth, mixMeth)
            outfile_other.write(f'{region_txt[:-1]}\t{fullMeth}\t{nullMeth}\t{mixMeth}\n')
    outfile_concordant.close()
    outfile_discordant.close()
    outfile_fullyMixed.close()
    outfile_other.close()
    logger.info(f'save to {[fn_concordant, fn_discordant, fn_fullyMixed, fn_other]}')

    logger.debug(f'meth_cnt={meth_cnt_dict}, unmeth_cnt={unmeth_cnt_dict}')

    ret = {'Concordant.5mC': meth_cnt_dict['Concordant'], 'Concordant.5C': unmeth_cnt_dict['Concordant'], 'Discordant.5mC': meth_cnt_dict['Discordant'], 'Discordant.5C': unmeth_cnt_dict['Discordant'], }
    return ret


def nonSingletonsPostprocessing2(referenceMeth, regionsBedFile, dataset, outdir=None, cutoff_meth=1.0):
    '''

    90% threshold, just set cutoff_meth = 0.9

    This function will take the input *.bed file from "NonSingletonsScanner" funtion, which corresponds with non-singletons.
    Next it will separate them into concordant non-singletons (i.e. fully methylated or fully unmethylated), and disconcordant (those with at least one CpG fully methylated and at least one fully unmethylated), or fully mixed (i.e. all CpGs in non-singletons have methylation level >0 and < 100)
    This kind of preprocessing will have to be done for each studied library separately.
    '''

    refMeth = BedTool(dict2txt(referenceMeth), from_string=True)
    refMeth = refMeth.sort()

    infn = os.path.join(nanocompare_basedir, "reports", regionsBedFile)
    # regions = BedTool(regionsBedFile)
    regions = BedTool(infn)

    regions = regions.sort()

    regions_refMeth = regions.intersect(refMeth, wa=True, wb=True)

    regions_refMeth_dict = {}  # {regionCoords : [methylation percentage list]}
    for ovr in regions_refMeth:
        regionKey = "{}\t{}\t{}\n".format(ovr[0], ovr[1], ovr[2])
        methKey = "{}\t{}\t{}\n".format(ovr[3], ovr[4], ovr[5])
        if regionKey not in regions_refMeth_dict:
            regions_refMeth_dict[regionKey] = []
        # regions_refMeth_dict value are a list of tuple as (sites, meth_val)
        regions_refMeth_dict[regionKey].append((methKey, referenceMeth[methKey]))

    outfile_prefix = regionsBedFile.replace(".bed", '')

    outfn = concat_dir_fn(outdir, f"{dataset}.{outfile_prefix}.concordant.bed")
    outfile_concordant = open(outfn, "w")

    outfn = concat_dir_fn(outdir, f"{dataset}.{outfile_prefix}.concordant.bed.sites")
    outfile_concordant_sites = open(outfn, "w")

    outfn = concat_dir_fn(outdir, f"{dataset}.{outfile_prefix}.discordant.bed")
    outfile_discordant = open(outfn, "w")

    outfn = concat_dir_fn(outdir, f"{dataset}.{outfile_prefix}.discordant.bed.sites")
    outfile_discordant_sites = open(outfn, "w")

    outfn = concat_dir_fn(outdir, f"{dataset}.{outfile_prefix}.fullyMixed.bed")
    outfile_fullyMixed = open(outfn, "w")

    outfn = concat_dir_fn(outdir, f"{dataset}.{outfile_prefix}.other.bed")
    outfile_other = open(outfn, "w")

    cnt_concordant_fully_methylated = 0
    cnt_concordant_unmethylated = 0

    cnt_concordant = 0
    cnt_discordant = 0
    cnt_fullyMixed = 0
    cnt_other = 0

    for region in regions_refMeth_dict:
        fullMeth = 0
        nullMeth = 0
        mixMeth = 0

        cur_sites = []

        # meth is a tuple of (sites, meth_val)
        for meth in regions_refMeth_dict[region]:
            if meth[1] >= (cutoff_meth - 1e-6):  # > 90% methlation level
                fullMeth = 1
                cur_sites.append(meth[0])
            elif meth[1] == 0:
                nullMeth = 1
                cur_sites.append(meth[0])
            else:
                mixMeth = 1

        if (fullMeth + nullMeth) == 1 and mixMeth == 0:
            #             print("Concordant")

            if fullMeth == 1:
                cnt_concordant_fully_methylated += 1
            else:
                cnt_concordant_unmethylated += 1

            outfile_concordant.write(region)
            cnt_concordant += 1

            for sites in cur_sites:
                outfile_concordant_sites.write(sites)

        elif (fullMeth + nullMeth) == 2:
            #             print("Discordant")
            outfile_discordant.write(region)
            cnt_discordant += 1

            for sites in cur_sites:
                outfile_discordant_sites.write(sites)

        elif (fullMeth + nullMeth) == 0 and mixMeth == 1:
            #             print("mixed")
            outfile_fullyMixed.write(region)
            cnt_fullyMixed += 1
        else:
            #             print("What do we have here? ", fullMeth, nullMeth, mixMeth)
            outfile_other.write("{}\t{}_{}_{}\n".format(region.strip(), fullMeth, nullMeth, mixMeth))
            cnt_other += 1

    outfile_concordant.close()
    outfile_discordant.close()
    outfile_fullyMixed.close()
    outfile_other.close()

    logger.info(f"DSNAME={dataset}, Non-singletons stats: concordant= {cnt_concordant} CpGs, discordant= {cnt_discordant} CpGs")
    logger.info(f"DSNAME={dataset}, Non-singletons stats: fully meth (Concordant)= {cnt_concordant_fully_methylated} CpGs, unmeth (Concordant)= {cnt_concordant_unmethylated} CpGs")

    ret = {'Nonsingletons'               : cnt_concordant + cnt_discordant,
            'Concordant'                 : cnt_concordant,
            'Discordant'                 : cnt_discordant,
            'Concordant.fully_methylated': cnt_concordant_fully_methylated,
            'Concordant.unmethylated'    : cnt_concordant_unmethylated}
    return ret


def singletonsPostprocessing(referenceMeth, singletonsBedFile, runPrefix, outdir, print_first=False):
    """

    Absoulute: In singleton area, CpG value is 0 or 1.

    BG-Truth using singleton as coordinate, output absolute, and mixed bed files.

    referenceMeth is key->value key=chr1  123  123  +, value=meth-freq

    This function will take the input *.bed file from "NonSingletonsScanner" funtion, which corresponds with singletons.
    Next it will separate them into absolute (i.e. fully methylated or fully unmethylated), and mixed (i.e. all CpGs in non-singletons have methylation level >0 and < 100)
    This kind of preprocessing will have to be done for each studied library separately.
    """

    logger.debug("singletonsPostprocessing")
    refMeth = BedTool(dict2txt(referenceMeth), from_string=True)
    refMeth = refMeth.sort()

    infnSingletons = os.path.join(data_base_dir, 'genome-annotation', singletonsBedFile)
    regionsSingletons = BedTool(infnSingletons)
    regionsSingletons = regionsSingletons.sort()

    refMeth_regions = refMeth.intersect(regionsSingletons, wa=True, u=True)

    outfile_prefix = singletonsBedFile.replace(".bed", '')

    absolute_fn = f"{outdir}/{runPrefix}.{outfile_prefix}.absolute.bed"
    mixed_fn = f"{outdir}/{runPrefix}.{outfile_prefix}.mixed.bed"
    outfile_absolute = open(absolute_fn, "w")
    outfile_mixed = open(mixed_fn, "w")

    absl_msite = 0
    absl_csite = 0
    for ovr in refMeth_regions:
        if print_first:
            logger.debug(f'ovr={ovr}')
            print_first = False
        # methKey = "{}\t{}\t{}\t{}\n".format(ovr[0], ovr[1], ovr[2], ovr[3])
        methKey = (ovr[0], int(ovr[1]), ovr[3])
        out_txt = '\t'.join([methKey[0], str(methKey[1]), str(methKey[1]), '.', '.', methKey[2]])

        if referenceMeth[methKey][0] > (1 - 1e-5) or referenceMeth[methKey][0] < 1e-5:  # TODO: due to strand info, we need to discuss how to add strand-info into the key
            outfile_absolute.write(f'{out_txt}\n')
            if referenceMeth[methKey][0] > (1 - 1e-5):
                absl_msite += 1
            else:
                absl_csite += 1

        else:
            outfile_mixed.write(f'{out_txt}\n')

    outfile_absolute.close()
    outfile_mixed.close()

    logger.info(f'In absolute singletons, msite={absl_msite:,}, csite={absl_csite:,}')

    logger.debug(f'save to {absolute_fn}')
    logger.debug(f'save to {mixed_fn}')

    ret = {'Singleton.5mc': absl_msite, 'Singleton.5C': absl_csite}
    return ret


# TODO: Deprecated, but can check before results if it is correct
def singletonsPostprocessing2(referenceMeth, singletonsBedFile, dataset, outdir=None, cutoff_meth=1.0):
    '''
    This function will take the input *.bed file from "NonSingletonsScanner" funtion, which corresponds with singletons.
    Next it will separate them into absolute (i.e. fully methylated or fully unmethylated), and mixed (i.e. all CpGs in non-singletons have methylation level >0 and < 100)
    This kind of preprocessing will have to be done for each studied library separately.
    '''

    refMeth = BedTool(dict2txt(referenceMeth), from_string=True)
    refMeth = refMeth.sort()

    infn = os.path.join(nanocompare_basedir, "reports", singletonsBedFile)

    # regions = BedTool(singletonsBedFile)
    regions = BedTool(infn)

    regions = regions.sort()

    refMeth_regions = refMeth.intersect(regions, wa=True, u=True)

    outfile_prefix = singletonsBedFile.replace(".bed", '')

    outfn = concat_dir_fn(outdir, f"{dataset}.{outfile_prefix}.absolute.bed")
    outfile_absolute = open(outfn, "w")

    outfn = concat_dir_fn(outdir, f"{dataset}.{outfile_prefix}.mixed.bed")
    outfile_mixed = open(outfn, "w")

    cnt_absolute = 0
    cnt_absolute_fully_methylated = 0
    cnt_absolute_unmethylated = 0

    cnt_mixed = 0
    for ovr in refMeth_regions:
        methKey = "{}\t{}\t{}\n".format(ovr[0], ovr[1], ovr[2])

        if referenceMeth[methKey] >= (cutoff_meth - 1e-6) or referenceMeth[methKey] == 0:
            outfile_absolute.write(methKey)
            cnt_absolute = cnt_absolute + 1
            if referenceMeth[methKey] == 0:
                cnt_absolute_unmethylated += 1
            else:
                cnt_absolute_fully_methylated += 1
        else:
            outfile_mixed.write(methKey)
            cnt_mixed = cnt_mixed + 1

    outfile_absolute.close()
    outfile_mixed.close()
    logger.info(f"DSNAME={dataset}, Singletons stats: absolute= {cnt_absolute} CpGs, mixed= {cnt_mixed} CpGs.")
    logger.info(f"DSNAME={dataset}, Singletons stats: unmethylated= {cnt_absolute_unmethylated} CpGs, fully_methylated= {cnt_absolute_fully_methylated} CpGs.")

    ret = {'Singletons'                  : cnt_absolute,
            'Singletons.fully_methylated': cnt_absolute_fully_methylated,
            'Singletons.unmethylated'    : cnt_absolute_unmethylated}

    return ret


def get_file_lines(infn):
    """
    Count number of lines in bedfile
    :param infn:
    :return:
    """
    num_lines = sum(1 for line in open(infn))
    return num_lines


def perRead2Frequency(inputDict, outfileName):
    """
    input file format = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}

    output file format:
    chr   start   end   methFreq    methylated_reads    unmethylated_reads
    """

    outfile = open(outfileName, 'w')

    for cpg in inputDict:
        outfile.write("{}\t{}\t{}\t{}\n".format(cpg.strip(), round(sum(inputDict[cpg]) / float(len(inputDict[cpg])) * 100), sum(inputDict[cpg]), len(inputDict[cpg])))

    outfile.close()
    logger.debug("###\tperRead2Frequency: completed frequency calculation for {} file".format(outfileName))


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


def load_deepmod_read_level_df(infn):
    """
    Load the DeepMod read level output results tsv into a dataframe
    :param infn:
    :return:
    """
    df = pd.read_csv(infn, sep='\t', header=None)
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


def get_dna_sequence_from_samfile(chr, start, end, bamfile):
    """
    Get the specific location DNA sequence at SAM files

    start is 0-based format

    Note: all functions of pysam need to be checked, see https://readthedocs.org/projects/pysam/downloads/pdf/latest/

    :param chr:
    :param start:
    :param end:
    :param bamfile:
    :return:
    """
    for read in bamfile.fetch(chr, start=start, end=end):

        alignedRefPositions = read.query_alignment_start
        # refStart = alignedRefPositions

        refStart = read.get_reference_positions()[0]

        # refSequence = read.get_reference_sequence()
        readSequence = read.query_alignment_sequence  # current use

        readSequence = read.query_sequence

        if readSequence is None:  # some read has no sequence, may return None, we only report the first read has sequence
            continue
        logger.debug(read.query_name)
        logger.debug(readSequence)

        # logger.debug(f'ref-start={refStart} len(align) = {len(readSequence)}, len(seq) = {len(readSequence1)} compare two:\n{readSequence}\n{readSequence1}')

        # logger.info(refStart)
        # logger.info(refSequence)
        # logger.info(readSequence)

        # logger.debug(readSequence[start - refStart:end - refStart])
        return readSequence[start - refStart:end - refStart]
    # if all reads return None, we report None
    return None


def get_dna_sequence_from_reference(chr, start, num_seq=5, ref_fasta=None):
    """
    Get the sequence from start-num_seq to start+num_seq, totally 2*num_seq+1

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


def get_ref_fasta():
    ref_fn = '/projects/li-lab/Ziwei/Nanopore/data/reference/hg38.fa'
    ref_fasta = SeqIO.to_dict(SeqIO.parse(open(ref_fn), 'fasta'))
    return ref_fasta


# Deal with bismark_bt2.CpG_report.txt.gz
def importGroundTruth_genome_wide_Bismark_Report(infn='/projects/li-lab/yang/results/2020-12-21/hl60-results-1/extractBismark/HL60_RRBS_ENCFF000MDF.Read_R1.Rep_2_trimmed_bismark_bt2.CpG_report.txt.gz', chr_col=0, start_col=1, strand_col=2, meth_col=3, unmeth_col=4, ccontect_col=5, covCutt=10, baseFormat=0, includeCov=True):
    """
    We are sure the input file is start using 1-based format.
    We use this format due to it contains strand info.
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
    :param covCutt:
    :param baseFormat:
    :param includeCov:
    :return:
    """
    df = pd.read_csv(infn, sep='\t', compression='gzip', header=None)

    df['cov'] = df.iloc[:, meth_col] + df.iloc[:, unmeth_col]
    df = df[df['cov'] >= covCutt]

    df = df[df.iloc[:, ccontect_col] == 'CG']

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

        if chr not in humanChrs:
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

    logger.info(f"###\timportGroundTruth_genome_wide_from_Bismark: loaded information for {len(cpgDict):,} CpGs from file {infn}")

    return cpgDict


def sanity_check_sequence(tlist):
    ret = get_ref_fasta()
    for site in tlist:
        ret_seq = get_dna_sequence_from_reference('chr8', site, ref_fasta=ret)
        logger.info(f'{site}:{ret_seq}')


def scatter_plot_x_y(xx, yy, tagname='', outdir=None):
    """
    :param xx:
    :param yy:
    :return:
    """
    if outdir is None:
        outdir = pic_base_dir

    import seaborn as sns
    import matplotlib.pyplot as plt

    plt.figure(figsize=(4, 4))

    ax = sns.scatterplot(x=xx, y=yy)
    ax.set_title(f"Scatter of {len(xx):,} CpG sites")

    plt.tight_layout()
    outfn = os.path.join(outdir, f'scatter-plot-{tagname}.jpg')
    plt.savefig(outfn, dpi=600, bbox_inches='tight')
    plt.show()
    plt.close()
    logger.info(f"save to {outfn}")

    pass


def scatter_plot_cov_compare_df(infn=None, df=None, outdir=None):
    if df is None:
        if infn is None:
            raise Exception("No data source specified.")

        if os.path.splitext(infn)[1] == '.csv':
            df = pd.read_csv(infn)
        elif os.path.splitext(infn)[1] == '.pkl':
            df = pd.read_pickle(infn)

    if outdir is None:
        outdir = pic_base_dir

    xx = df.iloc[:, 0]
    yy = df.iloc[:, 1]

    tagname = f'cov-of-{xx.name}-vs-{yy.name}'
    scatter_plot_x_y(xx, yy, tagname=tagname, outdir=outdir)


def scatter_analysis_cov(Tombo_calls1, Nanopolish_calls1, outdir, RunPrefix, tool1_name='Tombo', tool2_name='nanopolish'):
    tomboCpGs = set(Tombo_calls1.keys())
    intersectCpGs = tomboCpGs.intersection(Nanopolish_calls1.keys())
    key_list = []
    tombo_list = []
    tool_list = []
    for cpg in intersectCpGs:
        key_list.append(cpg[:-1])
        tombo_list.append(len(Tombo_calls1[cpg]))
        tool_list.append(len(Nanopolish_calls1[cpg]))

    scatterDf = pd.DataFrame(data={f'{tool1_name}-cov': tombo_list, f'{tool2_name}-cov': tool_list}, index=key_list)
    outfn = os.path.join(outdir, f'{RunPrefix}-{tool1_name}-{tool2_name}-scatter.csv')
    scatterDf.to_csv(outfn)

    outfn = os.path.join(outdir, f'{RunPrefix}-{tool1_name}-{tool2_name}-scatter.pkl')
    scatterDf.to_pickle(outfn)

    scatter_plot_cov_compare_df(df=scatterDf, outdir=outdir)


def get_high_cov_call1_low_cov_call2_df(call1, call2, call1_name='nanopolish', call2_name='tombo', low_cutoff=4, baseFormat=0):
    """
    Get the high cov of call1, but low cov of call2

    :param call1:
    :param call2:
    :param call1_name:
    :param call2_name:
    :param low_cutoff:
    :return:
    """
    interCpGs = set(call1.keys()).intersection(call2.keys())
    dataset = []
    for key in interCpGs:
        # kettext = key[:-1].split("\t")
        ret = {'chr': key[0], 'start': key[1], 'strand': key[2]}
        ret.update({call1_name: len(call1[key]), call2_name: len(call2[key])})
        dataset.append(ret)
    df = pd.DataFrame(dataset)

    logging.debug(f'Callings nocutoff {call1_name} and {call2_name} joined CpG sites={len(df)}')
    df = df[df[call2_name] < low_cutoff]
    df = df[df[call1_name] >= low_cutoff]
    logging.debug(f'After high {call1_name}-cov and low {call2_name}-cov: joined CpG sites={len(df)}')

    df['cov-diff'] = df[call1_name] - df[call2_name]
    df = df.sort_values(by='cov-diff', ascending=False)
    if baseFormat == 0:
        df['end'] = df['start'] + 1
    else:
        df['end'] = df['start']

    df = df[['chr', 'start', 'end', call1_name, call2_name, 'strand', 'cov-diff']]
    return df


def filter_non_human_chrs(cpgDict):
    """
    Filter out non-human chrs, only keep chr1-22 chrX chrY
    :param cpgDict:
    :return:
    """
    retDict = defaultdict(list)
    for key in cpgDict:
        if key[0] in humanChrs:
            retDict[key] = cpgDict[key]
    return retDict


def import_call(infn, callname, baseFormat=0):
    """
    Import fn based on callname and return key=(chr, start, strand) -> value=[0 0 1 1 1 1 ]
        call0   -   original results
    :param infn:
    :param callname:
    :return:
    """
    logger.debug(f"Start load {callname}")
    if callname == 'DeepSignal':
        calls0 = importPredictions_DeepSignal(infn, baseFormat=baseFormat)
    elif callname == 'Tombo':
        calls0 = importPredictions_Tombo(infn, baseFormat=baseFormat)
    elif callname == 'Nanopolish':
        calls0 = importPredictions_Nanopolish(infn, baseFormat=baseFormat, llr_cutoff=2.0)
    elif callname == 'DeepMod':
        calls0 = importPredictions_DeepMod_Read_Level(infn, baseFormat=baseFormat)
    elif callname == 'Megalodon':
        calls0 = importPredictions_Megalodon_Read_Level(infn, baseFormat=baseFormat)
    elif callname == 'DeepMod.Cluster':  # import DeepMod itself tool reports by cluster, key->value={'freq':68, 'cov':10}
        calls0 = importPredictions_DeepMod_clustered(infn, baseFormat=baseFormat)
    elif callname == 'DeepMod.C':  # import DeepMod itself tool reports by cluster, key->value={'freq':68, 'cov':10}
        calls0 = importPredictions_DeepMod(infn, baseFormat=baseFormat)
    else:
        raise Exception(f'Not support {callname} for file {infn} now')

    logger.debug(f'Import {callname} finished!\n')
    return calls0


def import_bgtruth(infn, encode, cov=10, baseFormat=0, includeCov=True):
    """
    Import bgtruth from file fn using encode, when use new dataset, MUST check input file start baseFormat and import functions are consistent!!!
    :param infn:
    :param encode:
    :param includeCov:  return [freq, cov] as value of each key=(chr, (int)start, strand)
    :return:
    """
    if encode == "encode":  # not used now, function need to check if use
        bgTruth = importGroundTruth_BedMethyl_from_Encode(infn, covCutt=cov, baseCount=baseFormat)
    elif encode == "oxBS_sudo":  # not used now, function need to check if use
        bgTruth = importGroundTruth_oxBS(infn, covCutt=cov, baseCount=baseFormat)
    elif encode == "bed":  # like K562 bg-truth, 0-based start input
        bgTruth = importGroundTruth_bed_file_format(infn, covCutt=cov, baseFormat=baseFormat, includeCov=includeCov)
    elif encode == "bismark_bedgraph":  # not used now, function need to check if use
        bgTruth = importGroundTruth_coverage_output_from_Bismark_BedGraph(infn, baseCount=baseFormat)
    elif encode == "bismark":  # for genome-wide Bismark Report ouput format results, such as HL60, etc. 1-based start input
        bgTruth = importGroundTruth_genome_wide_Bismark_Report(infn, covCutt=cov, baseFormat=baseFormat, includeCov=includeCov)
    else:
        raise Exception(f"encode={encode} is not supported yet, for inputfile={infn}")

    logger.debug(f'Import BG-Truth finished!\n')
    return bgTruth


def compare_cpg_key(item1, item2):
    """
    First compare chr, then start number
    usage:
    k1 = ('chr1', 123, '+')
    k2 = ('chr1', 124, '+')
    k3 = ('chr3', 1, '-')
    k4 = ('chr1', 124, '-')
    k5 = ('chr3', 1, '-')
    l = [k3, k1, k2, k4, k5]
    import functools
    l_sorted = sorted(l, key=functools.cmp_to_key(compare_cpg_key))
    logger.info(l_sorted)

    output:
    [('chr1', 123, '+'), ('chr1', 124, '+'), ('chr1', 124, '-'), ('chr3', 1, '-'), ('chr3', 1, '-')]

    :param item1:
    :param item2:
    :return:
    """
    chr1 = int(item1[0].replace('chr', '').replace('X', '23').replace('Y', '24'))
    chr2 = int(item2[0].replace('chr', '').replace('X', '23').replace('Y', '24'))

    # return chr1 < chr2 or (chr1 == chr2 and item1[1] < item2[1])

    if chr1 != chr2:
        return chr1 - chr2
    elif item1[1] != item2[1]:
        return item1[1] - item2[1]
    elif item1[2][0] != item2[2][0]:
        return ord(item1[2][0]) - ord(item2[2][0])
    return 0


if __name__ == '__main__':
    set_log_debug_level()

    # scatter_plot_cov_compare_df(infn='/projects/li-lab/yang/results/2020-12-28/K562_WGBS_Joined/K562_WGBS_Joinedtombo-nanopolish-scatter.pkl')

    # importGroundTruth_genome_wide_output_from_Bismark(covCutt=4)
    # sanity_check_sequence(tlist=[206309, 206316, 206318, 206494])
