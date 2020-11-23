#!/home/liuya/anaconda3/envs/nmf/bin/python

"""
Generate methy correlation plotting input data as a tsv

Sample usage:
    <pyfilename> $DeepSignal_calls $Tombo_calls $Nanopolish_calls \
			$DeepMod_calls $DeepMod_cluster_calls $bgTruth $parser $RunPrefix

All usedful functions are located in nanocompare.meth_stats.Universal_meth_stats_evaluation
"""
import sys

nanocompare_prj = "/projects/li-lab/yang/workspace/nano-compare/src"
sys.path.append(nanocompare_prj)

from nanocompare.meth_stats.meth_stats_common import *

from global_config import *

# from sklearn.utils.fixes import signature

import os
from sys import argv

# this script was based on this jupyter notebook: http://helix067:9912/notebooks/li-lab/NanoporeData/WR_ONT_analyses/ai/GitHub/nanoCompare/Scripts/UniversalMethStatsEvaluation.ipynb
# copied Wed Jul 10 14:27:33 EDT 2019
# example command run: python /projects/li-lab/NanoporeData/WR_ONT_analyses/ai/GitHub/nanoCompare/Scripts/UniversalMethStatsEvaluation.standalone_01.py $DeepSignal_calls $Tombo_calls $Nanopolish_calls $bgTruth $RunPrefix $parser

# where variable names came from:
# DeepSignal_calls = importPredictions_DeepSignal(argv[1])#"/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepSignal_runs/K562/K562_DeepSignal.MethCalls.Joined.chr20.tsv")
# Tombo_calls = importPredictions_Tombo(argv[2])#"/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/K562/K562_Tombo.perReadsStats.chr20.bed")
# Nanopolish_calls = importPredictions_Nanopolish_2(argv[3])#"/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.chr_chr20.tsv", IncludeNonSingletons = True)
# bgTruth = importGroundTruth_BedMethyl_from_Encode(argv[4])#"/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/K562/ENCSR765JPC_WGBS_hg38/ENCFF721JMB.bed", chrFilter="chr20")
# RunPrefix = argv[5]#"K562_WGBS_rep_ENCFF721JMB/K562_WGBS_rep_ENCFF721JMB"
# parser = "encode"

####################
# from nanocompare.meth_stats.Universal_meth_stats_evaluation import importPredictions_Tombo, importPredictions_DeepSignal, coverageFiltering, importPredictions_Nanopolish_2, importPredictions_DeepMod, DeepMod_clusteredResultParsing, importGroundTruth_BedMethyl_from_Encode, importGroundTruth_oxBS, importGroundTruth_coverage_output_from_Bismark, \
#     importGroundTruth_coverage_output_from_Bismark_BedGraph

#
# def report2dict(cr):
#     # solution of  "jolespin commented on Jan 19, 2017", + the updated by "HyungSeokPark commented on Mar 20, 2018"
#     # https://github.com/scikit-learn/scikit-learn/issues/7845
#     # i needed that because classification_report function does not recognize the "output_dict" parameter
#
#     # Parse rows
#     tmp = list()
#     for row in cr.split("\n"):
#         parsed_row = [x for x in row.split("  ") if len(x) > 0]
#         if len(parsed_row) > 0:
#             tmp.append(parsed_row)
#
#     # Store in dictionary
#     measures = tmp[0]
#     print(measures)
#
#     D_class_data = {}  # defaultdict(dict)
#     for row in tmp[1:]:
#         class_label = row[0].strip()
#         for j, m in enumerate(measures):
#             D_class_data[class_label][m.strip()] = float(row[j + 1].strip())
#     return D_class_data

#
# def importPredictions_NanoXGBoost(infileName, chr_col=0, start_col=1, meth_col=4, baseCount=1):
#     '''
#     Note that the function requires per read stats, not frequencies of methylation.
#     !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!
#
#     ### Parameters of the function:
#     chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
#     start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
#     meth_col - name (as header) of the column with methylation call (integer expected). If "header" variable == False, give integer number of the column.
#     baseCount - 0 or 1, standing for 0-based or 1-based, respectively
#     #header - True or False.
#
#     ### Example input format from my NanoXGBoost model:
#     chr20   15932019        15932019        /home/rosikw/fastscratch/APL_newSept/0/GXB01186_20180508_FAH83098_GA10000_sequencing_run_180508_18_li_001_GXB01186_001_12562_read_35498_ch_259_strand.fast5     1
#     chr20   15932898        15932898        /home/rosikw/fastscratch/APL_newSept/0/GXB01186_20180508_FAH83098_GA10000_sequencing_run_180508_18_li_001_GXB01186_001_12562_read_35498_ch_259_strand.fast5     0
#     chr20   15932997        15932997        /home/rosikw/fastscratch/APL_newSept/0/GXB01186_20180508_FAH83098_GA10000_sequencing_run_180508_18_li_001_GXB01186_001_12562_read_35498_ch_259_strand.fast5     1
#     chr20   15933019        15933019        /home/rosikw/fastscratch/APL_newSept/0/GXB01186_20180508_FAH83098_GA10000_sequencing_run_180508_18_li_001_GXB01186_001_12562_read_35498_ch_259_strand.fast5     1
#     Note: here there are no headers (probably this will change in the future)
#
#     ### Example input format from Nanopolish:
#     chromosome      start   end     read_name       log_lik_ratio   log_lik_methylated      log_lik_unmethylated    num_calling_strands     num_cpgs        sequence
#     chr20   106142  106142  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    2.32    -93.28  -95.61  1       1       CTCAACGTTTG
#     chr20   106226  106226  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    1.46    -193.91 -195.37 1       1       TGGCACGTGGA
#     chr20   104859  104859  107a0850-7500-443c-911f-4857424c889c    4.36    -163.26 -167.62 1       1       ATTCCCGAGAG
#     Note: Nanopolish output have header, yet the function needs to be universal enough to handle both cases (including whatever awaits in case of NanoMod, DeepSignal etc.)
#
#     ### Output format:
#     result = {"chr\tstart\tend\n" : [list of methylation calls]}
#     output coordinates are 1-based genome coordinates.
#
#     ============
#
#     Future changes:
#     This function cannot be for everything. Output from nanopolish is too different from mine (e.g. singletons and non-singletons), so it will need an independent parser. Maybe the same will go for other programs?
#     Nevertheless, all functions will have the same output, which is the most important part.
#     '''
#
#     infile = open(infileName, "r")
#     cpgDict = {}
#     count = 0
#
#     for row in infile:
#         tmp = row.strip().split("\t")
#         if baseCount == 1:
#             start = int(tmp[start_col])
#         elif baseCount == 0:
#             start = int(tmp[start_col]) + 1
#         else:
#             print("###\timportPredictions_NanoXGBoost InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#             exit
#         #         key = (tmp[chr_col], start)
#         key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)
#         if key not in cpgDict:
#             cpgDict[key] = []
#         cpgDict[key].append(int(tmp[meth_col]))
#         count += 1
#
#     infile.close()
#
#     print("###\timportPredictions_NanoXGBoost SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
#     return cpgDict
#
#
# def importPredictions_Nanopolish(infileName, chr_col=0, start_col=1, log_lik_ratio_col=4, sequence_col=-1, baseCount=0, logLikehoodCutt=2.5):
#     '''
#     !!! This function will be needed for NanoCompare project !!!
#
#     ### Example input format from Nanopolish:
#     chromosome      start   end     read_name       log_lik_ratio   log_lik_methylated      log_lik_unmethylated    num_calling_strands     num_cpgs        sequence
#     chr20   106142  106142  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    2.32    -93.28  -95.61  1       1       CTCAACGTTTG
#     chr20   106226  106226  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    1.46    -193.91 -195.37 1       1       TGGCACGTGGA
#     chr20   104859  104859  107a0850-7500-443c-911f-4857424c889c    4.36    -163.26 -167.62 1       1       ATTCCCGAGAG
#
#     ### Output format:
#     result = {"chr\tstart\tend\n" : [list of methylation calls]}
#     output coordinates are 1-based genome coordinates.
#
#     ###############
#
#     almost ready to use code is here: http://helix122:9912/edit/li-lab/NanoporeData/WR_ONT_analyses/ai/APL_nanopolishStats/automatedSingleReadPrecission_3.py
#     in function called "nanopolishMethCallsParser_2_NanopolishCalled" - I will add this later (or earlier if asked for it:)
#
#     Their current script for handling with conversion of calls to frequencies:
#     https://github.com/jts/nanopolish/blob/master/scripts/calculate_methylation_frequency.py
#
#     It looks that we both do the same, or not?
#
#
#     '''
#
#     cpgDict = {}
#     count = 0
#     infile = open(infileName, 'r')
#
#     for row in infile:
#         tmp = row.strip().split("\t")
#         if tmp[chr_col] != "chromosome":
#             if int(tmp[-2]) == 1:  # we have singleton, i.e. only one CpG within the area
#                 if baseCount == 0:
#                     key = "{}\t{}\t{}\n".format(tmp[chr_col], int(tmp[start_col]) + 1, int(tmp[start_col]) + 1)
#                 elif baseCount == 1:
#                     key = "{}\t{}\t{}\n".format(tmp[chr_col], int(tmp[start_col]), int(tmp[start_col]))
#                 else:
#                     print("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#                     exit
#
#                 if float(tmp[log_lik_ratio_col]) >= logLikehoodCutt:
#                     if key not in cpgDict:
#                         cpgDict[key] = []
#                     cpgDict[key].append(1)
#                     count += 1
#                 elif float(tmp[log_lik_ratio_col]) <= -logLikehoodCutt:
#                     if key not in cpgDict:
#                         cpgDict[key] = []
#                     cpgDict[key].append(0)
#                     count += 1
#
#             else:  # we deal with non-singleton
#                 firstCpgLoc = int(tmp[start_col]) - 5
#                 sequence = tmp[sequence_col]
#                 for cpg in re.finditer("CG", sequence):
#                     cpgStart = cpg.start() + firstCpgLoc + 1
#                     cpgEnd = cpgStart
#
#                     if baseCount == 0:
#                         key = "{}\t{}\t{}\n".format(tmp[chr_col], cpgStart + 1, cpgStart + 1)
#                     elif baseCount == 1:
#                         key = "{}\t{}\t{}\n".format(tmp[chr_col], cpgStart, cpgStart)
#                     else:
#                         print("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#                         exit
#
#                     if float(tmp[log_lik_ratio_col]) >= logLikehoodCutt:
#                         if key not in cpgDict:
#                             cpgDict[key] = []
#                         cpgDict[key].append(1)
#                         count += 1
#                     elif float(tmp[log_lik_ratio_col]) <= -logLikehoodCutt:
#                         if key not in cpgDict:
#                             cpgDict[key] = []
#                         cpgDict[key].append(0)
#                         count += 1
#
#     infile.close()
#
#     print("###\timportPredictions_Nanopolish SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
#     return cpgDict

#
# # Deprecated
# def e_importPredictions_Nanopolish_2(infileName, baseCount=0, logLikehoodCutt=2.5, IncludeNonSingletons=True):
#     '''
#     Nanopolish Parser function based on parsing script from Nanopolish:
#     https://github.com/jts/nanopolish/blob/master/scripts/calculate_methylation_frequency.py
#     Code was downloaded from their Github and modified for the purpose of this project on April 30th 2019.
#
#     Generally it gives exactly the same results as my own, but i think its better to use their code, so that nobody would be able to say that we did something differently
#
#
#     !!! This function will be needed for NanoCompare project !!!
#
#     ### Example input format from Nanopolish:
#     chromosome      start   end     read_name       log_lik_ratio   log_lik_methylated      log_lik_unmethylated    num_calling_strands     num_cpgs        sequence
#     chr20   106142  106142  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    2.32    -93.28  -95.61  1       1       CTCAACGTTTG
#     chr20   106226  106226  0b42b84b-c2a7-481b-be33-c555eb2e1fcf    1.46    -193.91 -195.37 1       1       TGGCACGTGGA
#     chr20   104859  104859  107a0850-7500-443c-911f-4857424c889c    4.36    -163.26 -167.62 1       1       ATTCCCGAGAG
#
#     ### Output format:
#     result = {"chr\tstart\tend\n" : [list of methylation calls]}
#     output coordinates are 1-based genome coordinates.
#
#
#
#     '''
#     count = 0
#     cpgDict = {}
#
#     infile = open(infileName, 'r')
#     csv_reader = csv.DictReader(infile, delimiter='\t')
#
#     for record in csv_reader:
#         #         print(record)
#         num_sites = int(record['num_cpgs'])
#         llr = float(record['log_lik_ratio'])
#
#         # Skip ambiguous call
#         if abs(llr) < logLikehoodCutt:
#             continue
#         sequence = record['sequence']
#
#         is_methylated = int(llr > 0)
#
#         # if this is a multi-cpg group and split_groups is set, break up these sites
#         if IncludeNonSingletons and num_sites > 1:
#             c = str(record['chromosome'])
#             s = int(record['start'])
#             e = int(record['end'])
#
#             # find the position of the first CG dinucleotide
#             sequence = record['sequence']
#             cg_pos = sequence.find("CG")
#             first_cg_pos = cg_pos
#             while cg_pos != -1:
#                 #                 key = (c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
#                 if baseCount == 0:
#                     key = "{}\t{}\t{}\n".format(c, s + cg_pos - first_cg_pos + 1, s + cg_pos - first_cg_pos + 1)
#                 elif baseCount == 1:
#                     key = "{}\t{}\t{}\n".format(c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
#                 else:
#                     print("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#                     exit
#
#                 if key not in cpgDict:
#                     cpgDict[key] = []
#                 cpgDict[key].append(is_methylated)
#                 count += 1
#
#                 cg_pos = sequence.find("CG", cg_pos + 1)
#         else:
#             if baseCount == 0:
#                 key = "{}\t{}\t{}\n".format(str(record['chromosome']), int(record['start']) + 1, int(record['end']) + 1)
#             elif baseCount == 1:
#                 key = "{}\t{}\t{}\n".format(str(record['chromosome']), int(record['start']), int(record['end']))
#             else:
#                 print("###\timportPredictions_Nanopolish InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#                 exit
#
#             if key not in cpgDict:
#                 cpgDict[key] = []
#             cpgDict[key].append(is_methylated)
#             count += 1
#
#     print("###\timportPredictions_Nanopolish SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
#     return cpgDict
#
#
# # Deprecated
# def e_importPredictions_DeepSignal(infileName, chr_col=0, start_col=1, meth_col=8, baseCount=0):
#     '''
#     Note that the function requires per read stats, not frequencies of methylation.
#     !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!
#
#     ### Parameters of the function:
#     chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
#     start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
#     meth_col - name (as header) of the column with methylation call (integer expected). If "header" variable == False, give integer number of the column. *7 is probability, 8 is binary call.
#     baseCount - 0 or 1, standing for 0-based or 1-based, respectively
#     #header - True or False.
#
#     ### Example input format from DeepSignal:
#     chr11   127715633       -       7370988 791c33e4-5b63-4cce-989c-186aff79db9b    t       0.16227172      0.83772826      1       AGGAAAATCGCTTGAAC
#     chr11   127715585       -       7371036 791c33e4-5b63-4cce-989c-186aff79db9b    t       0.69724774      0.3027523       0       TGCCACTGCGCTCTAGC
#     chr11   127715554       -       7371067 791c33e4-5b63-4cce-989c-186aff79db9b    t       0.786389        0.21361093      0       CAGAACTCCGTCTCAAA
#     chr11   127715423       -       7371198 791c33e4-5b63-4cce-989c-186aff79db9b    t       0.5939311       0.40606892      0       TTTTGTAGCGTTGTACA
#
#     ### Input file format description:
#     - chrom: the chromosome name
#     - pos: 0-based position of the targeted base in the chromosome
#     - strand: +/-, the aligned strand of the read to the reference
#     - pos_in_strand: 0-based position of the targeted base in the aligned strand of the chromosome
#     - readname: the read name
#     - read_strand: t/c, template or complement
#     - prob_0: [0, 1], the probability of the targeted base predicted as 0 (unmethylated)
#     - prob_1: [0, 1], the probability of the targeted base predicted as 1 (methylated)
#     - called_label: 0/1, unmethylated/methylated
#     - k_mer: the kmer around the targeted base
#
#     ### Output format:
#     result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
#     output coordinates are 1-based genome coordinates.
#
#     ** by default if this probability will be higher than 0.5, DeppSignal will tell that this is methylated site, lower, unmethylated
#
#     ============
#
#     '''
#
#     infile = open(infileName, "r")
#     cpgDict = {}
#     count = 0
#
#     for row in infile:
#         tmp = row.strip().split("\t")
#         if baseCount == 1:
#             start = int(tmp[start_col])
#         elif baseCount == 0:
#             start = int(tmp[start_col]) + 1
#         else:
#             print("###\timportPredictions_DeepSignal InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#             exit
#         #         key = (tmp[chr_col], start)
#         key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)
#         if key not in cpgDict:
#             cpgDict[key] = []
#         #         cpgDict[key].append(float(tmp[meth_col])) ##### uncomment this line to get probabilities instead of final, binary calls
#         cpgDict[key].append(int(tmp[meth_col]))
#         count += 1
#
#     infile.close()
#
#     print("###\timportPredictions_DeepSignal SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
#     return cpgDict
#
#
# # Deprecated
# def e_importPredictions_Tombo(infileName, chr_col=0, start_col=1, meth_col=4, baseCount=0, cutoff=2.5):
#     '''
#     Note that the function requires per read stats, not frequencies of methylation.
#     !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!
#
#     ### Parameters of the function:
#     chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
#     start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
#     meth_col - name (as header) of the column with methylation call (integer expected). If "header" variable == False, give integer number of the column.
#     baseCount - 0 or 1, standing for 0-based or 1-based, respectively
#     cutoff - sumilarly as in case of Nanopolish, here we have cutoff for the value from the statistical test. From this conversations (https://github.com/nanoporetech/tombo/issues/151), I know this value is by default 2.5.
#
#     ### Example input format from Tombo (v1.5):
#     chr1    66047   66047   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    6.057825558813564       +
#     chr1    66053   66053   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    -0.3359579051241508     +
#     chr1    66054   66054   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    0.1202407639936725      +
#     chr1    66055   66055   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    2.1077369345267907      +
#     chr1    66076   66076   ed4a12ec-e03a-4a0a-9d08-acf3c0ee11d4    0.8979673996582611      +
#
#     ### Output format:
#     result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
#     output coordinates are 1-based genome coordinates.
#
#     ** by default if this probability will be higher than 0.5, DeppSignal will tell that this is methylated site, lower, unmethylated
#
#     ============
#
#     '''
#
#     infile = open(infileName, "r")
#     cpgDict = {}
#     count = 0
#
#     for row in infile:
#         tmp = row.strip().split("\t")
#         if baseCount == 1:
#             start = int(tmp[start_col])
#         elif baseCount == 0:
#             start = int(tmp[start_col]) + 1
#         else:
#             print("###\timportPredictions_Tombo InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#             exit
#         #         key = (tmp[chr_col], start)
#         key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)
#
#         methCall = float(tmp[meth_col])
#         switch = 0
#         #         if methCall > cutoff:
#         if methCall < -cutoff:
#             methCall = 1
#             switch = 1
#         #         elif methCall < -cutoff:
#         elif methCall > cutoff:
#             methCall = 0
#             switch = 1
#
#         if switch == 1:
#             if key not in cpgDict:
#                 cpgDict[key] = []
#             cpgDict[key].append(methCall)
#         count += 1
#
#     infile.close()
#
#     print("###\timportPredictions_Tombo SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
#     return cpgDict
#
#
# # Deprecated
# def e_importPredictions_DeepMod(infileName, chr_col=0, start_col=1, meth_reads_col=12, coverage_col=10, clusteredResult=False, clustered_meth_freq_col=13, baseCount=0):
#     '''
#     Note that the function requires per read stats, not frequencies of methylation.
#     !!! Also, this note is now optimized for my NanoXGBoost output - nothing else. !!!
#
#     ### Parameters of the function:
#     chr_col - name (as header) of the column with chromosome. If "header" variable == False, give integer number of the column.
#     start_col - name (as header) of the column with start of CpG. If "header" variable == False, give integer number of the column.
#     meth_reads_col - name (as header) of the column with number of methylated reads mapped.
#     coverage_col - name (as header) of the column with coverage of the site.
#     [[TO DO]] clusteredResult - True / False. Input file is in the "clustered" format (additional post-processing step). False (default option) - standard output with calls.
#     [[TO DO]] clustered_meth_freq_col - column with the methylation frequency after additional postprocessing step.
#     baseCount - 0 or 1, standing for 0-based or 1-based, respectively
#
#     ### Example input format from DeepMod (standard):
#     chr2 110795922 110795923 C 4 -  110795922 110795923 0,0,0 4 75 3
#     chr2 110795929 110795930 C 3 -  110795929 110795930 0,0,0 3 66 2
#     chr2 110796453 110796454 C 4 -  110796453 110796454 0,0,0 4 25 1
#
#     Description (https://github.com/WGLab/DeepMod/blob/master/docs/Results_explanation.md):
#     The output is in a BED format like below. The first six columns are Chr,
#     Start pos, End pos, Base, Capped coverage, and Strand, and the last three
#     columns are Real coverage, Mehylation percentage and Methylation coverage.
#
#
#     ### Example input format from DeepMod (clustered - following Step 4 from "Example 3: Detect 5mC on Na12878" section; https://github.com/WGLab/DeepMod/blob/master/docs/Reproducibility.md):
#     chr2 241991445 241991446 C 3 -  241991445 241991446 0,0,0 3 100 3 69
#     chr2 241991475 241991476 C 3 -  241991475 241991476 0,0,0 3 33 1 75
#     chr2 241991481 241991482 C 2 -  241991481 241991482 0,0,0 2 50 1 76
#
#     Note: it is space-separated, not tab-separated file
#
#     ### Output format:
#     result = {"chr\tstart\tend\n" : [list of methylation calls (as a probability of methylation call**)]}
#     output coordinates are 1-based genome coordinates.
#
#     ============
#
#     '''
#
#     infile = open(infileName, "r")
#     cpgDict = {}
#     count = 0
#
#     for row in infile:
#         tmp = row.strip().split(" ")
#         if baseCount == 1:
#             start = int(tmp[start_col])
#         elif baseCount == 0:
#             start = int(tmp[start_col]) + 1
#         else:
#             print("###\timportPredictions_DeepMod InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#             exit
#         #         key = (tmp[chr_col], start)
#         key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)
#
#         methCalls = int(tmp[meth_reads_col])
#         coverage = int(tmp[coverage_col])
#
#         methCallsList = [1] * methCalls + [0] * (coverage - methCalls)
#         cpgDict[key] = methCallsList
#
#         count += len(methCallsList)
#
#     infile.close()
#
#     print("###\timportPredictions_DeepMod SUCCESS: {} methylation calls mapped to {} CpGs from {} file".format(count, len(cpgDict), infileName))
#     return cpgDict

#
# def e_importGroundTruth_oxBS(infileName, chr_col='#chromosome', start_col='start', meth_col="pmC", covCutt=4, baseCount=1, chrFilter=False):
#     '''
#     Note that this function was optimized to parse the data from my oxBS results. More specifically, the (sudo)format which I have created to be able to have the information from both BS and corresponding oxBS-seq.
#
#     ### Description of the columns in this format:
#     1. chromosome
#     2. start
#     3. end
#     4. pmC - percentage of methylated
#     5. artifact = 0
#     6. pC - percentage of unmethylated
#     7. artifact = 0
#     8. qA - No. of 5mC + 5hmC reads (oxBS-seq based // quadrant A)
#     9. qB - No. of C reads (oxBS-seq based // quadrant B)
#     10. artifact = 0
#     11. artifact = 0
#     12. artifact = 0
#
#     e.g.:
#     #chromosome	start	end	pmC	phmC	pC	err	qA	qB	qC	qD	N
#     chr1	10662	10663	1.0	0.0	0.0	0	4	0	4	0	1.0
#     chr1	10665	10666	1.0	0.0	0.0	0	6	0	4	0	1.5
#     chr1	10667	10668	1.0	0.0	0.0	0	6	0	4	0	1.5
#
#     ### Output files:
#     cpgDict = {"chr\tstart\tend\n" : methylation level (format: float (0-1))} # have all cpgs with specified cutoff (not only 100% 5mC or 100% 5C)
#     output coordinates are 1-based genome coordinates.
#
#     '''
#
#     cpgDict = {}
#
#     infile = open(infileName, 'r')
#     csvfile = csv.DictReader(infile, delimiter='\t')
#     for row in csvfile:
#         if baseCount == 1:
#             start = int(row[start_col])
#         elif baseCount == 0:
#             start = int(row[start_col]) + 1
#         else:
#             print("###\timportGroundTruth_oxBS InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#             exit
#
#         if chrFilter == False or chrFilter == row[chr_col]:
#             cov = int(row["qA"]) + int(row["qB"])
#             if cov >= covCutt:
#                 #                 key = (row[chr_col], start)
#                 key = "{}\t{}\t{}\n".format(row[chr_col], start, start)
#                 if key not in cpgDict:
#                     cpgDict[key] = float(row['pmC'])
#                 else:
#                     print("###\timportGroundTruth_oxBS SanityCheckError: One CpG should not have more than 1 entry")
#
#     infile.close()
#     return cpgDict
#
#
# def e_importGroundTruth_BedMethyl_from_Encode(infileName, chr_col=0, start_col=1, meth_col=10, cov_col=9, covCutt=10, baseCount=0, chrFilter=False):
#     '''
#
#     ### Description of the columns in this format (https://www.encodeproject.org/data-standards/wgbs/):
#
#     1. Reference chromosome or scaffold
#     2. Start position in chromosome (0-based position)
#     3. End position in chromosome (1-based position)
#     4. Name of item
#     5. Score from 0-1000. Capped number of reads
#     6. Strandedness, plus (+), minus (-), or unknown (.)
#     7. Start of where display should be thick (start codon)
#     8. End of where display should be thick (stop codon)
#     9. Color value (RGB)
#     10. Coverage, or number of reads
#     11. Percentage of reads that show methylation at this position in the genome
#
#     I think that this output is 0-based.
#
#     e.g.:
#     chr17   115761  115762  .       0       -       115761  115762  0,255,0 0       0
#     chr17   116083  116084  .       3       +       116083  116084  0,255,0 3       0
#     chr17   116084  116085  .       7       -       116084  116085  0,255,0 7       0
#     chr17   116353  116354  .       9       +       116353  116354  0,255,0 9       0
#     chr17   116354  116355  .       1       -       116354  116355  0,255,0 1       0
#     chr17   116445  116446  .       12      +       116445  116446  255,255,0       12      50
#     chr17   116446  116447  .       1       -       116446  116447  0,255,0 1       0
#     chr17   116703  116704  .       17      +       116703  116704  0,255,0 17      0
#     chr17   116704  116705  .       4       -       116704  116705  0,255,0 4       0
#
#     Note, that the first row above have coverage=0, so they list all CpGs (this is from WGBS data).
#     This is not the case for RRBS, where they only list covered sites:
#
#     e.g.:
#     chr1    9943211 9943212 K562_Rep3_RRBS  168     -       10003269        10003270        0,255,0 168     0
#     chr1    9943228 9943229 K562_Rep3_RRBS  168     -       10003286        10003287        0,255,0 168     1
#     chr1    9943239 9943240 K562_Rep3_RRBS  1       +       10003297        10003298        0,255,0 1       0
#     chr1    9943240 9943241 K562_Rep3_RRBS  168     -       10003298        10003299        0,255,0 168     4
#
#     ### Output files:
#     cpgDict = {"chr\tstart\tend\n" : methylation level (format: float (0-1))} # have all cpgs with specified cutoff (not only 100% 5mC or 100% 5C)
#     output coordinates are 1-based genome coordinates.
#
#     '''
#
#     cpgDict = {}
#
#     infile = open(infileName, 'r')
#
#     for row in infile:
#         tmp = row.strip().split("\t")
#         if baseCount == 1:
#             start = int(tmp[start_col])
#         elif baseCount == 0:
#             start = int(tmp[start_col]) + 1
#         else:
#             print("###\timportGroundTruth_BedMethyl_from_Encode InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#             exit
#
#         if chrFilter == False or chrFilter == tmp[chr_col]:
#             if int(tmp[cov_col]) >= covCutt:
#                 key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)
#                 if key not in cpgDict:
#                     cpgDict[key] = float(tmp[meth_col]) / 100.0
#                 else:
#                     print("###\timportGroundTruth_BedMethyl_from_Encode SanityCheckError: One CpG should not have more than 1 entry")
#
#     infile.close()
#     print("###\timportGroundTruth_BedMethyl_from_Encode: loaded information for {} CpGs".format(len(cpgDict)))
#     return cpgDict
#
#
# def e_importGroundTruth_coverage_output_from_Bismark(infileName, chr_col=0, start_col=1, meth_col=3, meth_reads_col=4, unmeth_reads_col=5, covCutt=10, baseCount=1, chrFilter=False, gzippedInput=True):
#     '''
#
#     ### Description of the columns in this format:
#
#     1. Reference chromosome or scaffold
#     2. Start position in chromosome (1-based)
#     3. End position in chromosome
#     4. methylation percentage (0-100)
#     5. methylated reads number
#     6. unmethylated reads number
#
#     I think that this output is by default 1-based. This is based on (https://www.bioinformatics.babraham.ac.uk/projects/bismark/):
#         bismark2bedGraph: This module does now produce these two output files:
#         (1) A bedGraph file, which now contains a header line: 'track type=bedGraph'. The genomic start coords are 0-based, the end coords are 1-based.
#         (2) A coverage file ending in .cov. This file replaces the former 'bedGraph --counts' file and is required to proceed with the subsequent step to generate a genome-wide cytosine report (the module doing this has been renamed to coverage2cytosine to reflect this file name change)
#     comparison of bedGraph and coverage files suggests that in the latter we deal with 1-based. Also, to get 0-based one have to activate appropriate flag in Bismark, which is not done in MINE, pipeline. I emphasize "mine" as files from somebody else might use different setting, so be careful!
#
#     e.g. structure of the coverage file:
#     chr8    205945  205945  100     2       0
#     chr8    206310  206310  100     10      0
#     chr8    206317  206317  100     10      0
#     chr8    206319  206319  90      9       1
#     chr8    206322  206322  80      8       2
#
#     ### Output files:
#     cpgDict = {"chr\tstart\tend\n" : methylation level (format: float (0-1))} # have all cpgs with specified cutoff (not only 100% 5mC or 100% 5C)
#     output coordinates are 1-based genome coordinates.
#
#     '''
#
#     cpgDict = {}
#
#     if gzippedInput:
#         infile = gzip.open(infileName, 'rb')
#     else:
#         infile = open(infileName, 'r')
#
#     for row in infile:
#         tmp = row.decode('ascii').strip().split("\t")
#         if baseCount == 1:
#             try:
#                 start = int(tmp[start_col])
#             except:
#                 logger.error(f" ### error when parse ground_truth row={row}")
#                 continue
#         elif baseCount == 0:
#             try:
#                 start = int(tmp[start_col]) + 1
#             except:
#                 logger.error(f" ### error when parse ground_truth row={row}")
#                 continue
#         else:
#             print("###\timportGroundTruth_coverage_output_from_Bismark InputValueError: baseCount value set to '{}'. It should be equal to 0 or 1".format(baseCount))
#             exit()
#
#         if chrFilter == False or chrFilter == tmp[chr_col]:
#             try:
#                 temp_meth_and_unmeth = int(tmp[meth_reads_col]) + int(tmp[unmeth_reads_col])
#             except:
#                 logger.error(f" ### Error parse gbTruth row = {row}")
#                 continue
#
#             if temp_meth_and_unmeth >= covCutt:
#                 try:
#                     key = "{}\t{}\t{}\n".format(tmp[chr_col], start, start)
#                     if key not in cpgDict:
#                         cpgDict[key] = float(tmp[meth_col]) / 100.0
#                     else:
#                         print("###\timportGroundTruth_coverage_output_from_Bismark SanityCheckError: One CpG should not have more than 1 entry")
#                 except:
#                     logger.error(f" ### Error parse gbTruth row = {row}")
#                     continue
#
#     infile.close()
#     logger.info("###\timportGroundTruth_coverage_output_from_Bismark: loaded information for {} CpGs".format(len(cpgDict)))
#     return cpgDict

#
# def plot_AUC_curve(scores, y, ax, title="", outfile=None):
#     try:
#         fpr, tpr, _ = roc_curve(y, scores)
#         roc_auc = auc(fpr, tpr)
#         lw = 2
#
#         plt.plot(fpr, tpr,
#                  lw=lw, label='{0} - ROC curve (area = {1:.4f})'.format(title, roc_auc))
#         plt.plot([0, 1], [0, 1], color='lightgrey', lw=lw, linestyle='--')
#         plt.xlim([0.0, 1.0])
#         plt.ylim([0.0, 1.05])
#         plt.xlabel('False Positive Rate')
#         plt.ylabel('True Positive Rate')
#         plt.title(title)
#         plt.legend(loc="lower right")
#     except ValueError:
#         print("###\tERROR for plot_AUC_curve: y:", y, "scores:", scores)


#
# def importGroundTruth_BS():
#     '''
#     !!! This function will be needed for NanoCompare project !!!
#     !!!!!!! its not a final version of the function, but with only some small changes it will be
#
#     note that this function was written to parse "bedMethyl" format from Encode. (description here: https://www.encodeproject.org/documents/964e2676-d0be-4b5d-aeec-f4f02310b221/@@download/attachment/WGBS%20pipeline%20overview.pdf)
#     additionally, the output file was preprocessed, such that only sites with coverage >= 10 reads and either 100% or 0% methylated, are included.
#     these will be taken to compute true positive, true negative etc.
#
#     #####
#     example usage here: http://helix122:9912/edit/li-lab/NanoporeData/WR_ONT_analyses/ai/APL_nanopolishStats/automatedSingleReadPrecission_3.py
#
#     '''
#
#     #     infile = open(infileName, 'r')
#
#     #     bsseqDict = {} # {"chr:start:end" : 100 or 0} - 100 in case if methylated, 0 if unmethylated
#     #     for row in infile:
#     #         tmp = row.strip().split("\t")
#     #         bsseqDict["{}:{}:{}".format(tmp[0], tmp[1], tmp[2])] = int(tmp[-1])
#     #     infile.close()
#     #     return bsseqDict
#     pass
#
#
# def importNarrowCpGsList():
#     pass

#
# def dict2txt(inputDict):
#     text = ""
#     for key in inputDict:
#         #         print(key)
#         text += key
#     return text
#
#
# def txt2dict(pybed):
#     d = {}
#     for t in pybed:
#         d[str(t)] = 1
#     return d

#
# def computePerReadStats(ontCalls, bsReference, title, bedFile=False, ontCutt_perRead=1, ontCutt_4corr=4, secondFilterBed=False, secondFilterBed_4Corr=False):
#     '''
#     ontCalls - dictionary of CpG coordinates with their per read methylation call (1 or 0) // Format: {"chr\tstart\tend\n" : [list of methylation calls]}
#     bsReference - dictionary of CpG coordinates with their methylation frequencies (range 0 - 1). This list is already prefiltered to meet minimal coverage (e.g. 4x) at this point. // Format: {"chr\tstart\tend\n" : methylation level (format: float (0-1))}
#     title - prefix of the analysis, output plots etc. - should be as short as possible, but unique in context of other analyses
#     bedFile - BED file which will be used to narrow down the list of CpGs for example to those inside CGIs or promoters etc.. By default "False" - which means no restrictions are done (i.e. genome wide)
#     secondFilterBed - these should be CpGs covered in some reference list. Format: BED
#
#
#     ============================================
#
#     Basically i want to fill in the table below:
#                Positive	Negative	Total
#      Presence	a	        b	        a+b
#      Absence	c	        d	        c+d
#      Total	    a+c	        b+d	        a+b+c+d
#
#     Nice summary also at wiki: https://en.wikipedia.org/wiki/F1_score
#     , where "Positive" and "Negative" corresponds with ONT based observations, while "Presence" and "Absence" is based on BS-Seq
#
#     '''
#
#     switch = 0
#     ontCalls_narrow = []
#     if bedFile != False:
#         ontCalls_bed = BedTool(dict2txt(ontCalls), from_string=True)
#         ontCalls_bed = ontCalls_bed.sort()
#         narrowBed = BedTool(bedFile)
#         narrowBed = narrowBed.sort()
#         ontCalls_intersect = ontCalls_bed.intersect(narrowBed, u=True, wa=True)
#         ontCalls_narrow = txt2dict(ontCalls_intersect)
#         suffix = bedFile
#     else:
#         switch = 1
#         suffix = "GenomeWide"
#
#     ## Second optional filter, organized in the same fashion as the first one. Designed to accomodate for example list of CpGs covered by second program
#     secondSwitch = 0
#     ontCalls_narrow_second = {}
#     if secondFilterBed != False:
#         infile = open(secondFilterBed, 'r')
#         secondFilterDict = {}
#         for row in infile:
#             secondFilterDict[row] = 0
#         infile.close()
#         ontCalls_narrow_second = dict.fromkeys(set(ontCalls.keys()).intersection(set(secondFilterDict.keys())), 0)
#     else:
#         secondSwitch = 1
#         suffix = "GenomeWide"
#
#     # Second optional filter, shoudl be used in combination with second optional filter above
#     secondSwitch_4corr = 0
#     ontCalls_narrow_second_4corr = {}
#     if secondFilterBed_4Corr != False:
#         infile = open(secondFilterBed_4Corr, 'r')
#         secondFilterDict = {}
#         for row in infile:
#             secondFilterDict[row] = 0
#         infile.close()
#         ontCalls_narrow_second_4corr = dict.fromkeys(set(ontCalls.keys()).intersection(set(secondFilterDict.keys())), 0)
#     else:
#         secondSwitch_4corr = 1
#         suffix = "GenomeWide"
#
#     suffixTMP = suffix.split("/")
#     if len(suffixTMP) > 1:
#         suffix = suffixTMP[-1]
#
#     TP_5mC = FP_5mC = FN_5mC = TN_5mC = TP_5C = FP_5C = FN_5C = TN_5C = 0
#     y = []
#     scores = []
#
#     ontSites = 0
#     mCsites = 0
#     Csites = 0
#     referenceCpGs = 0
#
#     ## four tuples for correlation:
#     ontFrequencies_4corr_mix = []  # mix(ed) are those, which in reference have methylation level >0 and <1
#     refFrequencies_4corr_mix = []
#
#     ontFrequencies_4corr_all = []  # all are all:) i.e. all CpGs with methyaltion level in refence in range 0-1
#     refFrequencies_4corr_all = []
#     leftovers = {}
#     leftovers1 = {}
#
#     for cpg_ont in ontCalls:
#         ##### for per read stats:
#         if cpg_ont in bsReference and len(ontCalls[cpg_ont]) >= ontCutt_perRead and (switch == 1 or cpg_ont in ontCalls_narrow) and (
#                 secondSwitch == 1 or cpg_ont in ontCalls_narrow_second):  # we should not take onCuttoffs for per read stats - shouldn't we? actually, we need to have the option to use this parameter, because at some point we may want to narrow down the per read stats to cover only the sites which were also covered by correlation with BS-Seq. Using this cutoff here is the easiest way to do just that
#             #         if cpg_ont in bsReference and (switch == 1 or cpg_ont in ontCalls_narrow):
#             if bsReference[cpg_ont] == 1 or bsReference[cpg_ont] == 0:  # we only consider absolute states here
#                 referenceCpGs += 1
#
#                 for ontCall in ontCalls[cpg_ont]:
#                     if bsReference[cpg_ont] == 1:
#                         mCsites += 1
#                     if bsReference[cpg_ont] == 0:
#                         Csites += 1
#                     ontSites += 1
#
#                     ### variables needed to compute precission, recall etc.:
#                     if ontCall == 1 and bsReference[cpg_ont] == 1:  # true positive
#                         TP_5mC += 1
#                     elif ontCall == 1 and bsReference[cpg_ont] == 0:  # false positive
#                         FP_5mC += 1
#                     elif ontCall == 0 and bsReference[cpg_ont] == 1:  # false negative
#                         FN_5mC += 1
#                     elif ontCall == 0 and bsReference[cpg_ont] == 0:  # true negative
#                         TN_5mC += 1
#
#                     if ontCall == 0 and bsReference[cpg_ont] == 0:  # true positive
#                         TP_5C += 1
#                     elif ontCall == 0 and bsReference[cpg_ont] == 1:  # false positive
#                         FP_5C += 1
#                     elif ontCall == 1 and bsReference[cpg_ont] == 0:  # false negative
#                         FN_5C += 1
#                     elif ontCall == 1 and bsReference[cpg_ont] == 1:  # true negative
#                         TN_5C += 1
#
#                     ### AUC related:
#                     scores.append(ontCall)
#                     y.append(bsReference[cpg_ont])
#         #             else:
#         #                 leftovers[cpg_ont] = ontCalls[cpg_ont]
#         #         else:
#         #             leftovers1[cpg_ont] = ontCalls[cpg_ont]
#
#         #         if cpg_ont in ontCalls_narrow_second:
#         #             leftovers1[cpg_ont] = ontCalls[cpg_ont]
#         #             if bsReference[cpg_ont] == 1 or bsReference[cpg_ont] == 0:
#         #                 leftovers[cpg_ont] = ontCalls[cpg_ont]
#
#         ##### for correlation stats:
#         if cpg_ont in bsReference and len(ontCalls[cpg_ont]) >= ontCutt_4corr and (switch == 1 or cpg_ont in ontCalls_narrow) and (secondSwitch_4corr == 1 or cpg_ont in ontCalls_narrow_second_4corr):
#             ontMethFreq = np.mean(ontCalls[cpg_ont])
#             ontFrequencies_4corr_all.append(ontMethFreq)
#             refFrequencies_4corr_all.append(bsReference[cpg_ont])
#             if bsReference[cpg_ont] > 0 and bsReference[cpg_ont] < 1:
#                 ontFrequencies_4corr_mix.append(ontMethFreq)
#                 refFrequencies_4corr_mix.append(bsReference[cpg_ont])
#
#     ### compute all per read stats:
#
#     #     Accuracy:
#     try:
#         accuracy = (TP_5mC + TN_5mC) / float(TP_5mC + FP_5mC + FN_5mC + TN_5mC)
#     except ZeroDivisionError:
#         accuracy = 0
#     #     print("Accuracy: {0:1.4f}".format(accuracy))
#
#     #     Positive predictive value (PPV), Precision = (TP) / E(Predicted condition positive)
#     try:
#         predicted_condition_positive_5mC = float(TP_5mC + FP_5mC)
#         precision_5mC = TP_5mC / predicted_condition_positive_5mC
#     except ZeroDivisionError:
#         precision_5mC = 0
#     #     print("Precision_5mC: {0:1.4f}".format(precision_5mC))
#     try:
#         predicted_condition_positive_5C = float(TP_5C + FP_5C)
#         precision_5C = TP_5C / predicted_condition_positive_5C
#     except ZeroDivisionError:
#         precision_5C = 0
#
#     #     print("Precision_5C: {0:1.4f}".format(precision_5C))
#
#     #     True positive rate (TPR), Recall, Sensitivity, probability of detection = (TP) / (TP+FN)
#     try:
#         recall_5mC = TP_5mC / float(TP_5mC + FN_5mC)
#     except ZeroDivisionError:
#         recall_5mC = 0
#     #     print("Recall_5mC: {0:1.4f}".format(recall_5mC))
#
#     try:
#         recall_5C = TP_5C / float(TP_5C + FN_5C)
#     except ZeroDivisionError:
#         recall_5C = 0
#     #     print("Recall_5C: {0:1.4f}".format(recall_5C))
#
#     #     F1 score:
#     try:
#         F1_5mC = 2 * ((precision_5mC * recall_5mC) / (precision_5mC + recall_5mC))
#     except ZeroDivisionError:
#         F1_5mC = 0
#     #     print("F1 score_5mC: {0:1.4f}".format(F1_5mC))
#
#     try:
#         F1_5C = 2 * ((precision_5C * recall_5C) / (precision_5C + recall_5C))
#     except ZeroDivisionError:
#         F1_5C = 0
#     #     print("F1 score_5C: {0:1.4f}".format(F1_5C))
#
#     #     print("ontSites:", ontSites)
#
#     ## plot AUC curve:
#     fig = plt.figure(figsize=(5, 5), dpi=300)
#
#     fprSwitch = 1
#     try:
#         fpr, tpr, _ = roc_curve(y, scores)
#     except ValueError:
#         print("###\tERROR for roc_curve: y:", y, "scores:", scores, "\nother settings:", title, bedFile, secondFilterBed, secondFilterBed_4Corr)
#         fprSwitch = 0
#         roc_auc = 0
#
#     if fprSwitch == 1:
#         roc_auc = auc(fpr, tpr)
#         #     print("AUC: {0:1.4f}".format(roc_auc))
#         #     print(title)
#
#         lw = 2
#
#         plt.plot(fpr, tpr, lw=lw, label='ROC curve (area = {0:.4f})'.format(roc_auc))
#         plt.plot([0, 1], [0, 1], color='lightgrey', lw=lw, linestyle='--')
#         plt.xlim([0.0, 1.0])
#         plt.ylim([0.0, 1.05])
#         plt.xlabel('False Positive Rate')
#         plt.ylabel('True Positive Rate')
#         plt.title(suffix)
#         plt.legend(loc="lower right")
#
#         fig.savefig("{}.{}.AUC.pdf".format(title, suffix), bbox_inches='tight', dpi=300)
#         plt.close()
#
#         ## Plot confusion matrix:
#         cnf_matrix = confusion_matrix(y, scores)
#         np.set_printoptions(precision=2)
#         plt.figure(figsize=(5, 5))
#         plot_confusion_matrix(cnf_matrix, classes=[0, 1], normalize=True, title=suffix)
#         plt.savefig("{}.{}.ConfusionMatrix.pdf".format(title, suffix), bbox_inches='tight', dpi=300)
#         plt.close()
#
#         ## plot Precission-recall:
#         average_precision = average_precision_score(y, scores)
#         #     print('Average precision-recall score: {0:0.2f}'.format(average_precision))
#         plt.figure(figsize=(5, 5))
#         precision, recall, _ = precision_recall_curve(y, scores)
#         # In matplotlib < 1.5, plt.fill_between does not have a 'step' argument
#         step_kwargs = ({'step': 'post'}
#                        if 'step' in signature(plt.fill_between).parameters
#                        else {})
#         plt.step(recall, precision, color='b', alpha=0.2,
#                  where='post')
#         plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)
#
#         plt.xlabel('Recall')
#         plt.ylabel('Precision')
#         plt.ylim([0.0, 1.05])
#         plt.xlim([0.0, 1.0])
#         plt.title('2-class Precision-Recall curve: AP={0:0.2f}\n{1}'.format(average_precision, suffix))
#         plt.savefig("{}.{}.PrecissionRecall.pdf".format(title, suffix), bbox_inches='tight', dpi=300)
#         plt.close()
#
#     ########################
#     # correlation based stats:
#     corrMix, pvalMix = pearsonr(ontFrequencies_4corr_mix, refFrequencies_4corr_mix)
#     corrAll, pvalAll = pearsonr(ontFrequencies_4corr_all, refFrequencies_4corr_all)
#
#     return accuracy, roc_auc, precision_5C, recall_5C, F1_5C, Csites, precision_5mC, recall_5mC, F1_5mC, mCsites, referenceCpGs, corrMix, len(ontFrequencies_4corr_mix), corrAll, len(ontFrequencies_4corr_all)  # , leftovers, leftovers1


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
#
# def plot_confusion_matrix(cm, classes,
#                           normalize=False,
#                           title="Confusion matrix",
#                           cmap=plt.cm.Blues):
#     """
#     This function prints and plots the confusion matrix.
#     Normalization can be applied by setting `normalize=True`.
#     """
#     if normalize:
#         cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
#     #         print("Normalized confusion matrix")
#     #     else:
#     #         print('Confusion matrix, without normalization')
#
#     #     print(cm)
#
#     plt.imshow(cm, interpolation='nearest', cmap=cmap)
#     plt.title(title)
#     plt.colorbar()
#     tick_marks = np.arange(len(classes))
#     plt.xticks(tick_marks, classes, rotation=45)
#     plt.yticks(tick_marks, classes)
#
#     fmt = '.2f' if normalize else 'd'
#     thresh = cm.max() / 2.
#     for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
#         plt.text(j, i, format(cm[i, j], fmt),
#                  horizontalalignment="center",
#                  color="white" if cm[i, j] > thresh else "black")
#
#     plt.ylabel('True label')
#     plt.xlabel('Predicted label')
#     plt.tight_layout()

#
# def combine2programsCalls(calls1, calls2, outfileName=False):
#     '''
#     call1 and call2 should have the following format:
#     {"chr\tstart\tend\n" : [list of methylation calls]}
#
#     result: dictionary of shared sites in the same format as input
#     '''
#
#     if outfileName != False:
#         tmp = dict.fromkeys(set(calls1.keys()).intersection(set(calls2.keys())), -1)
#         outfile = open(outfileName, 'w')
#         for key in tmp:
#             outfile.write(key)
#         outfile.close()
#         return dict.fromkeys(set(calls1.keys()).intersection(set(calls2.keys())), -1)
#     else:
#         return dict.fromkeys(set(calls1.keys()).intersection(set(calls2.keys())), -1)


# def combine2programsCalls_4Corr(calls1, calls2, cutt=4, outfileName=False):
#     '''
#     call1 and call2 should have the following format:
#     {"chr\tstart\tend\n" : [list of methylation calls]}
#
#     result: dictionary of shared sites in the same format as input with coverage of calls1 over calls 2 - use this with caution
#
#     '''
#
#     filteredCalls1 = {}
#     for call in calls1:
#         if isinstance(calls1[call], list):
#             if len(calls1[call]) >= cutt:
#                 filteredCalls1[call] = cutt  # calls1[call]
#         elif isinstance(calls1[call], int) or isinstance(calls1[call], float):
#             if calls1[call] >= cutt:
#                 filteredCalls1[call] = cutt  # calls1[call]
#     #         else:
#     #             print("WARNING ### combine2programsCalls_4Corr ### calls1[call]: {}".format(type(calls1[call])))
#
#     filteredCalls2 = {}
#     for call in calls2:
#         if isinstance(calls2[call], list):
#             if len(calls2[call]) >= cutt:
#                 filteredCalls2[call] = cutt  # calls2[call]
#         elif isinstance(calls2[call], int) or isinstance(calls2[call], float):
#             if calls2[call] >= cutt:
#                 filteredCalls2[call] = cutt  # calls2[call]
#     #         else:
#     #             print("WARNING ### combine2programsCalls_4Corr ### calls2[call]: {}".format(type(calls2[call])))
#     print(len(filteredCalls1), len(filteredCalls2))
#     if outfileName != False:
#         tmp = dict.fromkeys(set(filteredCalls1.keys()).intersection(set(filteredCalls2.keys())), cutt)
#         outfile = open(outfileName, 'w')
#         for key in tmp:
#             outfile.write(key)
#         outfile.close()
#
#         print("combine2programsCalls_4Corr DONE")
#         return dict.fromkeys(set(filteredCalls1.keys()).intersection(set(filteredCalls2.keys())), cutt)
#     else:
#         print("combine2programsCalls_4Corr DONE")
#         return dict.fromkeys(set(filteredCalls1.keys()).intersection(set(filteredCalls2.keys())), cutt)
#
#
# def combine_ONT_and_BS(ontCalls, bsReference, analysisPrefix, narrowedCoordinates=False, ontCutt=4, secondFilterBed=False, secondFilterBed_4Corr=False):
#     d = {"prefix"              : [],
#             "coord"            : [],
#             "accuracy"         : [],
#             "roc_auc"          : [],
#             "precision_5C"     : [],
#             "recall_5C"        : [],
#             "F1_5C"            : [],
#             "Csites"           : [],
#             "precision_5mC"    : [],
#             "recall_5mC"       : [],
#             "F1_5mC"           : [],
#             "mCsites"          : [],
#             "referenceCpGs"    : [],
#             "corrMix"          : [],
#             "Corr_mixedSupport": [],
#             "corrAll"          : [],
#             "Corr_allSupport"  : []}
#     for coord in narrowedCoordinates:
#         accuracy, roc_auc, precision_5C, recall_5C, F1_5C, Csites, precision_5mC, recall_5mC, F1_5mC, mCsites, referenceCpGs, corrMix, Corr_mixedSupport, corrAll, Corr_allSupport = computePerReadStats(ontCalls, bsReference, analysisPrefix, bedFile=coord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)
#         #         accuracy, roc_auc, precision_5C, recall_5C, F1_5C, Csites, precision_5mC, recall_5mC, F1_5mC, mCsites, referenceCpGs, corrMix, Corr_mixedSupport, corrAll, Corr_allSupport, leftovers, leftovers1 = computePerReadStats(ontCalls, bsReference, analysisPrefix, bedFile = coord, secondFilterBed = secondFilterBed, secondFilterBed_4Corr = secondFilterBed_4Corr)
#         d["prefix"].append(analysisPrefix)
#         d["coord"].append(coord)
#         d["accuracy"].append(accuracy)
#         d["roc_auc"].append(roc_auc)
#         d["precision_5C"].append(precision_5C)
#         d["recall_5C"].append(recall_5C)
#         d["F1_5C"].append(F1_5C)
#         d["Csites"].append(Csites)
#         d["precision_5mC"].append(precision_5mC)
#         d["recall_5mC"].append(recall_5mC)
#         d["F1_5mC"].append(F1_5mC)
#         d["mCsites"].append(mCsites)
#         d["referenceCpGs"].append(referenceCpGs)
#         d["corrMix"].append(corrMix)
#         d["Corr_mixedSupport"].append(Corr_mixedSupport)
#         d["corrAll"].append(corrAll)
#         d["Corr_allSupport"].append(Corr_allSupport)
#
#     #     df = pd.DataFrame.from_dict(d, orient='index')
#     df = pd.DataFrame.from_dict(d)
#     return df


#     return df, leftovers, leftovers1
#     accuracy, roc_auc, precision_5C, recall_5C, F1_5C, Csites, precision_5mC, recall_5mC, F1_5mC, mCsites, referenceCpGs = computePerReadStats(ontCalls, bsReference, analysisPrefix)

#     return accuracy, roc_auc, precision_5C, recall_5C, F1_5C, Csites, precision_5mC, recall_5mC, F1_5mC, mCsites, referenceCpGs
#
# def NonSingletonsScanner(referenceGenomeFile, outfileName_s, outfileName_ns):
#     '''
#     The output file is in 1-based coordinate system.
#     '''
#     reference = SeqIO.to_dict(SeqIO.parse(referenceGenomeFile, "fasta"))
#     print("###\tNonSingletonsScanner: {} reference genome parsed".format(referenceGenomeFile))
#
#     outfile_s = open(outfileName_s, "w")  # "s" stands for Singletons
#     outfile_ns = open(outfileName_ns, "w")  # "s" stands for Non-Singletons
#
#     for chromosome in list(reference.keys()):
#         idxs = re.finditer('CG', str(reference[chromosome].seq).upper())
#
#         singleton = -1  # 1 will stand for yes, 0 for no
#         for idx in idxs:
#             #             print(chromosome, idx, idx.start(), idx.end())
#             if singleton == -1:
#                 s = idx.start() + 1  # here 8: mock1 <_sre.SRE_Match object; span=(8, 10), match='CG'> 8 10
#                 e = idx.end()  # here 10: mock1 <_sre.SRE_Match object; span=(8, 10), match='CG'> 8 10
#                 singleton = 1
#             else:
#                 if (idx.start() - e) < 5:
#                     # we just found a non-singleton. I.e. accordingly to the Nanopolish approach, CGs closer than 5bp, are considered as non-singletons
#                     e = idx.end()
#                     singleton = 0
#                 else:
#                     # current CG is not part of non-singleton. It might mean that its not part of a big non-singleton or singleton upstream from it. We test which of these options below
#                     if singleton == 1:
#                         #                         print(chromosome, s, e, "SINGLETON")
#                         outfile_s.write("{}\t{}\t{}\n".format(chromosome, s, e))
#                     else:
#                         #                         print(chromosome, s, e, "NON-SINGLETON")
#                         outfile_ns.write("{}\t{}\t{}\n".format(chromosome, s, e))
#                     s = idx.start() + 1
#                     e = idx.end()
#                     singleton = 1
#
#         if singleton == 1:  # this code repetition takes care of the last instance in the long list of CG indexes
#             #             print(chromosome, s, e, "SINGLETON")
#             outfile_s.write("{}\t{}\t{}\n".format(chromosome, s, e))
#         else:
#             #             print(chromosome, s, e, "NON-SINGLETON")
#             outfile_ns.write("{}\t{}\t{}\n".format(chromosome, s, e))
#
#         print("###\tNonSingletonsScanner: chromosome {} processed".format(chromosome))
#
#     outfile_s.close()
#     outfile_ns.close()
#
#     print("###\tNonSingletonsScanner: {} file processed".format(referenceGenomeFile))
#
#
# def nonSingletonsPostprocessing(referenceMeth, regionsBedFile, dataset):
#     '''
#     This function will take the input *.bed file from "NonSingletonsScanner" funtion, which corresponds with non-singletons.
#     Next it will separate them into concordant non-singletons (i.e. fully methylated or fully unmethylated), and disconcordant (those with at least one CpG fully methylated and at least one fully unmethylated), or fully mixed (i.e. all CpGs in non-singletons have methylation level >0 and < 100)
#     This kind of preprocessing will have to be done for each studied library separately.
#     '''
#
#     refMeth = BedTool(dict2txt(referenceMeth), from_string=True)
#     refMeth = refMeth.sort()
#     regions = BedTool(regionsBedFile)
#     regions = regions.sort()
#
#     regions_refMeth = regions.intersect(refMeth, wa=True, wb=True)
#
#     regions_refMeth_dict = {}  # {regionCoords : [methylation percentage list]}
#     for ovr in regions_refMeth:
#         regionKey = "{}\t{}\t{}\n".format(ovr[0], ovr[1], ovr[2])
#         methKey = "{}\t{}\t{}\n".format(ovr[3], ovr[4], ovr[5])
#         if regionKey not in regions_refMeth_dict:
#             regions_refMeth_dict[regionKey] = []
#         regions_refMeth_dict[regionKey].append(referenceMeth[methKey])
#
#     outfile_prefix = regionsBedFile.replace(".bed", '')
#     outfile_concordant = open("{}.{}.concordant.bed".format(dataset, outfile_prefix), "w")
#     outfile_discordant = open("{}.{}.discordant.bed".format(dataset, outfile_prefix), "w")
#     outfile_fullyMixed = open("{}.{}.fullyMixed.bed".format(dataset, outfile_prefix), "w")
#     outfile_other = open("{}.{}.other.bed".format(dataset, outfile_prefix), "w")
#     for region in regions_refMeth_dict:
#         fullMeth = 0
#         nullMeth = 0
#         mixMeth = 0
#         for meth in regions_refMeth_dict[region]:
#             if meth == 1:
#                 fullMeth = 1
#             elif meth == 0:
#                 nullMeth = 1
#             else:
#                 mixMeth = 1
#
#         if (fullMeth + nullMeth) == 1 and mixMeth == 0:
#             #             print("Concordant")
#             outfile_concordant.write(region)
#         elif (fullMeth + nullMeth) == 2:
#             #             print("Discordant")
#             outfile_discordant.write(region)
#         elif (fullMeth + nullMeth) == 0 and mixMeth == 1:
#             #             print("mixed")
#             outfile_fullyMixed.write(region)
#         else:
#             #             print("What do we have here? ", fullMeth, nullMeth, mixMeth)
#             outfile_other.write("{}\t{}_{}_{}\n".format(region.strip(), fullMeth, nullMeth, mixMeth))
#     outfile_concordant.close()
#     outfile_discordant.close()
#     outfile_fullyMixed.close()
#     outfile_other.close()

#
# def singletonsPostprocessing(referenceMeth, singletonsBedFile, dataset):
#     '''
#     This function will take the input *.bed file from "NonSingletonsScanner" funtion, which corresponds with singletons.
#     Next it will separate them into absolute (i.e. fully methylated or fully unmethylated), and mixed (i.e. all CpGs in non-singletons have methylation level >0 and < 100)
#     This kind of preprocessing will have to be done for each studied library separately.
#     '''
#
#     refMeth = BedTool(dict2txt(referenceMeth), from_string=True)
#     refMeth = refMeth.sort()
#     regions = BedTool(singletonsBedFile)
#     regions = regions.sort()
#
#     refMeth_regions = refMeth.intersect(regions, wa=True, u=True)
#
#     outfile_prefix = singletonsBedFile.replace(".bed", '')
#     outfile_absolute = open("{}.{}.absolute.bed".format(dataset, outfile_prefix), "w")
#     outfile_mixed = open("{}.{}.mixed.bed".format(dataset, outfile_prefix), "w")
#     for ovr in refMeth_regions:
#         methKey = "{}\t{}\t{}\n".format(ovr[0], ovr[1], ovr[2])
#
#         if referenceMeth[methKey] == 1 or referenceMeth[methKey] == 0:
#             outfile_absolute.write(methKey)
#         else:
#             outfile_mixed.write(methKey)
#
#     outfile_absolute.close()
#     outfile_mixed.close()
#


if __name__ == '__main__':
    set_log_debug_level()

    # tool coverage cutoff 4
    minCovCutt = 4

    # bgtruth coverage cutoff 10
    bsCutt = 5

    logger.info(f'bgtruth cov={bsCutt}, tool cov={minCovCutt}')

    RunPrefix = argv[8]  # "K562_WGBS_rep_ENCFF721JMB"
    out_dir = os.path.join(pic_base_dir, RunPrefix)
    os.makedirs(out_dir, exist_ok=True)

    logger.debug(list(enumerate(argv)))
    # for argi, arg in enumerate(argv):
    #     logger.debug(f"{argi} = {arg}")

    logger.debug(f"Start load DeepSignal")
    if argv[1] == 'NO':
        DeepSignal_calls = None
    else:
        DeepSignal_calls = importPredictions_DeepSignal(argv[1])  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepSignal_runs/K562/K562_DeepSignal.MethCalls.Joined.tsv")
        DeepSignal_calls = coverageFiltering(DeepSignal_calls, minCov=minCovCutt)

    logger.debug(f"Start load Tombo")
    if argv[2] == 'NO':
        Tombo_calls = None
    else:
        Tombo_calls = importPredictions_Tombo(argv[2])  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/K562/K562_Tombo.batch_all.batch_0.perReadsStats.bed")
        Tombo_calls = coverageFiltering(Tombo_calls, minCov=minCovCutt)

    logger.debug(f"Start load Nanopolish")
    if argv[3] == 'NO':
        Nanopolish_calls = None
    else:
        # Nanopolish_calls = importPredictions_Nanopolish_v2(argv[3])  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.tsv", IncludeNonSingletons = True)
        Nanopolish_calls = importPredictions_Nanopolish(argv[3])
        Nanopolish_calls = coverageFiltering(Nanopolish_calls, minCov=minCovCutt)

    logger.debug(f"Start load DeepMod")
    if argv[4] == 'NO':
        DeepMod_calls = None
    else:
        DeepMod_calls = importPredictions_DeepMod(argv[4])  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/K562/K562.C.combined.bed")
        DeepMod_calls = coverageFiltering(DeepMod_calls, minCov=minCovCutt)

    logger.debug(f"Start load DeepMod_clusteredResultParsing")
    if argv[5] == 'NO':
        DeepMod_calls_clustered = None
    else:
        DeepMod_calls_clustered = importPredictions_DeepMod_clustered(argv[5])  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/K562/K562.C_clusterCpG.combined.bed")
        DeepMod_calls_clustered = coverageFiltering(DeepMod_calls_clustered, minCov=minCovCutt, byLength=False)

    # bgTruth = importGroundTruth_BedMethyl_from_Encode("/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/K562/ENCSR765JPC_WGBS_hg38/ENCFF721JMB.bed", covCutt = bsCutt)

    logger.debug(f"Start load bgTruth")

    # TODO all following functions removed to Universal_meth_stats_evaluation
    if argv[7] == "encode":
        bgTruth = importGroundTruth_BedMethyl_from_Encode(argv[6], covCutt=bsCutt)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/K562/ENCSR765JPC_WGBS_hg38/ENCFF721JMB.bed", chrFilter="chr20")
    elif argv[7] == "oxBS_sudo":
        bgTruth = importGroundTruth_oxBS(argv[6], covCutt=bsCutt)
    elif argv[7] == "bismark":
        bgTruth = importGroundTruth_coverage_output_from_Bismark(argv[6], covCutt=bsCutt)
    elif argv[7] == "bismark_bedgraph":
        bgTruth = importGroundTruth_coverage_output_from_Bismark_BedGraph(argv[6])
    else:
        logger.error("Methylation_correlation_plotting.py ERROR: Unknown bacground truth parser configuration. Aborting. FYI: currently supported are: encode, oxBS_sudo, bismark")
        exit(0)

    name_calls = ['DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'DeepMode_cluster']
    all_calls = [DeepSignal_calls, Tombo_calls, Nanopolish_calls, DeepMod_calls, DeepMod_calls_clustered]

    logger.debug(f"Study set intersection of each tool with bgtruth")
    dataset = []
    bgtruthCpGs = set(list(bgTruth.keys()))
    for call1, name1 in zip(all_calls, name_calls):
        if call1 is None:
            dataset.append({'bgtruth': len(bgtruthCpGs), 'tool': np.nan, 'joined': np.nan})
            continue
        overlapCpGs = bgtruthCpGs.intersection(set(list(call1.keys())))
        dataset.append({'bgtruth': len(bgtruthCpGs), 'tool': len(set(list(call1.keys()))), 'joined': len(overlapCpGs)})
        logger.info(f'BG-Truth join with {name1} get {len(overlapCpGs)} CpGs')
        outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-bgtruth-{name1}--bsCov{bsCutt}-minCov{minCovCutt}.bed')
        save_keys_to_bed(overlapCpGs, outfn)
    df = pd.DataFrame(dataset, index=name_calls)
    outfn = os.path.join(out_dir, f'{RunPrefix}-summary-bgtruth-tools-bsCov{bsCutt}-minCov{minCovCutt}.csv')
    df.to_csv(outfn)

    logger.debug(f"Start set intersection with all joined together (4+1 tools with bgtruth)")
    coveredCpGs = set(list(bgTruth.keys()))
    for call1, name1 in zip(all_calls, name_calls):
        if call1 is None:
            continue
        coveredCpGs = coveredCpGs.intersection(set(list(call1.keys())))
        logger.info(f'Join {name1} get {len(coveredCpGs)} CpGs')

    # coveredCpGs = set(list(bgTruth.keys())).intersection(set(list(DeepSignal_calls.keys()))).intersection(set(list(Tombo_calls.keys()))).intersection(set(list(Nanopolish_calls.keys()))).intersection(set(list(DeepMod_calls.keys()))).intersection(set(list(DeepMod_calls_clustered.keys())))
    logger.info(f"{len(coveredCpGs)} CpGs are covered by all tools and bgtruth")

    outfn = os.path.join(out_dir, f"Meth_corr_plot_data-{RunPrefix}-bsCov{bsCutt}-minCov{minCovCutt}-time-{current_time_str()}.tsv")
    logger.info(f"Start output results to {outfn}")

    outfile = open(outfn, 'w')
    outfile.write("chr\tstart\tend\tDeepSignal_freq\tDeepSignal_cov\tTombo_freq\tTombo_cov\tNanopolish_freq\tNanopolish_cov\tDeepMod_freq\tDeepMod_cov\tDeepMod_clust_freq\tDeepMod_clust_cov\tBSseq\n")

    for cpg in coveredCpGs:
        coords = cpg.strip().split("\t")
        outfile.write("{}\t{}\t{}\t".format(coords[0], coords[1], coords[2]))

        for call1 in all_calls:
            if call1 is None:
                outfile.write("\t\t")
                continue
            outfile.write(f"{call1[cpg][0]}\t{call1[cpg][1]}\t")

        # outfile.write("{}\t{}\t".format(DeepSignal_calls[cpg][0], DeepSignal_calls[cpg][1]))
        # outfile.write("{}\t{}\t".format(Tombo_calls[cpg][0], Tombo_calls[cpg][1]))
        # outfile.write("{}\t{}\t".format(Nanopolish_calls[cpg][0], Nanopolish_calls[cpg][1]))
        # outfile.write("{}\t{}\t".format(DeepMod_calls[cpg][0], DeepMod_calls[cpg][1]))
        # outfile.write("{}\t{}\t".format(DeepMod_calls_clustered[cpg][0], DeepMod_calls_clustered[cpg][1]))
        outfile.write(f"{bgTruth[cpg]}\n")

    outfile.close()
    logger.info(f"save to {outfn}")
