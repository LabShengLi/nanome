#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : read_level_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Read-level evaluation on methylation calls of tools, compute the performance results(F1, accuracy, ROC-AUC, etc.) in nanome paper.

This script will generate all per-read performance results.
"""

import argparse
import hashlib
import os.path
from concurrent.futures import as_completed

import pybedtools
from sklearn.metrics import confusion_matrix

from nanome.common.eval_common import *
from nanome.common.global_settings import NANOME_VERSION, perf_report_columns, \
    save_done_file, \
    region_filename_dict, region_tagname_dict, sing_tagname, nonsing_tagname, concord_tagname, discord_tagname, \
    load_genome_annotation_config


def calculate_meth_unmeth(bgTruth, keySet):
    """
    bgTruth is format of key->value, key=cpg, value=[freq, cov]
    :param bgTruth:
    :param keySet:
    :return:
    """
    num5c = num5mc = 0
    for key in keySet:
        if is_fully_meth(bgTruth[key][0]):
            num5mc += 1
        elif is_fully_unmeth(bgTruth[key][0]):
            num5c += 1
    return num5mc, num5c


def report_singleton_nonsingleton_bsseq_table(bgTruth, outfn, fn_concordant, fn_discordant,
                                              genome_annotation):
    """
    Report the number of fully-methylated or unmethylated sites of BS-seq as class of Singletons, Non-Singletons, Concordant and Discordant.
    :param bgTruth:
    :param outfn:
    :param fn_concordant:
    :param fn_discordant:
    :return:
    """
    logger.debug('Start report BS seq data, the number of sites in singletons and nonsingltons')
    ret = {}

    if isinstance(bgTruth, dict):
        bgTruthSet = set(bgTruth.keys())
    elif isinstance(bgTruth, set):
        bgTruthSet = bgTruth
    else:
        raise Exception(f"Not support type for bgTruth={type(bgTruth)}")

    singletonFileName = os.path.join(genome_annotation, singletonsFile)
    singletonSet = filter_cpgkeys_using_bedfile(bgTruthSet, singletonFileName)
    meth_unmeth = calculate_meth_unmeth(absoluteBGTruth, singletonSet)
    ret.update({'Singletons.5C': meth_unmeth[1], 'Singletons.5mC': meth_unmeth[0]})

    nonsingletonFilename = os.path.join(genome_annotation, nonsingletonsFile)
    nonsingletonSet = filter_cpgkeys_using_bedfile(bgTruthSet, nonsingletonFilename)
    meth_unmeth = calculate_meth_unmeth(absoluteBGTruth, nonsingletonSet)
    ret.update({'Non-Singletons.5C': meth_unmeth[1], 'Non-Singletons.5mC': meth_unmeth[0]})

    # Concordant
    concordantFileName = fn_concordant  # f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.concordant.bed"
    concordantSet = filter_cpgkeys_using_bedfile(bgTruthSet, concordantFileName)
    meth_unmeth = calculate_meth_unmeth(absoluteBGTruth, concordantSet)
    ret.update({'Concordant.5C': meth_unmeth[1], 'Concordant.5mC': meth_unmeth[0]})

    # Discordant
    discordantFileName = fn_discordant  # f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.discordant.bed"
    discordantSet = filter_cpgkeys_using_bedfile(bgTruthSet, discordantFileName)
    meth_unmeth = calculate_meth_unmeth(absoluteBGTruth, discordantSet)
    ret.update({'Discordant.5C': meth_unmeth[1], 'Discordant.5mC': meth_unmeth[0]})

    df = pd.DataFrame([ret], index=[f'{dsname}'])

    df['Singletons.sum'] = df['Singletons.5C'] + df['Singletons.5mC']
    df['Non-Singletons.sum'] = df['Non-Singletons.5C'] + df['Non-Singletons.5mC']
    df['Concordant.sum'] = df['Concordant.5C'] + df['Concordant.5mC']
    df['Discordant.sum'] = df['Discordant.5C'] + df['Discordant.5mC']

    df = df[['Singletons.5C', 'Singletons.5mC', 'Singletons.sum', 'Non-Singletons.5C', 'Non-Singletons.5mC',
             'Non-Singletons.sum', 'Concordant.5C', 'Concordant.5mC', 'Concordant.sum', 'Discordant.5C',
             'Discordant.5mC', 'Discordant.sum']]

    df.to_csv(outfn)
    logger.debug(f'save to {outfn}')


def report_per_read_performance(ontCalls, bgTruth, analysisPrefix, narrowedCoordinatesList=None,
                                secondFilterBedFileName=None, cutoff_meth=1.0, outdir=None, prefix_name=None):
    """
    Report performance results
    :param ontCalls: tool's call
    :param bgTruth:  BS seq results as bg-truth for evaluation
    :param analysisPrefix:
    :param narrowedCoordinatesList: The bed file list for evaluation performance at regions (Genome-wide, Singleton, non-singleton, etc.), it is a tuple (basefn, tagname, bed_of_region)
    :param secondFilterBedFileName: None for bgTruth or Joined bed files
    :param cutoff_meth:
    :return:
    """
    d = defaultdict(list)
    bar = tqdm(narrowedCoordinatesList)
    for coord_tuple in bar:
        bar.set_description(
            f"Read-level-{prefix_name}-{genome_wide_tagname if coord_tuple is None else coord_tuple[1]}")
        if coord_tuple[1] != genome_wide_tagname:
            if not args.large_mem:
                eval_coord_tuple = \
                    get_region_bed_tuple(coord_tuple[0],
                                         enable_base_detection_bedfile=not args.disable_bed_check,
                                         enable_cache=args.enable_cache,
                                         using_cache=args.using_cache,
                                         cache_dir=ds_cache_dir)
            else:
                eval_coord_tuple = coord_tuple
        else:  # genome-wide setting
            eval_coord_tuple = coord_tuple

        if eval_coord_tuple[1] != genome_wide_tagname and eval_coord_tuple[2] is None:
            logger.debug(
                f"Bed region tagname={eval_coord_tuple[1]} is not found, not evaluated, check genome-annotaion dir={args.genome_annotation} for file {eval_coord_tuple[0]}")
            continue
        accuracy, roc_auc, ap, f1_macro, f1_micro, \
        precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, \
        recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, \
        F1_5mC, mCalls, referenceCpGs, cSites_BGTruth, mSites_BGTruth, tagname = \
            computePerReadPerfStats(ontCalls, bgTruth, analysisPrefix, coordBedFileName=eval_coord_tuple,
                                    secondFilterBedFileName=secondFilterBedFileName,
                                    cutoff_fully_meth=cutoff_meth, outdir=outdir,
                                    prefix_name=prefix_name, save_curve_data=args.save_curve_data)

        if tagname is None:
            continue

        d["prefix"].append(analysisPrefix)
        d["coord"].append(tagname)
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

        ## save a temp file for each region, for preview of results
        tmpdf = pd.DataFrame.from_dict(d)
        tmpfn = os.path.join(outdir, 'performance.report.tmp.csv')
        tmpdf.to_csv(tmpfn)
    df = pd.DataFrame.from_dict(d)
    logger.info(f"Memory report: {get_current_memory_usage()}")
    return df


def report_per_read_performance_mpi(ontCalls, bgTruth, analysisPrefix, narrowedCoordinatesList=None,
                                    secondFilterBedFileName=None, cutoff_meth=1.0, outdir=None, prefix_name=None):
    """
    Report performance results, multi-thread version
    :param ontCalls: tool's call
    :param bgTruth:  BS seq results as bg-truth for evaluation
    :param analysisPrefix:
    :param narrowedCoordinatesList: The bed file list for evaluation performance at regions (Genome-wide, Singleton, non-singleton, etc.), it is a tuple (basefn, tagname, bed_of_region)
    :param secondFilterBedFileName: None for bgTruth or Joined bed files
    :param cutoff_meth:
    :return:
    """
    executor = ThreadPoolExecutor(max_workers=args.processors)
    all_tasks = []

    global progress_bar_global_read
    progress_bar_global_read = tqdm(total=len(narrowedCoordinatesList))
    progress_bar_global_read.set_description(f"MT-Read-level-{analysisPrefix}")

    for coord_tuple in narrowedCoordinatesList:
        ## Load bed
        if coord_tuple[1] != genome_wide_tagname:
            if not args.large_mem:
                eval_coord_tuple = \
                    get_region_bed_tuple(coord_tuple[0],
                                         enable_base_detection_bedfile=not args.disable_bed_check,
                                         enable_cache=args.enable_cache,
                                         using_cache=args.using_cache,
                                         cache_dir=ds_cache_dir)
            else:  # already in memory
                eval_coord_tuple = coord_tuple
        else:  # genome-wide, (None, 'Genome-wide', None)
            eval_coord_tuple = coord_tuple
        ## skip for None bed regions, except for genome-wide
        if eval_coord_tuple[1] != genome_wide_tagname and eval_coord_tuple[2] is None:
            logger.debug(
                f"Bed region tagname={eval_coord_tuple[1]} is not found, not evaluated, check genome-annotaion dir={args.genome_annotation} for file {eval_coord_tuple[0]}")
            continue
        ## compute read-level performance for a bed tuple of a region
        future = executor.submit(computePerReadPerfStats, ontCalls, bgTruth, analysisPrefix,
                                 coordBedFileName=eval_coord_tuple,
                                 secondFilterBedFileName=secondFilterBedFileName,
                                 cutoff_fully_meth=cutoff_meth, outdir=outdir,
                                 prefix_name=prefix_name, save_curve_data=args.save_curve_data)
        future.add_done_callback(update_progress_bar_read_level)
        all_tasks.append(future)
    executor.shutdown()
    progress_bar_global_read.close()

    datasets = []
    for future in all_tasks:
        ret = future.result()
        ret_dict = unpack_read_level_perf_ret_dict(ret)
        if ret_dict is None:  # skip None tagname results
            continue
        ret_dict.update({"prefix": analysisPrefix})
        datasets.append(ret_dict)
    df = pd.DataFrame(datasets)
    logger.info(f"Memory report: {get_current_memory_usage()}")
    return df


def report_ecoli_metro_paper_evaluations(ontCallDict, evalCPGSet, threshold=0.2):
    """
    We now simply check results performance on positive data for E. coli.
    :param ontCallDict:
    :return:
    """
    per_read_dataset = defaultdict(list)
    per_base_dataset = defaultdict(list)

    for callname in ontCallDict:
        call = ontCallDict[callname]

        nsites = 0
        numcalls = 0
        methcalls = 0
        unmethcalls = 0
        if evalCPGSet:
            cpgSet = evalCPGSet
        else:
            cpgSet = set(call.keys())
        ## Read level evaluation
        for cpg in cpgSet:
            meth_indicator_list = [tt[0] for tt in call[cpg]]
            nsites += 1
            methcalls += sum(meth_indicator_list)
            unmethcalls += len(meth_indicator_list) - sum(meth_indicator_list)
            numcalls += len(meth_indicator_list)
        per_read_dataset['Method'].append(callname)
        per_read_dataset['#Base'].append(nsites)
        per_read_dataset['#Call'].append(numcalls)
        per_read_dataset['#Pos'].append(methcalls)
        per_read_dataset['#Neg'].append(unmethcalls)
        per_read_dataset['#Pos/#Call'].append(methcalls / numcalls)

        ## Base level evaluation
        ylabel = []
        ypred = []
        for cpg in cpgSet:
            meth_indicator_list = [tt[0] for tt in call[cpg]]
            meth_percentage = sum(meth_indicator_list) / len(meth_indicator_list)
            ylabel.append(1)
            ypred.append(1 if meth_percentage >= threshold else 0)

        npmatrix = confusion_matrix(ylabel, ypred)

        precision = precision_score(ylabel, ypred)
        recall = recall_score(ylabel, ypred)

        per_base_dataset['Method'].append(callname)
        per_base_dataset['#Base'].append(nsites)
        per_base_dataset[f'Methylated base >= {threshold:.2f}'].append(sum(ypred))
        per_base_dataset[f'Unmethylated base'].append(len(ypred) - sum(ypred))
        per_base_dataset['Precision'].append(precision)
        per_base_dataset['Recall'].append(recall)

    df = pd.DataFrame.from_dict(per_read_dataset)
    logger.info(df)

    if evalCPGSet:
        tag = "joined_sets"
    else:
        tag = "no_joined-sets"
    outfn = os.path.join(out_dir, f'report_ecoli_metro_paper_evaluations_on_{tag}.read.level.xlsx')
    df.to_excel(outfn)

    df = pd.DataFrame.from_dict(per_base_dataset)
    logger.info(df)

    if evalCPGSet:
        tag = "joined_sets"
    else:
        tag = "no_joined-sets"
    outfn = os.path.join(out_dir, f'report_ecoli_metro_paper_evaluations_on_{tag}.base.level.xlsx')
    df.to_excel(outfn)

    pass


def import_ont_calls_for_read_level(toolname, call_encode, callfn, score_cutoff, absoluteBGTruthCov,
                                    multi_processor=False):
    """
    Import read level of ont calls
    Args:
        call_encode:
        callfn:
        absoluteBGTruthCov:
        multi_processor:

    Returns:

    """
    ## MUST import read-level results, and include score for plot ROC curve and PR curve
    ## ont_call0 is raw ont-calls, too large, it will be cut to only with bs-seq, named ont_call1
    ont_call0 = import_call(callfn, call_encode, baseFormat=baseFormat, include_score=True, siteLevel=False,
                            filterChr=args.chrSet, using_cache=using_cache, enable_cache=enable_cache,
                            cache_dir=ds_cache_dir, toolname=toolname, score_cutoff=score_cutoff)
    sites_summary = {'Dataset': dsname,
                     'Method': toolname,
                     'Sites': len(ont_call0),
                     f'BSseq-cov{args.min_bgtruth_cov}-certain': len(
                         absoluteBGTruthCov) if absoluteBGTruthCov else None,
                     }

    if absoluteBGTruthCov:  # Filter out and keep only bg-truth cpgs, due to memory out of usage on NA19240
        logger.debug(f'Filter out CpG sites not in bgtruth for {toolname}')
        ont_call1 = filter_cpg_dict(ont_call0,
                                    absoluteBGTruthCov,
                                    toolname=toolname)  # using absoluteBGTruthCov for even fewer sites
        del ont_call0
        sites_summary.update({f'Join-tool-cov1-with-BSseq-cov{args.min_bgtruth_cov}-certain': len(ont_call1)})
        logger.debug(f'{toolname}:{call_encode} left only sites={len(ont_call1):,}')
    else:
        ont_call1 = ont_call0
    ont_call2 = filter_cpg_dict_by_cov(ont_call1, coverage=args.toolcov_cutoff)

    if not multi_processor:  # sequencial/multithreading running, in same process, can directly return object
        return (toolname, ont_call2, sites_summary,)

    raise Exception("Under development below")
    # Save to temp, for main process use at multi-processing
    outfn = f"tmp_ont_calls_for_read_level_{args.runid}_{os.path.basename(callfn)}.pkl"
    outfnmd5 = os.path.join(args.bedtools_tmp, "tmp_mp_ont_" + hashlib.md5(outfn.encode('utf-8')).hexdigest() + ".pkl")
    logger.debug(f"[MP message] '{outfn}' is encoded into '{outfnmd5}', transfer to main process")
    with open(outfnmd5, 'wb') as handle:
        pickle.dump(ont_call1, handle)
    del ont_call1
    logger.debug(f"Memory report: {get_current_memory_usage()}")
    return (toolname, outfnmd5, sites_summary,)


def import_bsseq_for_read_level(infn, encode, multi_processor=False):
    """
    Import bs-seq data
    Args:
        infn:
        encode:

    Returns:

    """
    # import if cov >= 1 firstly, then after join two replicates step, remove low coverage
    # bgTruth1 is dict of key->value, key=(chr, start, strand), and value=[meth.freq, cov]
    bg1 = import_bgtruth(infn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True,
                         using_cache=using_cache, enable_cache=enable_cache, cache_dir=ds_cache_dir)
    if not multi_processor:  # in same process, directly return
        return bg1
    # Save to temp, for main process use, multi-processing
    outfn = f"tmp_bsseq_for_read_level_{args.runid}_{os.path.basename(infn)}.pkl"
    outfnmd5 = os.path.join(args.bedtools_tmp,
                            "tmp_mp_bsseq_" + hashlib.md5(outfn.encode('utf-8')).hexdigest() + ".pkl")
    logger.debug(f"[MP message] '{outfn}' is encoded into '{outfnmd5}', transfer to main process")
    with open(outfnmd5, 'wb') as handle:
        pickle.dump(bg1, handle)

    del bg1
    return (infn, outfnmd5,)


def update_progress_bar_read_level(*a):
    """
    Update progress for multiprocessing
    :param a:
    :return:
    """
    global progress_bar_global_read
    progress_bar_global_read.update()


def compute_dist_at_region_mp(joined_bed, four_region_bed_list, region_tuple):
    logging.debug(f"Distribution analysis for Region={region_tuple}")
    singleton_bed, nonsingleton_bed, concordant_bed, discordant_bed = four_region_bed_list

    coordFn = None  # will be used for intersection later, None for genome-wide
    if region_tuple[1] != genome_wide_tagname:
        if region_tuple[1] in [sing_tagname, nonsing_tagname, concord_tagname, discord_tagname]:
            return None
        if args.large_mem:
            (coordFn, tagname, coordBed) = region_tuple
        else:
            (coordFn, tagname, coordBed) = \
                get_region_bed_tuple(region_tuple[0],
                                     enable_base_detection_bedfile=not args.disable_bed_check,
                                     enable_cache=args.enable_cache,
                                     using_cache=args.using_cache,
                                     cache_dir=ds_cache_dir)
        if coordBed is None:
            logger.debug(f"genomic region {tagname} is not found, not evaluated.")
            return None
        intersect_coord_bed = intersect_bed_regions(joined_bed, coordBed, coordFn)
    else:  # Genome-wide results, keep using joined
        intersect_coord_bed = joined_bed
        tagname = region_tuple[1]

    if tagname is None:
        logger.debug(f"ERROR: No region registered in config file for region_tuple={region_tuple}")
        return None

    logger.debug(f"Start study tagname={tagname}, coordFn={coordFn}")
    num_total = len(intersect_coord_bed)
    intersect_singleton_bed = intersect_coord_bed.intersect(singleton_bed, u=True, wa=True)
    num_singleton = len(intersect_singleton_bed)
    del intersect_singleton_bed

    intersect_nonsingleton_bed = intersect_coord_bed.intersect(nonsingleton_bed, u=True, wa=True)
    num_nonsingleton = len(intersect_nonsingleton_bed)
    del intersect_nonsingleton_bed

    intersect_concordant_bed = intersect_coord_bed.intersect(concordant_bed, u=True, wa=True)
    num_concordant = len(intersect_concordant_bed)
    del intersect_concordant_bed

    intersect_discordant_bed = intersect_coord_bed.intersect(discordant_bed, u=True, wa=True)
    num_discordant = len(intersect_discordant_bed)
    del intersect_discordant_bed
    del intersect_coord_bed

    ret = {
        'Dataset': dsname,
        'Coord': tagname,
        'Total': num_total,
        sing_tagname: num_singleton,
        nonsing_tagname: num_nonsingleton,
        concord_tagname: num_concordant,
        discord_tagname: num_discordant
    }

    logger.debug(
        f"Coord={tagname}, Total={num_total:,}, Singletons={num_singleton:,}, Non-singletons={num_nonsingleton:,}, Total(Sing+Nonsing)={num_singleton + num_nonsingleton:,}, Concordant={num_concordant:,}, Discordant={num_discordant:,}, Total(Con+Disc)={num_concordant + num_discordant:,}")
    # Sanity check sums
    if num_total != num_singleton + num_nonsingleton:
        logger.debug(f"WARN: Found incorrect sums at {tagname}: for num_total != num_singleton + num_nonsingleton")
    if num_nonsingleton != num_concordant + num_discordant:
        logger.debug(
            f"WARN: Found incorrect sums at {tagname}: for num_nonsingleton != num_concordant+ num_discordant")
    return ret


def unpack_read_level_perf_ret_dict(ret_list):
    """
    parse all outputs of read-level performance function
    Args:
        ret_list:

    Returns:

    """
    accuracy, roc_auc, ap, f1_macro, f1_micro, \
    precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, \
    recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, \
    F1_5mC, mCalls, referenceCpGs, cSites_BGTruth, mSites_BGTruth, tagname = ret_list

    if tagname is None:
        return None

    ret_dict = {
        "coord": tagname,
        "Accuracy": accuracy,
        "Average-Precision": ap,
        "Macro-F1": f1_macro,
        "Micro-F1": f1_micro,
        "Macro-Precision": precision_macro,
        "Micro-Precision": precision_micro,
        "Macro-Recall": recall_macro,
        "Micro-Recall": recall_micro,
        "ROC-AUC": roc_auc,
        "Precision_5C": precision_5C,
        "Recall_5C": recall_5C,
        "F1_5C": F1_5C,
        "Csites_called": cCalls,
        "Csites": cSites_BGTruth,
        "Precision_5mC": precision_5mC,
        "Recall_5mC": recall_5mC,
        "F1_5mC": F1_5mC,
        "mCsites_called": mCalls,
        "mCsites": mSites_BGTruth,
        "referenceCpGs": referenceCpGs
    }
    return ret_dict


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='read_level_eval (NANOME)',
                                     description='Read-level performance evaluation in nanome paper')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str,
                        help="running prefix/output folder name, such as MethPerf-Dataset_WGBS_2Reps",
                        required=True)
    parser.add_argument('--calls', nargs='+',
                        help='all ONT call results <tool-name>:<file-encode>:<file-name> seperated by space, tool-name/file-encode can be Nanopolish, Megalodon, DeepSignal, Guppy, Tombo, METEORE, DeepMod, NANOME',
                        required=True)
    parser.add_argument('--bgtruth', type=str,
                        help="background truth file <encode-type>:<file-name1>;<file-name2>, encode-type can be 'encode' or 'bismark'",
                        default=None)
    parser.add_argument('--genome-annotation', type=str,
                        help='genome annotation dir, contain BED files',
                        default=None)
    parser.add_argument('--min-bgtruth-cov', type=int, help="min bg-truth coverage cutoff, default is 5", default=5)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff for coverage in nanopore tools, default is >=1",
                        default=1)
    parser.add_argument('--processors', type=int, help="number of processors used, default is 1", default=1)
    parser.add_argument('--report-no-join', action='store_true', help="true if report not on joined sets")
    parser.add_argument('--chrSet', nargs='+', help='chromosome list, default is human chr1-22, X and Y',
                        default=HUMAN_CHR_SET)
    parser.add_argument('-o', type=str, help=f"output base dir, default is {pic_base_dir}", default=pic_base_dir)
    parser.add_argument('--enable-cache', help="if enable cache functions", action='store_true')
    parser.add_argument('--using-cache', help="if use cache files", action='store_true')
    parser.add_argument('--distribution', help="if report singleton/nonsingleton distributions at all regions",
                        action='store_true')
    parser.add_argument('--bsseq-report', help="if report singleton/nonsingleton in bs-seq", action='store_true')
    parser.add_argument('--analysis', type=str, help='special analysis specifications for ecoli', default="")
    parser.add_argument('--save-curve-data', help="if save pred/truth points for curve plot", action='store_true')
    parser.add_argument('--large-mem', help="if using large memory (>100GB) for speed up", action='store_true')
    parser.add_argument('--bedtools-tmp', type=str, help=f'bedtools temp dir, default is {global_temp_dir}',
                        default=global_temp_dir)
    parser.add_argument('--cache-dir', type=str,
                        help=f'cache dir used for loading calls/bs-seq (speed up running), default is {global_cache_dir}',
                        default=global_cache_dir)
    parser.add_argument('--disable-bed-check',
                        help="if disable auto-checking the 0/1 base format for genome annotations",
                        action='store_true')
    parser.add_argument('--mpi',
                        help="if using multi-processing/threading for evaluation, it can speed-up but may need more memory",
                        action='store_true')
    parser.add_argument('--mpi-import',
                        help="if using multi-processing/threading for import, it can speed-up, only for small size data",
                        action='store_true')
    parser.add_argument('--config', help="if print out config file for genome annotation", action='store_true')
    parser.add_argument('--verbose', help="if output verbose info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    if args.verbose:
        set_log_debug_level()
    else:
        set_log_info_level()

    dsname = args.dsname

    ## Set tmp dir for bedtools, each process use a bed tmp dir
    ## because the tmp dir files may be cleaned by the end of the process
    bed_temp_dir = os.path.join(args.bedtools_tmp, f"{dsname}_perf")
    os.makedirs(bed_temp_dir, exist_ok=True)
    pybedtools.helpers.set_tempdir(bed_temp_dir)

    ## Set cache dir for each dataset
    if args.enable_cache or args.using_cache:
        ds_cache_dir = os.path.join(args.cache_dir, dsname)
        # os.makedirs(ds_cache_dir, exist_ok=True)
    else:
        ds_cache_dir = None

    # We use coverage >= args.min_bgtruth_cov for bg-truth, but 1x coverage for ONT calls
    cutoffBGTruth = args.min_bgtruth_cov

    # If enable cache, loaded results will be saved to cache
    enable_cache = args.enable_cache

    # If enable and using, import functions can read from cache
    using_cache = args.using_cache

    # runid is always like 'MethPerf-K562_WGBS_2Reps', remove first word as RunPrefix like K562_WGBS_2Reps
    RunPrefix = args.runid.replace('MethPerf-', '')

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir: {out_dir}')

    # Add logging files also to result output dir
    add_logging_file(os.path.join(out_dir, 'run-results.log'))
    logger.debug(args)

    if args.config:
        load_genome_annotation_config(verbose=True)

    if sing_tagname in region_tagname_dict and nonsing_tagname in region_tagname_dict:
        singletonsFile = region_tagname_dict[sing_tagname][0]
        nonsingletonsFile = region_tagname_dict[nonsing_tagname][0]
        logger.debug(f"singletonsFile={singletonsFile}, nonsingletonsFile={nonsingletonsFile}")
    else:
        singletonsFile = nonsingletonsFile = "None-file"

    ## Test if singleton/nonsinglton BED file exists
    exists_singleton_or_nonsingleton = True
    if args.genome_annotation is None or \
            not os.path.exists(os.path.join(args.genome_annotation, singletonsFile)) or \
            not os.path.exists(os.path.join(args.genome_annotation, nonsingletonsFile)):
        exists_singleton_or_nonsingleton = False
    logger.debug(
        f"Detection of singleton/nonsingleton files: exists_singleton_or_nonsingleton={exists_singleton_or_nonsingleton}")
    # Note: all bed files (Singleton and NonSingleton) are 1-based start, even for "chr1  123  124" (This will conform for + or - strand).
    # So we must import as 1-based format for our tool or bgtruth, DO NOT USE baseFormat=0
    # We import and report 1-based results in our project
    baseFormat = 1

    if args.bgtruth:
        # We firstly parse and import bg-truth
        encode, fnlist = args.bgtruth.split(':')
        fnlist = fnlist.split(';')
        logger.debug(f'We are going to import BS-seq data from fnlist={fnlist}, encode={encode}')

        # Load bgtruth one/two replicates, using multiprocessing
        bgTruthList = []
        arg_list = []
        for fn in fnlist:
            arg_list.append((fn, encode,))
        if args.mpi:
            executor = ThreadPoolExecutor(max_workers=args.processors)
            all_task = [executor.submit(import_bsseq_for_read_level, *arg) for arg in arg_list]
            executor.shutdown()
            for future in as_completed(all_task):
                bg1 = future.result()
                bgTruthList.append(bg1)
            # with Pool(args.processors) as pool:
            #     # ret_list is list of (infn, outfn)
            #     ret_list = pool.starmap(import_bsseq_for_read_level, arg_list)
            #     for ret in ret_list:
            #         # read output from a sub-process
            #         with open(ret[1], 'rb') as handle:
            #             logger.debug(f"Load bs-seq encode={ret[0]} from file={ret[1]}, it is from subprocess")
            #             bg1 = pickle.load(handle)
            #         bgTruthList.append(bg1)
            logger.debug("MT import for bs-seq finished")
        else:
            for arg in arg_list:
                bg1 = import_bsseq_for_read_level(*arg)
                bgTruthList.append(bg1)
            logger.debug("SEQ import for bs-seq finished")

        logger.debug(f"bgTruthList={len(bgTruthList)}")
        if len(bgTruthList) == 0:
            logger.info(f"Found no BS-seq data provided, no need to run evalutation")
            sys.exit(0)

        # Combine multiple bgtruth together for analysis
        # We use union of two replicates as BG-Truth
        combineBGTruth = combineBGTruthList(bgTruthList, covCutoff=1)

        # Clean bgTruthList
        del bgTruthList

        logger.debug("\n\n########################\n\n")

        logger.debug(f'Start find absolute state (100% or 0% level), take times')
        absoluteBGTruth = {key: combineBGTruth[key] for key in combineBGTruth if
                           satisfy_fully_meth_or_unmeth(combineBGTruth[key][0])}
        logger.debug(
            f'Combined bgtruth sites={len(combineBGTruth):,}, Absolute bgtruth (100% and 0% level) sites={len(absoluteBGTruth):,}')

        ## Clean all level combineBGTruth
        del combineBGTruth

        logger.debug(f'Start cutoff on absolute bg-truth, take times')
        # This is the smallest sites we use for evaluation
        absoluteBGTruthCov = {key: absoluteBGTruth[key] for key in absoluteBGTruth if
                              absoluteBGTruth[key][1] >= cutoffBGTruth}

        logger.debug(f'After apply cutoff={cutoffBGTruth}, bgtruth sites={len(absoluteBGTruthCov):,}')

        ## add additional two region files based on bgtruth (Concordant, Discordant):
        ## file name is like: K562_WGBS_2Reps.hg38_nonsingletons.concordant.bed
        ## let nonsingletonsFilePrefix = "hg38_nonsingletons"
        nonsingletonsFilePrefix = nonsingletonsFile.replace('.bed.gz',
                                                            '')  # TODO: check all conco/disco file name and tag mapping

        # concordant and discordant file are dataset dependent
        # file like: HL60_RRBS_2Reps_HL60.hg38_nonsingletons_10bp.concordant.bed.gz
        fn_concordant = f"{out_dir}/{RunPrefix}_{args.dsname}.{nonsingletonsFilePrefix}.concordant.bed.gz"
        fn_discordant = f"{out_dir}/{RunPrefix}_{args.dsname}.{nonsingletonsFilePrefix}.discordant.bed.gz"

        if exists_singleton_or_nonsingleton:
            # Define concordant and discordant based on bg-truth (only 100% and 0% sites in BG-Truth) with cov>=1
            # Classify concordant and discordant based on cov>=1 bgtruth
            logger.info(f"Scan concordant and discordant regions based on BS-seq")
            nonSingletonsPostprocessing(absoluteBGTruth, nonsingletonsFile, nsConcordantFileName=fn_concordant,
                                        nsDisCordantFileName=fn_discordant,
                                        genome_annotation_dir=args.genome_annotation)

        if args.bsseq_report and exists_singleton_or_nonsingleton:
            # Report singletons vs non-singletons of bgtruth with cov cutoff >= 1
            outfn = os.path.join(out_dir, f'{RunPrefix}.summary.bsseq.singleton.nonsingleton.cov1.csv')
            logger.debug(f"For coverage >= 1")
            report_singleton_nonsingleton_bsseq_table(absoluteBGTruth, outfn, fn_concordant=fn_concordant,
                                                      fn_discordant=fn_discordant,
                                                      genome_annotation=args.genome_annotation)

            # Report singletons vs non-singletons of bgtruth with cov cutoff >= 5
            if cutoffBGTruth > 1:
                outfn = os.path.join(out_dir,
                                     f'{RunPrefix}.summary.bsseq.singleton.nonsingleton.cov{cutoffBGTruth}.table.s2.csv')
                logger.debug(f"For coverage >= {cutoffBGTruth}")
                report_singleton_nonsingleton_bsseq_table(absoluteBGTruthCov, outfn, fn_concordant=fn_concordant,
                                                          fn_discordant=fn_discordant,
                                                          genome_annotation=args.genome_annotation)
        logger.debug("\n\n########################\n\n")
        logger.debug(f"Memory report: {get_current_memory_usage()}")
        logger.info(f"Import BS-seq data done for fnlist={fnlist}")
    else:
        absoluteBGTruth = None
        absoluteBGTruthCov = None
        logger.info("WARN: Can not get BS-seq")

    ## Narrow down to BG-Truth if there BG-Truth is available
    ontCallWithinBGTruthDict = defaultdict()  # name->call

    arg_list = []
    for callstr in args.calls:
        try:
            if len(callstr.split(':')) == 3:
                toolname, callencode, callfn = callstr.split(':')
                score_cutoff = None
            elif len(callstr.split(':')) == 5:
                toolname, callencode, callfn, cutoff1, cutoff2 = callstr.split(':')
                cutoff1 = float(cutoff1)
                cutoff2 = float(cutoff1)
                score_cutoff = (cutoff1, cutoff2)
        except:
            raise Exception(f"--calls params is not correct: {callstr}")

        if len(callfn.strip()) == 0:  # skip empty filename
            continue
        arg_list.append((toolname, callencode, callfn, score_cutoff, absoluteBGTruthCov, False,))

    ## Dataframe for cpgs in tool and each tool joined with BS-seq
    sites_database_list = []  # list of dict
    ontCallWithinBGTruthDict = dict()  # callname-> call-dict

    if args.mpi_import:
        # Multi-thread may lead to large memory (out-of-memory) for big ont-calls
        executor = ThreadPoolExecutor(max_workers=args.processors)
        all_task = [executor.submit(import_ont_calls_for_read_level, *arg) for arg in arg_list]
        executor.shutdown()
        for future in all_task:
            ret = future.result()
            ontCallWithinBGTruthDict[ret[0]] = ret[1]
            sites_database_list.append(ret[2])
        logger.debug("MT parallel import ont calls finished")
        logger.debug(f"Memory report: {get_current_memory_usage()}")
    else:
        # Sequencial import will be benefit for large data
        for arg in arg_list:
            ret = import_ont_calls_for_read_level(*arg)
            logger.debug(f"Memory report: {get_current_memory_usage()}")
            ontCallWithinBGTruthDict[ret[0]] = ret[1]
            sites_database_list.append(ret[2])

    logger.info(f"Import tools's calling done for toolist={list(ontCallWithinBGTruthDict.keys())}")
    logger.debug(f"sites_database_list={sites_database_list}")
    logger.info(f"Memory report: {get_current_memory_usage()}")

    ## Report each tool (cov>=1) joined with BS-seq cov>=5 certain sites(0%, 100%)
    df = pd.DataFrame(sites_database_list)
    outfn = os.path.join(out_dir,
                         f'{dsname}.tools.cov1.join.with.bsseq.cov{args.min_bgtruth_cov}.read.level.report.xlsx')
    df.to_excel(outfn)

    logger.debug("\n\n########################\n\n")

    # Study the joined CpG sites by all tools (cov>=1) with BG-Truth (cov>=5, 0% or 100% sites),
    # evaluation on joined CpG by default
    if absoluteBGTruthCov:
        joinedCPG = set(absoluteBGTruthCov.keys())
    else:
        joinedCPG = None

    for toolname in ontCallWithinBGTruthDict:
        if joinedCPG is None:
            joinedCPG = set(ontCallWithinBGTruthDict[toolname].keys())
            continue
        joinedCPG = joinedCPG.intersection(set(ontCallWithinBGTruthDict[toolname].keys()))
        logger.debug(f'After joined with {toolname}, cpgs={len(joinedCPG):,}')

    # this file is the all tool (cov>=1) joined with BS-seq (cov>=5) 0%, 100% together sites BED file, for evaluation on joined sites
    bedfn_tool_join_bgtruth = f"{out_dir}/{RunPrefix}_{args.dsname}.Tools_Certain_BGTruth_cov{cutoffBGTruth}_Joined_baseFormat1.bed.gz"
    save_keys_to_single_site_bed(joinedCPG, outfn=bedfn_tool_join_bgtruth, callBaseFormat=baseFormat, outBaseFormat=1)

    ## Sort the bed file
    bedfn_tool_join_bgtruth_sorted = f"{out_dir}/{RunPrefix}_{args.dsname}.Tools_Certain_BGTruth_cov{cutoffBGTruth}_Joined_baseFormat1.sorted.bed.gz"
    sort_bed_file(infn=bedfn_tool_join_bgtruth, outfn=bedfn_tool_join_bgtruth_sorted)

    ## Delete not sorted file
    os.remove(bedfn_tool_join_bgtruth)

    if args.bsseq_report and exists_singleton_or_nonsingleton:
        ## Report joined CpGs in each regions, this is the really read level evaluation sites, Table S2
        outfn = os.path.join(out_dir,
                             f'{RunPrefix}.summary.bsseq.cov{cutoffBGTruth}.joined.tools.singleton.nonsingleton.table.like.s2.csv')
        logger.debug("Start report joined CPGs of BS-seq for singleton/nonsingleton")
        report_singleton_nonsingleton_bsseq_table(joinedCPG, outfn, fn_concordant=fn_concordant,
                                                  fn_discordant=fn_discordant, genome_annotation=args.genome_annotation)

    ## Note all tools using cov>=1 for evaluation read-leval performance
    logger.debug(
        f"Data points for joined all tools with bg-truth (if any, cov>={cutoffBGTruth}) sites={len(joinedCPG):,}\n\n")

    if "ecoli_metropaper_sanity" in args.analysis:
        report_ecoli_metro_paper_evaluations(ontCallWithinBGTruthDict, joinedCPG)
        logger.debug(f'Analysis:[{args.analysis}] DONE')

        report_ecoli_metro_paper_evaluations(ontCallWithinBGTruthDict, None)
        logger.debug(f'Analysis:[{args.analysis}] DONE')

        sys.exit(0)

    # Next extract sites in joined set only
    certainJoinedBGTruth = {}  # all joined bg-truth
    cnt5C = 0
    cnt5mC = 0
    for key in joinedCPG:
        if satisfy_fully_meth_or_unmeth(absoluteBGTruthCov[key][0]):
            # only add joined CpG keys
            certainJoinedBGTruth[key] = absoluteBGTruthCov[key]
            if is_fully_meth(absoluteBGTruthCov[key][0]):
                cnt5mC += 1
            else:
                cnt5C += 1

    logger.info(
        f'Joined BGTruth (cov>={cutoffBGTruth}) = {len(certainJoinedBGTruth):,}  (5C={cnt5C:,}, 5mC={cnt5mC:,}) for performance comparison on fully meth or unmeth sites with tools')

    certainBGTruth = dict(absoluteBGTruthCov)  # not used now
    cntNoJoined5C = 0
    cntNoJoined5mC = 0
    for key in certainBGTruth:
        if is_fully_meth(certainBGTruth[key][0]):
            cntNoJoined5mC += 1
        else:
            cntNoJoined5C += 1

    logger.info(
        f'No-Joined BGTruth (cov>={cutoffBGTruth}) = {len(certainBGTruth):,}  (5C={cntNoJoined5C:,}, 5mC={cntNoJoined5mC:,}) for performance comparison on fully meth or unmeth sites with tools')

    logger.debug("\n\n############\n\n")

    if args.report_no_join:
        # only based on bgtruth joined with a tool
        perf_dir = os.path.join(out_dir, 'performance-results-nojoined')
        os.makedirs(perf_dir, exist_ok=True)
        eval_bgTruth = certainBGTruth  # evaluation bs-seq
        secondBedFileName = None
        raise Exception("Under development")
    else:
        # Joined all together sites for evaluation
        perf_dir = os.path.join(out_dir, 'performance-results')
        os.makedirs(perf_dir, exist_ok=True)
        eval_bgTruth = certainJoinedBGTruth  # evaluation bs-seq
        # params passed for joined sets evaluation, may be remove, due to bgtruth is now joined
        secondBedFileName = bedfn_tool_join_bgtruth_sorted

    # Evaluated all region filename lists,
    # assume all genome annotations are in args.genome_annotation dir
    annot_dir = args.genome_annotation if args.genome_annotation is not None else '.'

    # region file path from genome-wide, singletons, to genic/intergenic, cg-density, and repetitive, and then concordant and discordant
    regions_full_filepath = [None] + [os.path.join(annot_dir, cofn) for cofn in region_filename_dict.keys()] + \
                            [fn_concordant, fn_discordant]

    # Create the bed list for evaluation, save time for every loading of bed region
    # Genome-wide, singelton, non-singletons, ...
    if args.large_mem:
        eval_region_tuple_list = get_region_bed_pairs_list_mp(
            regions_full_filepath, processors=args.processors,
            enable_base_detection_bedfile=not args.disable_bed_check,
            enable_cache=args.enable_cache,
            using_cache=args.using_cache,
            cache_dir=ds_cache_dir)
        logger.info(f"Memory report: {get_current_memory_usage()}")
    else:
        eval_region_tuple_list = [(infn, get_region_tagname(infn), None,)
                                  for infn in regions_full_filepath]

    if args.distribution and exists_singleton_or_nonsingleton:
        logger.debug("Report singletons/non-singletons in each genomic context regions in Fig.3 and 4")

        joined_bed = BedTool(bedfn_tool_join_bgtruth_sorted).sort()

        singleton_tuple = eval_region_tuple_list[1]
        nonsingleton_tuple = eval_region_tuple_list[2]
        concordant_tuple = eval_region_tuple_list[-2]
        discordant_tuple = eval_region_tuple_list[-1]

        ## Preload singleton, non-singleton, concordant and discordant for all threads usage
        if args.large_mem:  # in memory already
            singleton_bed = singleton_tuple[2]
            nonsingleton_bed = nonsingleton_tuple[2]
            concordant_bed = concordant_tuple[2]
            discordant_bed = discordant_tuple[2]
        else:  # load on demand for limit memory
            singleton_bed = \
                get_region_bed_tuple(singleton_tuple[0],
                                     enable_base_detection_bedfile=not args.disable_bed_check,
                                     enable_cache=args.enable_cache,
                                     using_cache=args.using_cache,
                                     cache_dir=ds_cache_dir)[2]
            nonsingleton_bed = \
                get_region_bed_tuple(nonsingleton_tuple[0],
                                     enable_base_detection_bedfile=not args.disable_bed_check,
                                     enable_cache=args.enable_cache,
                                     using_cache=args.using_cache,
                                     cache_dir=ds_cache_dir)[2]
            concordant_bed = \
                get_region_bed_tuple(concordant_tuple[0],
                                     enable_base_detection_bedfile=not args.disable_bed_check,
                                     enable_cache=args.enable_cache,
                                     using_cache=args.using_cache,
                                     cache_dir=ds_cache_dir)[2]
            discordant_bed = \
                get_region_bed_tuple(discordant_tuple[0],
                                     enable_base_detection_bedfile=not args.disable_bed_check,
                                     enable_cache=args.enable_cache,
                                     using_cache=args.using_cache,
                                     cache_dir=ds_cache_dir)[2]
        four_region_bed_list = (singleton_bed, nonsingleton_bed, concordant_bed, discordant_bed)
        ## Check not None for single/nonsingl, concord and discord
        for bedk in four_region_bed_list:
            if bedk is None:
                raise Exception(f"Not enough bed data, dist_region_tuple_list={four_region_bed_list}")

        ret_list = []

        executor = ThreadPoolExecutor(max_workers=args.processors)

        global progress_bar_global_read
        progress_bar_global_read = tqdm(total=len(eval_region_tuple_list))
        progress_bar_global_read.set_description("MT-Distribution(singl/nonsingle)")
        all_task = []
        for reg_tuple in eval_region_tuple_list:
            future = executor.submit(compute_dist_at_region_mp, joined_bed, four_region_bed_list, reg_tuple)
            future.add_done_callback(update_progress_bar_read_level)
            all_task.append(future)
        executor.shutdown()
        progress_bar_global_read.close()

        datasets = []  # list of dict for dataframe
        for future in all_task:
            ret = future.result()
            if ret is None:
                continue
            datasets.append(ret)

        df = pd.DataFrame(datasets)
        outfn = os.path.join(out_dir,
                             f'{dsname}.bgtruth.certain.sites.distribution.sing.nonsing.each.genomic.cov{cutoffBGTruth}.table.s6.xlsx')
        df.to_excel(outfn)
        logger.debug(f"save to {outfn}")

    if args.mpi:  # Using mpi may cause error, not fixed, but fast running
        logger.debug('Using multi-threading function for evaluations')

    for tool in ontCallWithinBGTruthDict:
        tmpPrefix = f'{RunPrefix}.{tool}'
        logger.info(f'Evaluating per-read performance: {tmpPrefix}')

        if args.mpi:
            df = report_per_read_performance_mpi(ontCallWithinBGTruthDict[tool], eval_bgTruth, tmpPrefix,
                                                 narrowedCoordinatesList=eval_region_tuple_list,
                                                 secondFilterBedFileName=secondBedFileName, outdir=perf_dir,
                                                 prefix_name=tmpPrefix)
        else:
            df = report_per_read_performance(ontCallWithinBGTruthDict[tool], eval_bgTruth, tmpPrefix,
                                             narrowedCoordinatesList=eval_region_tuple_list,
                                             secondFilterBedFileName=secondBedFileName, outdir=perf_dir,
                                             prefix_name=tmpPrefix)

            # This file will always report intermediate results, after for each tool, remove temp file
            tmpfn = os.path.join(perf_dir, 'performance.report.tmp.csv')
            if os.path.exists(tmpfn):
                os.remove(tmpfn)

        df['Tool'] = tool
        df['Dataset'] = dsname

        # Rename function need to be checked
        df["Location"] = df["coord"]

        # Select columns to save
        df = df[perf_report_columns]

        outfn = os.path.join(perf_dir, f"{RunPrefix}.{tool.replace('.', '_')}.performance.report.csv")
        df.to_csv(outfn)
        logger.info(f"save to {outfn}")

    save_done_file(out_dir)
    logger.info(f"Memory report: {get_current_memory_usage()}")
    logger.info("### Read level performance analysis DONE.")
