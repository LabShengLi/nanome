"""
Read-level evaluation on methylation calls of tools, compute the performance results(F1, accuracy, ROC-AUC, etc.)

This script will generate all per-read performance results, bed files of singleton, non-singleton, based on results by Nanopore methylation calling tool related to BGTruth results
"""

import argparse
from multiprocessing import Manager, Pool

from nanocompare.global_settings import nonsingletonsFile
from nanocompare.global_settings import rename_location_from_coordinate_name, perf_report_columns, get_tool_name, singletonFileExtStr
from nanocompare.meth_stats.meth_stats_common import *


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Performance evaluation task')
    parser.add_argument('--min-bgtruth-cov', type=int, help="min bg-truth coverage cutoff", default=5)
    # parser.add_argument('--min-tool-cov', type=int, help="min tool coverage cutoff", default=3)
    parser.add_argument('--dsname', type=str, help="dataset name", default='DS')
    parser.add_argument('--processors', type=int, help="multi-processors", default=8)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--report-joined', action='store_true', help="True if report on only joined sets")
    parser.add_argument('--test', action='store_true', help="True if only test for short time running")
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by space', required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name>;<file-name>", required=True)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--enable-cache', action='store_true')
    parser.add_argument('--using-cache', action='store_true')
    parser.add_argument('--mpi', action='store_true')
    return parser.parse_args()


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


def report_singleton_nonsingleton_table(bgTruth, outfn, fn_concordant, fn_discordant):
    logger.info('Start report sites in singletons and nonsingltons')
    ret = {}
    combineBGTruthSet = set(bgTruth.keys())

    singletonFileName = narrowCoordFileList[1]  # Singleton file path
    singletonSet = filter_cpgkeys_using_bedfile(combineBGTruthSet, singletonFileName)

    meth_unmeth = calculate_meth_unmeth(combineBGTruth, singletonSet)
    ret.update({'Singletons.5C': meth_unmeth[1], 'Singletons.5mC': meth_unmeth[0]})

    nonsingletonFilename = narrowCoordFileList[2]  # Non-Singleton file path
    nonsingletonSet = filter_cpgkeys_using_bedfile(combineBGTruthSet, nonsingletonFilename)
    meth_unmeth = calculate_meth_unmeth(combineBGTruth, nonsingletonSet)
    ret.update({'Non-Singletons.5C': meth_unmeth[1], 'Non-Singletons.5mC': meth_unmeth[0]})

    # Concordant
    concordantFileName = fn_concordant  # f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.concordant.bed"
    concordantSet = filter_cpgkeys_using_bedfile(combineBGTruthSet, concordantFileName)
    meth_unmeth = calculate_meth_unmeth(combineBGTruth, concordantSet)
    ret.update({'Concordant.5C': meth_unmeth[1], 'Concordant.5mC': meth_unmeth[0]})

    # Discordant
    discordantFileName = fn_discordant  # f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.discordant.bed"
    discordantSet = filter_cpgkeys_using_bedfile(combineBGTruthSet, discordantFileName)
    meth_unmeth = calculate_meth_unmeth(combineBGTruth, discordantSet)
    ret.update({'Discordant.5C': meth_unmeth[1], 'Discordant.5mC': meth_unmeth[0]})

    df = pd.DataFrame([ret], index=[f'{dsname}'])

    df['Singletons.sum'] = df['Singletons.5C'] + df['Singletons.5mC']
    df['Non-Singletons.sum'] = df['Non-Singletons.5C'] + df['Non-Singletons.5mC']
    df['Concordant.sum'] = df['Concordant.5C'] + df['Concordant.5mC']
    df['Discordant.sum'] = df['Discordant.5C'] + df['Discordant.5mC']

    df = df[['Singletons.5C', 'Singletons.5mC', 'Singletons.sum', 'Non-Singletons.5C', 'Non-Singletons.5mC', 'Non-Singletons.sum', 'Concordant.5C', 'Concordant.5mC', 'Concordant.sum', 'Discordant.5C', 'Discordant.5mC', 'Discordant.sum']]

    df.to_csv(outfn)
    logger.info(f'save to {outfn}')


def report_per_read_performance(ontCalls, bgTruth, analysisPrefix, narrowedCoordinatesList=None, secondFilterBedFileName=None, cutoff_meth=1.0, outdir=None, tagname=None, test=False):
    """
    Report performance results
    :param ontCalls: tool's call
    :param bgTruth:  BS seq results as bg-truth for evaluation
    :param analysisPrefix:
    :param narrowedCoordinatesList: The bed file list for evaluation performance at regions (Genome-wide, Singleton, non-singleton, etc.)
    :param secondFilterBedFileName: None for bgTruth or Joined bed files
    :param cutoff_meth:
    :return:
    """
    d = defaultdict(list)

    for coord_fn in tqdm(narrowedCoordinatesList):
        accuracy, roc_auc, ap, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, cSites_BGTruth, mSites_BGTruth = \
            computePerReadStats(ontCalls, bgTruth, analysisPrefix, coordBedFileName=coord_fn, secondFilterBedFileName=secondFilterBedFileName,
                                cutoff_fully_meth=cutoff_meth, outdir=outdir, tagname=tagname)

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

        tmpdf = pd.DataFrame.from_dict(d)
        tmpfn = os.path.join(outdir, 'performance.report.tmp.csv')
        tmpdf.to_csv(tmpfn)

        if test:
            break
    df = pd.DataFrame.from_dict(d)
    return df


def report_per_read_performance_mp(ontCalls, bgTruth, analysisPrefix, narrowedCoordinatesList=None, secondFilterBedFileName=None, cutoff_fully_meth=1.0, outdir=None, tagname=None, processors=10):
    """
    New performance evaluation by Yang
    referenceCpGs is number of all CpGs that is fully-methylated (>=cutoff_meth) or unmethylated in BG-Truth

    :param ontCalls:
    :param bgTruth:
    :param analysisPrefix:
    :param narrowedCoordinatesList: coord such as Genome-wide, Singletons, etc.
    :param ontCutt:
    :param secondFilterBedFileName: Joined of 4 tools and bgtruth bed file, or None for not joined
    :param cutoff_fully_meth:
    :return:
    """

    ret_list = []

    with Manager() as manager:

        ontCalls = manager.dict(ontCalls)
        bgTruth = manager.dict(bgTruth)

        with Pool(processes=processors) as pool:
            for coord_fn in narrowedCoordinatesList:
                # accuracy, roc_auc, ap, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, corrMix, Corr_mixedSupport, corrAll, Corr_allSupport, cSites_BGTruth, mSites_BGTruth = \
                #     computePerReadStats(ontCalls, bgTruth, analysisPrefix, coordBedFileName=coord_fn, secondFilterBedFile=secondFilterBed,
                #                         secondFilterBed_4CorrFile=secondFilterBed_4Corr,
                #                         cutoff_meth=cutoff_meth, outdir=outdir, tagname=tagname)
                coord = os.path.basename(f'{coord_fn if coord_fn else "x.x.Genome-wide"}')
                ret = pool.apply_async(computePerReadStats, (ontCalls, bgTruth, analysisPrefix,), kwds={'coordBedFileName': coord_fn, 'secondFilterBedFileName': secondFilterBedFileName, 'cutoff_fully_meth': cutoff_fully_meth, 'outdir': outdir, 'tagname': tagname})
                ret_list.append((coord, ret))

            pool.close()
            pool.join()

        d = defaultdict(list)
        for k in range(len(ret_list)):
            coord, ret = ret_list[k]

            # logger.info(coord)
            # logger.info(ret_list[k])
            # logger.info(ret.get())
            # accuracy, roc_auc, ap, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, corrMix, Corr_mixedSupport, corrAll, Corr_allSupport, cSites_BGTruth, mSites_BGTruth = ret.get()
            accuracy, roc_auc, ap, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, cSites_BGTruth, mSites_BGTruth = ret.get()
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

    df = pd.DataFrame.from_dict(d)
    return df


if __name__ == '__main__':
    set_log_debug_level()

    args = parse_arguments()

    dsname = args.dsname
    cutoffBGTruth = args.min_bgtruth_cov
    # tool_cov_cutoff = args.min_tool_cov # No use of tool coverage due to per-read level evaluation
    report_joined = args.report_joined
    enable_cache = args.enable_cache
    using_cache = args.using_cache

    RunPrefix = args.runid.replace('MethPerf-', '')

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir: {out_dir}')

    # Add logging files also to result output dir
    add_logging_file(os.path.join(out_dir, 'run-results.log'))

    logger.debug(args)

    # Note: all bed files (Singleton and NonSingleton) are 1-based start, even for "chr1  123  124" (This will conform for + or - strand).
    # So we must import as 1-based format for our tool or bgtruth, DO NOT USE baseFormat=0
    baseFormat = 1

    # We firstly import bg-truth, then each tool results, and remove non-bg-truth cpgs for memory usage
    encode, fnlist = args.bgtruth.split(':')
    fnlist = fnlist.split(';')
    logger.debug(f'BGTruth fnlist={fnlist}, encode={encode}')

    bgTruthList = []
    for fn in fnlist:
        # import if cov >= 1 firstly, then after join two replicates step, remove low coverage
        bgTruth1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True, using_cache=using_cache, enable_cache=enable_cache)
        bgTruthList.append(bgTruth1)

    # Combine multiple bgtruth together for analysis
    # We use joined of two replicates as BG-Truth
    combineBGTruth = combineBGTruthList(bgTruthList, covCutoff=1)

    # Define concordant and discordant based on bg-truth (only 100% and 0% sites in BG-Truth)
    logger.info(f'Start find absolute state (100% or 0% level), take times')
    absoluteBGTruth = {key: combineBGTruth[key] for key in combineBGTruth if satisfy_fully_meth_or_unmeth(combineBGTruth[key][0])}
    logger.info(f'Combined bgtruth sites={len(combineBGTruth):,}, Absolute bgtruth (100% and 0% level) sites={len(absoluteBGTruth):,}')

    logger.info(f'Start cutoff on absolute bg-truth, take times')
    absoluteBGTruthCov = {key: absoluteBGTruth[key] for key in absoluteBGTruth if absoluteBGTruth[key][1] >= cutoffBGTruth}
    logger.info(f'After apply cutoff={cutoffBGTruth}, bgtruth sites={len(absoluteBGTruthCov):,}')

    callfn_dict = defaultdict()  # callname -> filename
    ontCallWithinBGTruthDict = defaultdict()  # name->call
    callname_list = []  # [DeepSignal, DeepMod, etc.]

    for callstr in args.calls:
        call_encode, callfn = callstr.split(':')

        call_name = get_tool_name(call_encode)
        callname_list.append(call_name)
        callfn_dict[call_name] = callfn

        ## MUST import read-level results, and include score for plot ROC curve and PR curve
        call0 = import_call(callfn, call_encode, baseFormat=baseFormat, include_score=True, deepmod_cluster_freq_cov_format=False, using_cache=using_cache, enable_cache=enable_cache)
        # Filter out and keep only bg-truth cpgs, due to memory out of usage on NA19240
        logger.info('Filter out CpG sites not in bgtruth')
        ontCallWithinBGTruthDict[call_name] = filter_cpg_dict(call0, absoluteBGTruth)
        logger.info(f'Left only sites={len(ontCallWithinBGTruthDict[call_name]):,}')

    # logger.debug(callfn_dict)

    # The all coordinate file list (full path) in this runs
    relateCoord = list(narrowCoordFileList)  # copy the basic coordinate

    ## add additional two region files based on bgtruth (Concordant, Discordant):
    nonsingletonsFilePrefix = nonsingletonsFile.replace(singletonFileExtStr, '')
    fn_concordant = f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.concordant.bed"
    fn_discordant = f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.discordant.bed"

    relateCoord.append(fn_concordant)
    relateCoord.append(fn_discordant)

    # Current results version
    # nonSingletonsPostprocessing(absoluteBGTruth, nonsingletonsFile, nsConcordantFileName=fn_concordant, nsDisCordantFileName=fn_discordant, print_first=True)

    nonSingletonsPostprocessing_same_deepsignal(absoluteBGTruth, nonsingletonsFile, nsConcordantFileName=fn_concordant, nsDisCordantFileName=fn_discordant, print_first=True)

    # Report singletons vs non-singletons of bgtruth with cov cutoff >= 1
    outfn = os.path.join(out_dir, f'{RunPrefix}.summary.singleton.nonsingleton.cov1.csv')
    report_singleton_nonsingleton_table(absoluteBGTruth, outfn, fn_concordant=fn_concordant, fn_discordant=fn_discordant)

    # Report singletons vs non-singletons of bgtruth with cov cutoff >= 5
    outfn = os.path.join(out_dir, f'{RunPrefix}.summary.singleton.nonsingleton.cov{cutoffBGTruth}.csv')
    report_singleton_nonsingleton_table(absoluteBGTruthCov, outfn, fn_concordant=fn_concordant, fn_discordant=fn_discordant)

    logger.info("\n\n########################\n\n")

    # this file is the all tool joined together sites
    bedfn_tool_join_bgtruth = f"{out_dir}/{RunPrefix}.Tools_BGTruth_cov{cutoffBGTruth}_Joined.bed"

    joinedCPG = set(absoluteBGTruthCov.keys())
    for toolname in ontCallWithinBGTruthDict:
        joinedCPG = joinedCPG.intersection(set(ontCallWithinBGTruthDict[toolname].keys()))

    save_keys_to_single_site_bed(joinedCPG, outfn=bedfn_tool_join_bgtruth, callBaseFormat=baseFormat, outBaseFormat=1)

    logger.info(f"Data points for joined all tools with bg-truth (cov>={cutoffBGTruth}) sites={len(joinedCPG):,}\n\n")

    # Next calculate fully methylated and unmethylated sites
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

    logger.info(f'The fully meth or unmeth sites of Joined BGTruth = {len(certainJoinedBGTruth):,}  (5C={cnt5C:,}, 5mC={cnt5mC:,}) for performance comparison')

    certainBGTruth = dict(absoluteBGTruthCov)  # not used now
    cntNoJoined5C = 0
    cntNoJoined5mC = 0
    for key in absoluteBGTruthCov:
        if is_fully_meth(absoluteBGTruthCov[key][0]):
            cntNoJoined5mC += 1
        else:
            cntNoJoined5C += 1

    logger.info(f'The fully meth or unmeth sites of No-Joined BGTruth (cov>={cutoffBGTruth}) = {len(certainBGTruth):,}  (5C={cntNoJoined5C:,}, 5mC={cntNoJoined5mC:,}) for performance comparison if No-joined')

    logger.info("\n\n############\n\n")

    if report_joined:  # Joined together evaluation
        perf_dir = os.path.join(out_dir, 'performance-results')
        os.makedirs(perf_dir, exist_ok=True)
        bgTruth = certainJoinedBGTruth
        secondBedFileName = bedfn_tool_join_bgtruth
    else:  # only based on bgtruth
        perf_dir = os.path.join(out_dir, 'performance-results-nojoined')
        os.makedirs(perf_dir, exist_ok=True)
        bgTruth = certainBGTruth
        secondBedFileName = None

    for tool in ontCallWithinBGTruthDict:
        tmpPrefix = f'{RunPrefix}.{tool}'
        logger.info(f'Evaluating: {tmpPrefix}')

        # Note: narrowedCoordinatesList - all singleton (absolute and mixed) and non-singleton generated bed. ranges
        #       secondFilterBedFileName - joined sites of four tools and bg-truth. points

        if args.mpi:  # Using mpi may cause error, not fixed, but fast running
            df = report_per_read_performance_mp(ontCallWithinBGTruthDict[tool], bgTruth, tmpPrefix, narrowedCoordinatesList=relateCoord, secondFilterBedFileName=secondBedFileName, outdir=perf_dir, tagname=tmpPrefix, processors=args.processors)
        else:
            df = report_per_read_performance(ontCallWithinBGTruthDict[tool], bgTruth, tmpPrefix, narrowedCoordinatesList=relateCoord, secondFilterBedFileName=secondBedFileName, outdir=perf_dir, tagname=tmpPrefix)

        df['Tool'] = tool
        df['Dataset'] = dsname

        # Rename function need to be checked
        df = rename_location_from_coordinate_name(df)

        # logger.debug(df)
        # logger.debug(df.columns)

        # Select columns to save
        df = df[perf_report_columns]

        outfn = os.path.join(perf_dir, f"{RunPrefix}.{tool}.performance.report.csv")
        df.to_csv(outfn)
        logger.info(f"save to {outfn}")

        # This file will always report intermediate results
        # tmpfn = os.path.join(perf_dir, 'performance.report.tmp.csv')
        # os.remove(tmpfn)

        if args.test:
            break
    logger.info("Meth stats performance evaluation results generation DONE.")
