#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : read_level_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

"""
Read-level evaluation on methylation calls of tools, compute the performance results(F1, accuracy, ROC-AUC, etc.) in nanome paper.

This script will generate all per-read performance results.
"""

import argparse
import os.path
from multiprocessing import Manager

import pybedtools
from sklearn.metrics import confusion_matrix

from nanocompare.eval_common import *
from nanocompare.global_settings import nonsingletonsFile, save_done_file, singletonsFile, narrowCoordNameList, \
    cg_density_coord_name_list, rep_coord_name_list, nanome_version
from nanocompare.global_settings import perf_report_columns, singletonFileExtStr


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


def report_singleton_nonsingleton_table(bgTruth, outfn, fn_concordant, fn_discordant,
                                        genome_annotation_dir=os.path.join(data_base_dir, 'genome-annotation')):
    """
    Report the number of fully-methylated or unmethylated sites in Singletons, Non-Singletons, Concordant and Discordant.
    :param bgTruth:
    :param outfn:
    :param fn_concordant:
    :param fn_discordant:
    :return:
    """
    logger.debug('Start report BS seq data, the number of sites in singletons and nonsingltons')
    ret = {}

    if isinstance(bgTruth, dict):
        combineBGTruthSet = set(bgTruth.keys())
    elif isinstance(bgTruth, set):
        combineBGTruthSet = bgTruth
    else:
        raise Exception(f"Not support type for bgTruth={type(bgTruth)}")

    singletonFileName = os.path.join(genome_annotation_dir, singletonsFile)
    singletonSet = filter_cpgkeys_using_bedfile(combineBGTruthSet, singletonFileName)
    meth_unmeth = calculate_meth_unmeth(combineBGTruth, singletonSet)
    ret.update({'Singletons.5C': meth_unmeth[1], 'Singletons.5mC': meth_unmeth[0]})

    nonsingletonFilename = os.path.join(genome_annotation_dir, nonsingletonsFile)
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

    df = df[['Singletons.5C', 'Singletons.5mC', 'Singletons.sum', 'Non-Singletons.5C', 'Non-Singletons.5mC',
             'Non-Singletons.sum', 'Concordant.5C', 'Concordant.5mC', 'Concordant.sum', 'Discordant.5C',
             'Discordant.5mC', 'Discordant.sum']]

    df.to_csv(outfn)
    logger.debug(f'save to {outfn}')


def report_per_read_performance(ontCalls, bgTruth, analysisPrefix, narrowedCoordinatesList=None,
                                secondFilterBedFileName=None, cutoff_meth=1.0, outdir=None, tagname=None):
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

    for coord_tuple in tqdm(narrowedCoordinatesList):
        if coord_tuple is not None:
            if not args.large_mem:
                eval_coord_tuple = get_region_bed_tuple(coord_tuple[0],
                                                        enable_base_detection_bedfile=not args.disable_bed_check)
            else:
                eval_coord_tuple = coord_tuple
        else:
            eval_coord_tuple = None

        if eval_coord_tuple is not None and eval_coord_tuple[2] is None:
            logger.warning(
                f"Bed region tagname={eval_coord_tuple[1]} is not found, not evaluated, check genome-annotaion dir={args.genome_annotation} for file {eval_coord_tuple[0]}")
            continue
        accuracy, roc_auc, ap, f1_macro, f1_micro, \
        precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, \
        recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, \
        F1_5mC, mCalls, referenceCpGs, cSites_BGTruth, mSites_BGTruth = \
            computePerReadPerfStats(ontCalls, bgTruth, analysisPrefix, coordBedFileName=eval_coord_tuple,
                                    secondFilterBedFileName=secondFilterBedFileName,
                                    cutoff_fully_meth=cutoff_meth, outdir=outdir,
                                    tagname=tagname, save_curve_data=args.save_curve_data)

        coord_basefn = os.path.basename(f'{eval_coord_tuple[0] if eval_coord_tuple is not None else "x.x.Genome-wide"}')
        if coord_tuple is not None and not args.large_mem:
            if eval_coord_tuple is not None:
                del eval_coord_tuple  # clean memory if possible

        d["prefix"].append(analysisPrefix)
        d["coord"].append(coord_basefn)
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


# Will be deprecated
def report_per_read_performance_mp(ontCalls, bgTruth, analysisPrefix, narrowedCoordinatesList=None,
                                   secondFilterBedFileName=None, cutoff_fully_meth=1.0, outdir=None, tagname=None,
                                   processors=10):
    """
    Multi-processor version of per-read performance evaluation.

    :param ontCalls:    ONT methylation call results
    :param bgTruth:     1 or 2 replicates as bgtruth
    :param analysisPrefix:  output file tagname for this runs
    :param narrowedCoordinatesList: coordinate BED files for regions such as Genome-wide, Singletons, Promoters, etc.
    :param secondFilterBedFileName: Joined of tools and bgtruth bed file, or None for not joined
    :param cutoff_fully_meth: 1.0 means 1 is fully-meth, or 0.9 for >=0.9 is fully-meth
    :param outdir:  output dir for curve data
    :return: Dataframe of results

    results DF:
        referenceCpGs is number of all CpGs that is fully-methylated or unmethylated in BG-Truth
    """

    ret_list = []

    with Manager() as manager:
        ontCalls = manager.dict(ontCalls)
        bgTruth = manager.dict(bgTruth)

        with Pool(processes=processors) as pool:
            for coord_tuple in narrowedCoordinatesList:
                # accuracy, roc_auc, ap, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, corrMix, Corr_mixedSupport, corrAll, Corr_allSupport, cSites_BGTruth, mSites_BGTruth = \
                #     computePerReadStats(ontCalls, bgTruth, analysisPrefix, coordBedFileName=coord_fn, secondFilterBedFile=secondFilterBed,
                #                         secondFilterBed_4CorrFile=secondFilterBed_4Corr,
                #                         cutoff_meth=cutoff_meth, outdir=outdir, tagname=tagname)
                coord_basefn = os.path.basename(f'{coord_tuple[0] if coord_tuple is not None else "x.x.Genome-wide"}')
                ret = pool.apply_async(computePerReadPerfStats, (ontCalls, bgTruth, analysisPrefix,),
                                       kwds={'coordBedFileName': coord_tuple,
                                             'secondFilterBedFileName': secondFilterBedFileName,
                                             'cutoff_fully_meth': cutoff_fully_meth, 'outdir': outdir,
                                             'tagname': tagname})
                ret_list.append((coord_basefn, ret))

            pool.close()
            pool.join()

        d = defaultdict(list)
        for k in range(len(ret_list)):
            coord_basefn, ret = ret_list[k]

            # logger.info(coord)
            # logger.info(ret_list[k])
            # logger.info(ret.get())
            # accuracy, roc_auc, ap, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, corrMix, Corr_mixedSupport, corrAll, Corr_allSupport, cSites_BGTruth, mSites_BGTruth = ret.get()
            accuracy, roc_auc, ap, f1_macro, f1_micro, precision_macro, precision_micro, recall_macro, recall_micro, precision_5C, recall_5C, F1_5C, cCalls, precision_5mC, recall_5mC, F1_5mC, mCalls, referenceCpGs, cSites_BGTruth, mSites_BGTruth = ret.get()
            d["prefix"].append(analysisPrefix)
            d["coord"].append(coord_basefn)
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


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='read_level_eval (NANOME)',
                                     description='Read-level performance evaluation in nanome paper')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{nanome_version}')
    parser.add_argument('--dsname', type=str, help="dataset name", default='DS')
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by space',
                        required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name>;<file-name>",
                        default=None)
    parser.add_argument('--min-bgtruth-cov', type=int, help="min bg-truth coverage cutoff", default=5)
    parser.add_argument('--processors', type=int, help="multi-processors", default=1)
    parser.add_argument('--report-joined', action='store_true', help="true if report on only joined sets")
    parser.add_argument('--chrSet', nargs='+', help='chromosome list', default=humanChrSet)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--enable-cache', help="if enable cache functions", action='store_true')
    parser.add_argument('--using-cache', help="if use cache files", action='store_true')
    parser.add_argument('--distribution', help="report singleton/nonsingleton distributions", action='store_true')
    parser.add_argument('--mpi', help="if using multi-processing for evaluation", action='store_true')
    parser.add_argument('--analysis', type=str, help='special analysis specifications for ecoli', default="")
    parser.add_argument('--bedtools-tmp', type=str, help='bedtools temp dir', default=temp_dir)
    parser.add_argument('--genome-annotation', type=str, help='genome annotation dir',
                        default=os.path.join(data_base_dir, 'genome-annotation'))
    parser.add_argument('--reference-genome', type=str, help='genome annotation dir',
                        default=referenceGenomeFile)
    parser.add_argument('--save-curve-data', help="if save pred/truth points for curve plot", action='store_true')
    parser.add_argument('--large-mem', help="if using large memory (>100GB)", action='store_true')
    parser.add_argument('--disable-bed-check', help="if disable checking the 0/1 base format for genome annotations",
                        action='store_true')
    parser.add_argument('--debug', help="if output debug info", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    if args.debug:
        set_log_debug_level()
    else:
        set_log_info_level()

    ## Set tmp dir for bedtools
    # bedtool_tmp_dir = "/fastscratch/liuya/nanocompare/bedtools_tmp"
    os.makedirs(args.bedtools_tmp, exist_ok=True)
    pybedtools.helpers.set_tempdir(args.bedtools_tmp)

    dsname = args.dsname

    # We use coverage >= args.min_bgtruth_cov for bg-truth, but 1x coverage for ONT calls
    cutoffBGTruth = args.min_bgtruth_cov

    # Report performance on joined sites
    report_joined = args.report_joined

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

    # Note: all bed files (Singleton and NonSingleton) are 1-based start, even for "chr1  123  124" (This will conform for + or - strand).
    # So we must import as 1-based format for our tool or bgtruth, DO NOT USE baseFormat=0
    # We import and report 1-based results in our project
    baseFormat = 1

    if args.bgtruth:
        # We firstly parse and import bg-truth
        encode, fnlist = args.bgtruth.split(':')
        fnlist = fnlist.split(';')
        logger.debug(f'We are going to import BS-seq data from fnlist={fnlist}, encode={encode}')

        # Load bgtruth one/two replicates
        bgTruthList = []
        for fn in fnlist:
            # import if cov >= 1 firstly, then after join two replicates step, remove low coverage
            # bgTruth1 is dict of key->value, key=(chr, start, strand), and value=[meth.freq, cov]
            bgTruth1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True,
                                      using_cache=using_cache, enable_cache=enable_cache)
            bgTruthList.append(bgTruth1)

        # Combine multiple bgtruth together for analysis
        # We use union of two replicates as BG-Truth
        combineBGTruth = combineBGTruthList(bgTruthList, covCutoff=1)
        logger.debug("\n\n########################\n\n")

        logger.debug(f'Start find absolute state (100% or 0% level), take times')
        absoluteBGTruth = {key: combineBGTruth[key] for key in combineBGTruth if
                           satisfy_fully_meth_or_unmeth(combineBGTruth[key][0])}
        logger.debug(
            f'Combined bgtruth sites={len(combineBGTruth):,}, Absolute bgtruth (100% and 0% level) sites={len(absoluteBGTruth):,}')

        logger.debug(f'Start cutoff on absolute bg-truth, take times')
        # This is the smallest sites we use for evaluation
        absoluteBGTruthCov = {key: absoluteBGTruth[key] for key in absoluteBGTruth if
                              absoluteBGTruth[key][1] >= cutoffBGTruth}
        logger.debug(f'After apply cutoff={cutoffBGTruth}, bgtruth sites={len(absoluteBGTruthCov):,}')

        # Load all coordinate file list (full path) in this runs
        # relateCoord = list(narrowCoordFileList)  # copy the basic coordinate

        ## add additional two region files based on bgtruth (Concordant, Discordant):
        ## file name is like: K562_WGBS_2Reps.hg38_nonsingletons.concordant.bed
        ## let nonsingletonsFilePrefix = "hg38_nonsingletons"
        nonsingletonsFilePrefix = nonsingletonsFile.replace(singletonFileExtStr, '')
        fn_concordant = f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.concordant.bed.gz"
        fn_discordant = f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.discordant.bed.gz"

        # Define concordant and discordant based on bg-truth (only 100% and 0% sites in BG-Truth) with cov>=1
        # Classify concordant and discordant based on cov>=1 bgtruth
        logger.info(f"Generate concordant and discordant regions based on BS-seq")
        nonSingletonsPostprocessing(absoluteBGTruth, nonsingletonsFile, nsConcordantFileName=fn_concordant,
                                    nsDisCordantFileName=fn_discordant, genome_annotation_dir=args.genome_annotation)

        # Report singletons vs non-singletons of bgtruth with cov cutoff >= 1
        outfn = os.path.join(out_dir, f'{RunPrefix}.summary.bsseq.singleton.nonsingleton.cov1.csv')
        report_singleton_nonsingleton_table(absoluteBGTruth, outfn, fn_concordant=fn_concordant,
                                            fn_discordant=fn_discordant, genome_annotation_dir=args.genome_annotation)

        # Report singletons vs non-singletons of bgtruth with cov cutoff >= 5
        if cutoffBGTruth > 1:
            outfn = os.path.join(out_dir,
                                 f'{RunPrefix}.summary.bsseq.singleton.nonsingleton.cov{cutoffBGTruth}.table.s2.csv')
            report_singleton_nonsingleton_table(absoluteBGTruthCov, outfn, fn_concordant=fn_concordant,
                                                fn_discordant=fn_discordant,
                                                genome_annotation_dir=args.genome_annotation)

        logger.debug("\n\n########################\n\n")
        logger.info(f"Import BS-seq data done for fnlist={fnlist}")
    else:
        absoluteBGTruth = None
        absoluteBGTruthCov = None

    filterChrSet = args.chrSet

    ## Narrow down to BG-Truth if there BG-Truth is available
    ontCallWithinBGTruthDict = defaultdict()  # name->call

    ## Dataframe for cpgs in tool and each tool joined with BS-seq
    sitesDataset = defaultdict(list)

    for callstr in args.calls:
        call_encode, callfn = callstr.split(':')

        if len(callfn.strip()) == 0:  # skip empty filename
            continue

        # Only DeepMod.C/DeepMod.Cluster will always named as DeepMod
        # call_name = get_tool_name(call_encode)
        call_name = call_encode.replace('.', '_')

        # We do now allow import DeepMod.Cluster for read level evaluation
        if call_encode == 'DeepMod.Cluster':
            raise Exception(f'{call_encode} is not allowed for read level evaluation, please use DeepMod.C file here')

        ## MUST import read-level results, and include score for plot ROC curve and PR curve
        ont_call0 = import_call(callfn, call_encode, baseFormat=baseFormat, include_score=True, siteLevel=False,
                                using_cache=using_cache, enable_cache=enable_cache, filterChr=filterChrSet)
        sitesDataset['Dataset'].append(dsname)
        sitesDataset['Method'].append(call_name)
        sitesDataset['Sites'].append(len(ont_call0))
        sitesDataset['BSseq-cov5-certain'].append(len(absoluteBGTruthCov))

        if absoluteBGTruthCov:  # Filter out and keep only bg-truth cpgs, due to memory out of usage on NA19240
            logger.debug(f'Filter out CpG sites not in bgtruth for {call_name}')
            ontCallWithinBGTruthDict[call_name] = filter_cpg_dict(ont_call0,
                                                                  absoluteBGTruthCov)  # using absoluteBGTruthCov for even fewer sites
            sitesDataset[f'Join-tool-cov1-with-BSseq-cov{args.min_bgtruth_cov}-certain'].append(
                len(ontCallWithinBGTruthDict[call_name]))
            logger.debug(f'{call_name} left only sites={len(ontCallWithinBGTruthDict[call_name]):,}')
            # Clean large size object
            del ont_call0
        else:
            ontCallWithinBGTruthDict[call_name] = ont_call0

    logger.info(f"Import tools's calling done for toolist={list(ontCallWithinBGTruthDict.keys())}")
    logger.info(f"Memory report: {get_current_memory_usage()}")

    ## Report each tool (cov>=1) joined with BS-seq cov>=5 certain sites(0%, 100%)
    df = pd.DataFrame.from_dict(sitesDataset)
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
    bedfn_tool_join_bgtruth = f"{out_dir}/{RunPrefix}.Tools_BGTruth_cov{cutoffBGTruth}_Joined_baseFormat1.bed.gz"
    save_keys_to_single_site_bed(joinedCPG, outfn=bedfn_tool_join_bgtruth, callBaseFormat=baseFormat, outBaseFormat=1)

    ## Sort the bed file
    bedfn_tool_join_bgtruth_sorted = f"{out_dir}/{RunPrefix}.Tools_BGTruth_cov{cutoffBGTruth}_Joined_baseFormat1.sorted.bed.gz"
    sort_bed_file(infn=bedfn_tool_join_bgtruth, outfn=bedfn_tool_join_bgtruth_sorted)

    ## Delete not sorted file
    os.remove(bedfn_tool_join_bgtruth)

    ## Report joined CpGs in each regions, this is the really read level evaluation sites, Table S2
    outfn = os.path.join(out_dir,
                         f'{RunPrefix}.summary.bsseq.cov{cutoffBGTruth}.joined.tools.singleton.nonsingleton.table.like.s2.csv')
    report_singleton_nonsingleton_table(joinedCPG, outfn, fn_concordant=fn_concordant,
                                        fn_discordant=fn_discordant)

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

    if report_joined:  # Joined all together sites for evaluation
        perf_dir = os.path.join(out_dir, 'performance-results')
        os.makedirs(perf_dir, exist_ok=True)
        bgTruth = certainJoinedBGTruth
        secondBedFileName = bedfn_tool_join_bgtruth_sorted  # params passed for joined sets evaluation, may be remove, due to bgtruth is now joined
    else:  # only based on bgtruth joined with a tool
        perf_dir = os.path.join(out_dir, 'performance-results-nojoined')
        os.makedirs(perf_dir, exist_ok=True)
        bgTruth = certainBGTruth
        secondBedFileName = None
        raise Exception("Currently only support joined sets evaluation.")

    # Evaluated all region filename lists,
    # assume all genome annotations are in args.genome_annotation dir
    regions_full_filepath = [os.path.join(args.genome_annotation, cofn) for cofn in narrowCoordNameList[1:]] + \
                            [os.path.join(args.genome_annotation, cofn) for cofn in cg_density_coord_name_list] + \
                            [os.path.join(args.genome_annotation, cofn) for cofn in rep_coord_name_list] + \
                            [fn_concordant, fn_discordant]
    # logger.debug(f"Evaluated regions: {regions_full_filepath}")

    # Create the bed list for evaluation, save time for every loading of bed region
    # None, singelton, non-singletons, ...
    if args.large_mem:
        region_bedtuple_list = [None] + get_region_bed_pairs_list_mp(regions_full_filepath, processors=args.processors,
                                                                     enable_base_detection_bedfile=not args.disable_bed_check)
        logger.info(f"Memory report: {get_current_memory_usage()}")

    else:
        region_bedtuple_list = [None] + [(infn, map_region_fn_to_name(infn), None,)
                                         for infn in regions_full_filepath]

    if args.distribution:
        logger.debug("Report singletons/non-singletons in each genomic context regions in Fig.3 and 4")

        joined_bed = BedTool(bedfn_tool_join_bgtruth_sorted).sort()
        singleton_tuple = region_bedtuple_list[1]
        nonsingleton_tuple = region_bedtuple_list[2]
        concordant_tuple = region_bedtuple_list[-2]
        discordant_tuple = region_bedtuple_list[-1]

        logger.debug(
            f"Distribution based region bed files: {[singleton_tuple, nonsingleton_tuple, concordant_tuple, discordant_tuple]}")

        if args.large_mem:  # in memory already
            singleton_bed = singleton_tuple[2]
            nonsingleton_bed = nonsingleton_tuple[2]
            concordant_bed = concordant_tuple[2]
            discordant_bed = discordant_tuple[2]
        else:  # load on demand for limit memory
            singleton_bed = \
                get_region_bed_tuple(singleton_tuple[0], enable_base_detection_bedfile=not args.disable_bed_check)[2]
            nonsingleton_bed = \
                get_region_bed_tuple(nonsingleton_tuple[0], enable_base_detection_bedfile=not args.disable_bed_check)[2]
            concordant_bed = \
                get_region_bed_tuple(concordant_tuple[0], enable_base_detection_bedfile=not args.disable_bed_check)[2]
            discordant_bed = \
                get_region_bed_tuple(discordant_tuple[0], enable_base_detection_bedfile=not args.disable_bed_check)[2]

        datasets = defaultdict(list)
        sum_cg = 0
        for coord_tuple in region_bedtuple_list[:]:
            coordFn = None
            if coord_tuple is not None:  # get genomic region results
                if coord_tuple[1] in ['Singletons', 'Non-singletons', 'Concordant', 'Discordant']:
                    continue
                if args.large_mem:
                    (coordFn, tagname, coordBed) = coord_tuple
                else:
                    (coordFn, tagname, coordBed) = get_region_bed_tuple(coord_tuple[0],
                                                                        enable_base_detection_bedfile=not args.disable_bed_check)

                if coordBed is None:
                    logger.warning(f"genomic region {tagname} is not found, not evaluated.")
                    continue
                intersect_coord_bed = intersect_bed_regions(joined_bed, coordBed, coordFn)
            else:  # Genome-wide results
                intersect_coord_bed = joined_bed
                tagname = "Genome-wide"
            logger.debug(f"Start study tagname={tagname}, coordFn={coordFn}")

            num_total = len(intersect_coord_bed)

            if tagname.startswith('CG_'):
                sum_cg += num_total
                logger.debug(f"sanity check sum_cg={sum_cg}")

            intersect_singleton_bed = intersect_coord_bed.intersect(singleton_bed, u=True, wa=True)
            num_singleton = len(intersect_singleton_bed)

            intersect_nonsingleton_bed = intersect_coord_bed.intersect(nonsingleton_bed, u=True, wa=True)
            num_nonsingleton = len(intersect_nonsingleton_bed)

            intersect_concordant_bed = intersect_coord_bed.intersect(concordant_bed, u=True, wa=True)
            num_concordant = len(intersect_concordant_bed)

            intersect_discordant_bed = intersect_coord_bed.intersect(discordant_bed, u=True, wa=True)
            num_discordant = len(intersect_discordant_bed)

            datasets['Dataset'].append(dsname)
            datasets['Coord'].append(tagname)
            datasets['Total'].append(num_total)
            datasets['Singletons'].append(num_singleton)
            datasets['Non-singletons'].append(num_nonsingleton)
            datasets['Concordant'].append(num_concordant)
            datasets['Discordant'].append(num_discordant)

            logger.debug(
                f"Coord={tagname}, Total={num_total:,}, Total={num_total:,}, Singletons={num_singleton:,}, Non-singletons={num_nonsingleton:,}, Concordant={num_concordant:,}, Discordant={num_discordant:,}")
            # Sanity check sums
            if num_total != num_singleton + num_nonsingleton:
                logger.error(f"Found incorrect sums at {tagname}: for num_total != num_singleton + num_nonsingleton")
            if num_nonsingleton != num_concordant + num_discordant:
                logger.error(
                    f"Found incorrect sums at {tagname}: for num_nonsingleton != num_concordant+ num_discordant")

        df = pd.DataFrame.from_dict(datasets)
        outfn = os.path.join(out_dir,
                             f'{dsname}.bgtruth.certain.sites.distribution.sing.nonsing.each.genomic.cov{cutoffBGTruth}.table.s6.xlsx')
        df.to_excel(outfn)
        logger.debug(f"save to {outfn}")

    if args.mpi:  # Using mpi may cause error, not fixed, but fast running
        logger.debug('Using multi-processor function for evaluations:')

    for tool in ontCallWithinBGTruthDict:
        tmpPrefix = f'{RunPrefix}.{tool}'
        logger.info(f'Evaluating: {tmpPrefix}')

        if args.mpi:  # Using mpi may cause memory error for large data, not fixed, but fast running on small data
            # Note: narrowedCoordinatesList - all singleton (absolute and mixed) and non-singleton generated bed. ranges
            #       secondFilterBedFileName - joined sites of four tools and bg-truth. points
            df = report_per_read_performance_mp(ontCallWithinBGTruthDict[tool], bgTruth, tmpPrefix,
                                                narrowedCoordinatesList=region_bedtuple_list,
                                                secondFilterBedFileName=secondBedFileName, outdir=perf_dir,
                                                tagname=tmpPrefix, processors=args.processors)

        else:
            df = report_per_read_performance(ontCallWithinBGTruthDict[tool], bgTruth, tmpPrefix,
                                             narrowedCoordinatesList=region_bedtuple_list,
                                             secondFilterBedFileName=secondBedFileName, outdir=perf_dir,
                                             tagname=tmpPrefix)

            # This file will always report intermediate results, after for each tool, remove temp file
            tmpfn = os.path.join(perf_dir, 'performance.report.tmp.csv')
            os.remove(tmpfn)

        df['Tool'] = tool
        df['Dataset'] = dsname

        # Rename function need to be checked
        # df = rename_location_from_coordinate_name(df)
        df["Location"] = df["coord"].apply(map_region_fn_to_name)

        # Select columns to save
        df = df[perf_report_columns]

        outfn = os.path.join(perf_dir, f"{RunPrefix}.{tool.replace('.', '_')}.performance.report.csv")
        df.to_csv(outfn)
        logger.info(f"save to {outfn}")

    save_done_file(out_dir)
    logger.info(f"Memory report: {get_current_memory_usage()}")
    logger.info("### Read level performance analysis DONE.")
