#!/usr/bin/env python3

import argparse
from multiprocessing import Pool

from sklearn.metrics import confusion_matrix

from nanocompare.eval_common import *
from nanocompare.global_settings import nonsingletonsFile
from nanocompare.global_settings import singletonFileExtStr


def evaluate_roc_pr_values(call_fn):
    """
    Sample file name is: HL60.Guppy.per_site.32.160.sorted.bed
    :param call_fn:
    :return:
    """
    m = re.search(r'Guppy\.per_site\.(\d+)\.(\d+)\.sorted\.bed', call_fn)
    if m:
        cutoff1 = int(m.group(1))
        cutoff2 = int(m.group(2))
    else:
        raise Exception(f"{call_fn} is not match with pattern")

    logger.debug(f'absoluteBGTruthCov = {len(absoluteBGTruthCov)}')
    logger.debug(f'evalSet = {len(evalSet)}')

    guppyCall = import_call(call_fn, "Guppy.ZW", include_score=False, siteLevel=False, enable_cache=enable_cache, using_cache=using_cache)
    logger.debug(f'guppyCall = {len(guppyCall)}')

    logger.debug("Start calculate precision, recall, etc.")
    y_label = []
    y_pred = []
    num_cpgs = 0
    num_methcpgs = num_unmethcpgs = 0
    num_methcalls = num_unmethcalls = 0
    for cpg in evalSet:
        label_bgtruth_cpg = 1 if is_fully_meth(absoluteBGTruthCov[cpg][0]) else 0
        if cpg not in guppyCall:
            continue
        num_cpgs += 1
        num_methcpgs += label_bgtruth_cpg
        for pred_cpg in guppyCall[cpg]:
            if pred_cpg == 1:
                num_methcalls += 1
            else:
                num_unmethcalls += 1

            y_label.append(label_bgtruth_cpg)
            y_pred.append(pred_cpg)
    num_unmethcpgs = num_cpgs - num_methcpgs

    logger.info(f'len ylabel={len(y_label)}, ylabel={sum(y_label)}, len y_pred={len(y_pred)}, ypred={sum(y_pred)}')
    ## Calculate precision, recall, tpr, fpr, etc.
    precision = precision_score(y_label, y_pred)
    recall = recall_score(y_label, y_pred)

    tn, fp, fn, tp = confusion_matrix(y_label, y_pred).ravel()

    logger.debug(f'Confusion matrix={confusion_matrix(y_label, y_pred, labels=[0, 1])}')

    tpr = tp / (tp + fn)
    fpr = fp / (tn + fp)

    ret = {"cutoff1"   : cutoff1, "cutoff2": cutoff2,
            "precision": precision, "recall": recall,
            "tpr"      : tpr, "fpr": fpr,
            "5mC.sites": num_methcpgs, "5C.sites": num_unmethcpgs,
            "5mC.calls": num_methcalls, "5C.calls": num_unmethcalls
            }
    logger.info(f'Finished for {call_fn}, ret={ret}')
    return ret


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Read level performance evaluation in Nanocompare paper')
    parser.add_argument('--min-bgtruth-cov', type=int, help="min bg-truth coverage cutoff", default=5)
    parser.add_argument('--dsname', type=str, help="dataset name", default='DS')
    parser.add_argument('--processors', type=int, help="multi-processors", default=8)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--report-joined', action='store_true', help="True if report on only joined sets")
    parser.add_argument('--test', action='store_true', help="True if only test for short time running")
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by space')
    parser.add_argument('--chrSet', nargs='+', help='chromosome list', default=humanChrSet)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name>;<file-name>", default=None)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--enable-cache', action='store_true')
    parser.add_argument('--using-cache', action='store_true')
    parser.add_argument('--joined-set-fn', type=str, help="joined sets file name", default=None)  # --calls-dir
    parser.add_argument('--calls-dir', type=str, help="guppy call dir include cutoffs", default=None)
    parser.add_argument('-mpi', action='store_true')
    parser.add_argument('--analysis', type=str, help='special analysis specifications', default="")
    return parser.parse_args()


if __name__ == '__main__':
    set_log_debug_level()
    # set_log_info_level()

    args = parse_arguments()

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
    RunPrefix = args.runid.replace('ExpPRData-', '')

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir: {out_dir}')

    # Add logging files also to result output dir
    add_logging_file(os.path.join(out_dir, 'run-results.log'))

    logger.debug(args)

    # Note: all bed files (Singleton and NonSingleton) are 1-based start, even for "chr1  123  124" (This will conform for + or - strand).
    # So we must import as 1-based format for our tool or bgtruth, DO NOT USE baseFormat=0
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
            bgTruth1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True, using_cache=using_cache, enable_cache=enable_cache)
            bgTruthList.append(bgTruth1)

        # Combine multiple bgtruth together for analysis
        # We use union of two replicates as BG-Truth
        combineBGTruth = combineBGTruthList(bgTruthList, covCutoff=1)
        logger.info("\n\n########################\n\n")

        logger.info(f'Start find absolute state (100% or 0% level), take times')
        absoluteBGTruth = {key: combineBGTruth[key] for key in combineBGTruth if satisfy_fully_meth_or_unmeth(combineBGTruth[key][0])}
        logger.info(f'Combined bgtruth sites={len(combineBGTruth):,}, Absolute bgtruth (100% and 0% level) sites={len(absoluteBGTruth):,}')

        logger.info(f'Start cutoff on absolute bg-truth, take times')
        # This is the smallest sites we use for evaluation
        absoluteBGTruthCov = {key: absoluteBGTruth[key] for key in absoluteBGTruth if absoluteBGTruth[key][1] >= cutoffBGTruth}
        logger.info(f'After apply cutoff={cutoffBGTruth}, bgtruth sites={len(absoluteBGTruthCov):,}')

        # Load all coordinate file list (full path) in this runs
        relateCoord = list(narrowCoordFileList)  # copy the basic coordinate

        ## add additional two region files based on bgtruth (Concordant, Discordant):
        ## file name is like: K562_WGBS_2Reps.hg38_nonsingletons.concordant.bed
        nonsingletonsFilePrefix = nonsingletonsFile.replace(singletonFileExtStr, '')
        fn_concordant = f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.concordant.bed"
        fn_discordant = f"{out_dir}/{RunPrefix}.{nonsingletonsFilePrefix}.discordant.bed"

        relateCoord.append(fn_concordant)
        relateCoord.append(fn_discordant)

        # Define concordant and discordant based on bg-truth (only 100% and 0% sites in BG-Truth) with cov>=1
        # Classify concordant and discordant based on cov>=1 bgtruth
        # TODO: if we use cov>= 5 BS-seq define them, please update absoluteBGTruth with absoluteBGTruthCov
        nonSingletonsPostprocessing(absoluteBGTruth, nonsingletonsFile, nsConcordantFileName=fn_concordant, nsDisCordantFileName=fn_discordant, print_first=False)

        logger.info("\n\n########################\n\n")
    else:
        absoluteBGTruth = None
        absoluteBGTruthCov = None

    joined_set_filename = args.joined_set_fn
    evalSet = load_single_sites_bed_as_set(joined_set_filename)

    fnlist = glob.glob(os.path.join(args.calls_dir, "*.Guppy.per_site.*.*.sorted.bed"))
    logger.info(fnlist)

    logger.info(f'absoluteBGTruthCov = {len(absoluteBGTruthCov)}')
    logger.info(f'evalSet = {len(evalSet)}')

    processors = args.processors
    with Pool(processes=processors) as pool:
        retList = pool.map(evaluate_roc_pr_values, fnlist[:1])
    df = pd.DataFrame(retList)

    df = df.sort_values(by=['cutoff1'])
    logger.info(df)

    outfn = os.path.join(out_dir, f'{RunPrefix}.roc.pr.data.across.cutoff.csv')
    df.to_csv(outfn, index=False)
    logger.info(f'save to {outfn}')

    logger.info("DONE")
