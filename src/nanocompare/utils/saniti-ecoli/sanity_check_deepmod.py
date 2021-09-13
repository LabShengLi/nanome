#!/usr/bin/env python3

"""
Sanity check NA12878 for DeepMod.
"""
import argparse
from multiprocessing import Pool

from nanocompare.eval_common import *
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score

def report_performance_deepmod(ontCall, bgTruth, threshold=0.5):
    """
    Join with bgtruth, then evaluate
    :param ontCall:
    :param bgTruth:
    :return:
    """
    evalSet = set(ontCall.keys()).intersection(set(bgTruth.keys()))

    y_truth = []
    y_pred = []
    for cpg in evalSet:
        y_truth.append(bgTruth[cpg][0])
        y_pred.append(1 if ontCall[cpg][0] > threshold else 0)
    methCnt = sum(y_truth)
    unmethCnt = len(y_truth) - methCnt
    prec = precision_score(y_truth, y_pred)
    recall = recall_score(y_truth, y_pred)

    f1 = f1_score(y_truth, y_pred)
    accuracy = accuracy_score(y_truth, y_pred)
    ret = {"Precision": prec, "Recall": recall,
            "F1"      : f1, "Accuracy": accuracy, "Threshold": threshold,
            "Unmeth"  : unmethCnt, "Meth": methCnt}

    return ret


def report_meth_unmeth_table_by_chr(chr):
    logger.debug(f'Start for chr={chr}')
    chrSet = [chr]
    # Combine one/two replicates using DeepMod methods
    bgTruth1, ret1 = combineBGTruthList_by_DeepModPaper(bgTruthList, covCutoff=1, filterChrs=chrSet)

    bgTruth5 = {key: bgTruth1[key] for key in bgTruth1 if bgTruth1[key][1] >= 5}
    methCnt5 = sum([bgTruth5[key][0] for key in bgTruth5])
    ret5 = (methCnt5, len(bgTruth5) - methCnt5)

    bgTruth10 = {key: bgTruth1[key] for key in bgTruth1 if bgTruth1[key][1] >= 10}
    methCnt10 = sum([bgTruth10[key][0] for key in bgTruth10])
    ret10 = (methCnt10, len(bgTruth10) - methCnt10)

    ret = (chr, ret1[2], ret1[1], ret1[0], ret5[1], ret5[0], ret10[1], ret10[0])
    logger.debug(f'Finished for chr={chr}, ret={ret}')
    return ret


def report_meth_unmeth_table():
    chrSet = list(humanChrSet)
    chrSet.remove('chrY')
    chrSet.remove('chr22')
    logger.debug(f'chrSet={chrSet}')

    with Pool(args.processors) as pool:
        retList = pool.map(report_meth_unmeth_table_by_chr, chrSet)

    dataset = defaultdict(list)
    for ret in retList:
        dataset['chr'].append(ret[0])
        dataset['Not-used.cov1'].append(ret[1])
        dataset['Unmeth.cov1'].append(ret[2])
        dataset['Meth.cov1'].append(ret[3])

        dataset['Unmeth.cov5'].append(ret[4])
        dataset['Meth.cov5'].append(ret[5])

        dataset['Unmeth.cov10'].append(ret[6])
        dataset['Meth.cov10'].append(ret[7])

    df = pd.DataFrame.from_dict(dataset)
    logger.info(df)
    outfn = os.path.join(out_dir, 'NA12878.all.chr.table.csv')
    df.to_csv(outfn)
    logger.info(f'save to {outfn}')
    pass


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Sanity check')
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by spaces', required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name1>;<file-name1>", required=True)
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--beddir', type=str, help="base dir for bed files", default=None)  # need perform performance evaluation before, then get concordant, etc. bed files, like '/projects/li-lab/yang/results/2021-04-01'
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=' ')
    parser.add_argument('--chr', type=str, help="chromosome to evaluate", default=None)
    parser.add_argument('--processors', type=int, help="running processors", default=8)
    parser.add_argument('--pred-threshold', type=float, help="prediction threshold", default=0.5)
    parser.add_argument('--min-bgtruth-cov', type=int, help="cutoff of coverage in bg-truth", default=1)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff of coverage in nanopore calls", default=1)
    parser.add_argument('--baseFormat', type=int, help="base format after imported", default=1)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--enable-cache', action='store_true')
    parser.add_argument('--using-cache', action='store_true')
    parser.add_argument('--is-report', action='store_true')
    parser.add_argument('--deepmod-paper-results', action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    set_log_debug_level()

    args = parse_arguments()

    # cache functions
    enable_cache = args.enable_cache
    using_cache = args.using_cache

    if enable_cache:
        os.makedirs(cache_dir, exist_ok=True)

    RunPrefix = args.runid.replace('SanityCheck-', '')

    ## Now, we are exporting >=1 cov results, the cov = 3,5, etc. can be later used before plotting
    # tool coverage cutoff 1, or 3, 5
    minToolCovCutt = args.toolcov_cutoff
    minToolCovCutt = 1

    # bgtruth coverage cutoff 1, or 5, 10  --min-bgtruth-cov
    bgtruthCutt = args.min_bgtruth_cov
    bgtruthCutt = 1

    # load into program format 0-base or 1-base
    baseFormat = args.baseFormat
    # Currently we only use 1-base start format, for BED of singletons, non-singletons are use 1-base format
    baseFormat = 1

    # output csv seperator: , or tab
    sep = args.sep

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir:{out_dir}')

    logger.debug(args)
    logger.info(f'\n\n####################\n\n')

    # we import multiple (1 or 2) replicates and join them
    encode, fnlist = args.bgtruth.split(':')
    fnlist = fnlist.split(';')

    if len(fnlist) > 2:
        raise Exception(f'Currently only support bgtruth with upto two, but found more: {fnlist}')

    logger.debug(f'BGTruth fnlist={fnlist}, encode={encode}')

    bgTruthList = []
    for fn in fnlist:
        if len(fn) == 0:  # incase of input like 'bismark:/a/b/c;'
            continue
        # import if cov >= 1 firstly, then after join two replicates step, remove low coverage
        bgTruth_the1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True, using_cache=using_cache, enable_cache=enable_cache)
        bgTruthList.append(bgTruth_the1)

    logger.info(f'\n\n####################\n\n')
    if args.is_report:
        report_meth_unmeth_table()
        sys.exit(0)

    if args.deepmod_paper_results:
        ## Evaluate for chr
        logger.debug(f'Start check paper results for chr={args.chr}')
        chrSet = [args.chr]
        # Combine one/two replicates using DeepMod methods
        bgTruth1, ret1 = combineBGTruthList_by_DeepModPaper(bgTruthList, covCutoff=1, filterChrs=chrSet)

        # value is (meth_indicator: int, cov)
        bgTruth5 = {key: bgTruth1[key] for key in bgTruth1 if bgTruth1[key][1] >= 5}

        for callstr in args.calls:
            callencode, callfn = callstr.split(':')
            if len(callfn) == 0:
                continue
            if callencode != 'DeepMod.Cluster':
                raise Exception("Only Deepmod paper results")
            # value is (meth_freq, cov)
            ontCall = import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache,
                                      using_cache=using_cache, include_score=False, siteLevel=True)
        callSet = set(ontCall.keys())
        bsSet = set(bgTruth5.keys())
        joinedSet = callSet.intersection(bsSet)
        logger.info(f"BS-seq(cov>=5)={len(bgTruth5):,}, DeepMod={len(ontCall):,}, intersect={len(joinedSet):,}")

        ## joined sets performance report
        y_truth = []
        y_pred = []
        for key in joinedSet:
            bs_label=bgTruth5[key][0]
            pred_label = 1 if ontCall[key][0] >= args.pred_threshold else 0
            y_truth.append(bs_label)
            y_pred.append(pred_label)
        ## Calculate precision and recall
        prec = precision_score(y_truth, y_pred)
        recall = recall_score(y_truth, y_pred)
        f1 = f1_score(y_truth, y_pred)
        accuracy = accuracy_score(y_truth, y_pred)
        logger.info(f"Joined report: chr={args.chr} (Meth={sum(y_truth):,}, Unmeth={len(y_truth)-sum(y_truth)}), precision={prec}, recall={recall}, f1={f1}, accuracy={accuracy}")
        print(f"Joined report\t{args.chr}\t{sum(y_truth)}\t{len(y_truth)-sum(y_truth)}\t{prec}\t{recall}\t{f1}\t{accuracy}\n")

        y_truth = []
        y_pred = []
        for key in bgTruth5:
            bs_label = bgTruth5[key][0]
            if key in ontCall:
                pred_label = 1 if ontCall[key][0] >= args.pred_threshold else 0
            else:
                pred_label = 0
            y_truth.append(bs_label)
            y_pred.append(pred_label)
        ## Calculate precision and recall
        prec = precision_score(y_truth, y_pred)
        recall = recall_score(y_truth, y_pred)
        f1 = f1_score(y_truth, y_pred)
        accuracy = accuracy_score(y_truth, y_pred)
        logger.info(f"BS-seq report: chr={args.chr} (Meth={sum(y_truth):,}, Unmeth={len(y_truth)-sum(y_truth):,}), precision={prec}, recall={recall}, f1={f1}, accuracy={accuracy}")
        print(
            f"Bs-seq report\t{args.chr}\t{sum(y_truth)}\t{len(y_truth) - sum(y_truth)}\t{prec}\t{recall}\t{f1}\t{accuracy}\n")
        sys.exit(0)

    ## Evalutate on chromsome
    logger.info(f'\n\n####################\n\n')
    logger.info(f'Evalutate on chr={args.chr}')
    chrSet = [args.chr]
    # Combine one/two replicates using DeepMod methods
    bgTruth1, ret1 = combineBGTruthList_by_DeepModPaper(bgTruthList, covCutoff=1, filterChrs=chrSet)
    bgTruth5 = {key: bgTruth1[key] for key in bgTruth1 if bgTruth1[key][1] >= 5}
    bgTruth10 = {key: bgTruth1[key] for key in bgTruth1 if bgTruth1[key][1] >= 10}

    ## Load and evalutate on calls
    dataset = []
    for callstr in args.calls:
        callencode, callfn = callstr.split(':')
        if len(callfn) == 0:
            continue
        if callencode in ['DeepMod.C', 'DeepMod.Cluster', 'Guppy']:
            ontCall = import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache, using_cache=using_cache, include_score=False, siteLevel=True)
        else:
            ontCall_ReadLevel = import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache, using_cache=using_cache, include_score=False, siteLevel=False)
            ontCall = readLevelToSiteLevelWithCov(ontCall_ReadLevel, minCov=1, toolname=callencode)

        ret = report_performance_deepmod(ontCall, bgTruth5, threshold=args.pred_threshold)
        logger.info(f'ret={ret}')
        ret.update({"coverage": 5, "chr": args.chr, "tool": callencode})
        dataset.append(ret)

    df = pd.DataFrame(dataset)
    logger.info(df)

    outfn = os.path.join(out_dir, f'sanity.check.using.DeepMod.paper.performance.chr_{args.chr}.threshold_{args.pred_threshold:.2f}.csv')
    df.to_csv(outfn)
    logger.info(f'save to {outfn}')
    logger.info("### Sanity check DONE")
