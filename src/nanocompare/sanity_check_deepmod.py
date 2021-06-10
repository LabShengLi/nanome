#!/usr/bin/env python3

"""
Sanity check NA12878 for DeepMod.
"""
import argparse
from multiprocessing import Pool

from nanocompare.eval_common import *


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
    from sklearn.metrics import accuracy_score
    from sklearn.metrics import f1_score

    f1_score = f1_score(y_truth, y_pred)
    accuracy = accuracy_score(y_truth, y_pred)
    ret = {"Precision": prec, "Recall": recall,
            "F1"      : f1_score, "Accuracy": accuracy, "Threshold": threshold,
            "Unmeth"  : unmethCnt, "Meth": methCnt}

    return ret


def report_meth_unmeth_table_by_chr(chr):
    logger.debug(f'Start for chr={chr}')
    chrSet = [chr]
    # Combine one/two replicates using DeepMod methods
    bgTruth1, ret1 = combineBGTruthListUsingDeepMod(bgTruthList, covCutoff=1, filterChrs=chrSet)

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

    # for chr in humanChrSet:
    #     if chr in ['chr22', 'chrY']:
    #         continue
    #     chrSet = [chr]
    #     # Combine one/two replicates using DeepMod methods
    #     bgTruth1, ret1 = combineBGTruthListUsingDeepMod(bgTruthList, covCutoff=1, filterChrs=chrSet)
    #
    #     bgTruth5 = {bgTruth1[key] for key in bgTruth1 if bgTruth1[key][1] >= 5}
    #     methCnt5 = sum([bgTruth5[key][0] for key in bgTruth5])
    #     ret5 = (methCnt5, len(bgTruth5) - methCnt5)
    #
    #     bgTruth10 = {bgTruth1[key] for key in bgTruth1 if bgTruth1[key][1] >= 10}
    #     methCnt10 = sum([bgTruth10[key][0] for key in bgTruth10])
    #     ret10 = (methCnt10, len(bgTruth10) - methCnt10)

    # bgTruth5, ret5 = combineBGTruthListUsingDeepMod(bgTruthList, covCutoff=5, filterChrs=chrSet)
    # bgTruth10, ret10 = combineBGTruthListUsingDeepMod(bgTruthList, covCutoff=10, filterChrs=chrSet)

    # dataset['chr'].append(chr)
    # dataset['Not-used.cov1'].append(ret1[2])
    # dataset['Unmeth.cov1'].append(ret1[1])
    # dataset['Meth.cov1'].append(ret1[0])
    #
    # dataset['Unmeth.cov5'].append(ret5[1])
    # dataset['Meth.cov5'].append(ret5[0])
    #
    # dataset['Unmeth.cov10'].append(ret10[1])
    # dataset['Meth.cov10'].append(ret10[0])

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

    ## Evalutate on chromsome
    logger.info(f'\n\n####################\n\n')
    logger.info(f'Evalutate on chr={args.chr}')
    chrSet = [args.chr]
    # Combine one/two replicates using DeepMod methods
    bgTruth1, ret1 = combineBGTruthListUsingDeepMod(bgTruthList, covCutoff=1, filterChrs=chrSet)
    bgTruth5 = {key: bgTruth1[key] for key in bgTruth1 if bgTruth1[key][1] >= 5}
    bgTruth10 = {key: bgTruth1[key] for key in bgTruth1 if bgTruth1[key][1] >= 10}

    ## Load and evalutate on calls
    dataset = []
    for callstr in args.calls:
        callencode, callfn = callstr.split(':')
        if len(callfn) == 0:
            continue
        ontCall = import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache, using_cache=using_cache, include_score=False, siteLevel=True)

        ret = report_performance_deepmod(ontCall, bgTruth5, threshold=args.pred_threshold)
        logger.info(f'ret={ret}')
        ret.update({"chr": args.chr, "tool": callencode})
        dataset.append(ret)

    df = pd.DataFrame(dataset)
    logger.info(df)

    outfn = os.path.join(out_dir, f'DeepMod.performance.chr_{args.chr}.threshold_{args.pred_threshold:.2f}.csv')
    df.to_csv(outfn)
    logger.info(f'save to {outfn}')

    logger.info("DONE")
    sys.exit(0)

    # Combine one/two replicates, using cutoff=1 or 5
    bgTruth = combineBGTruthList(bgTruthList, covCutoff=bgtruthCutt)

    outfn = os.path.join(out_dir, f'{RunPrefix}.tss.bgtruth.cov{bgtruthCutt}.bed')

    logger.info(f'Combined BS-seq data (cov>={bgtruthCutt}), all methylation level sites={len(bgTruth):,}')

    output_dict_to_bed(bgTruth, outfn)

    logger.info(f'\n\n####################\n\n')

    callfn_dict = defaultdict()  # callname -> filename

    # callname -> [call0, call1], call0 is no-filter results, call1 is filter by cutoff, and convert to [meth-freq, meth-cov] results.
    callresult_dict = defaultdict()
    loaded_callname_list = []

    for callstr in args.calls:
        callencode, callfn = callstr.split(':')

        if len(callfn) == 0:
            continue

        callname = get_tool_name(callencode)
        callfn_dict[callname] = callfn

        # We do now allow import DeepMod.Cluster for read level evaluation
        if callencode == 'DeepMod.C':
            raise Exception(f'{callencode} is not allowed for site level evaluation, please use DeepMod.Cluster file here')

        loaded_callname_list.append(callname)

        # For site level evaluation, only need (freq, cov) results, no score needed. Especially for DeepMod, we must import as freq and cov format from DeepMod.Cluster encode
        # TODO: cov=1 will lead to too large size of dict objects, do we really report cov=1 results?
        # Do not filter bgtruth, because we use later for overlapping (without bg-truth)

        ontCall = import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache, using_cache=using_cache, include_score=False, siteLevel=True)
        ontCallWithCov = readLevelToSiteLevelWithCov(ontCall, minCov=minToolCovCutt, toolname=callname)
        callresult_dict[callname] = ontCallWithCov
        outfn = os.path.join(out_dir, f'{RunPrefix}.tss.{callname}.cov{minToolCovCutt}.bed')
        output_dict_to_bed(ontCallWithCov, outfn)

    logger.info(f'\n\n####################\n\n')
    logger.info("TSS DONE")
    sys.exit(0)

    logger.info('Overlapping analysis start:')
    logger.info(f"Start gen venn data for each tool (cov>={minToolCovCutt})")

    # Study five set venn data, no join with bgtruth, tool-cov > tool-cutoff=1 or 3
    if len(loaded_callname_list) >= 5:
        cpg_set_dict = defaultdict()
        for callname in ToolNameList:
            cpg_set_dict[callname] = set(callresult_dict[callname][1].keys())  # .intersection(set(bgTruth.keys()))
        gen_venn_data(cpg_set_dict, namelist=ToolNameList, outdir=out_dir, tagname=f'{RunPrefix}.{args.dsname}.five.tools.cov{minToolCovCutt}')

    logger.info(f"Start gen venn data for TOP3 tools (cov>={minToolCovCutt})")
    # Study top3 tool's venn data, no join with bgtruth, tool-cov > tool-cutoff=3
    top3_cpg_set_dict = defaultdict()
    for callname in Top3ToolNameList:
        top3_cpg_set_dict[callname] = set(callresult_dict[callname][1].keys())
    gen_venn_data(top3_cpg_set_dict, namelist=Top3ToolNameList, outdir=out_dir, tagname=f'{RunPrefix}.{args.dsname}.top3.cov{minToolCovCutt}')

    logger.info(f'\n\n####################\n\n')

    logger.info(f"Start getting intersection (all joined) sites by tools and bgtruth")
    coveredCpGs = set(list(bgTruth.keys()))
    for name in loaded_callname_list:
        coveredCpGs = coveredCpGs.intersection(set(list(callresult_dict[name][1].keys())))
        logger.info(f'Join {name} get {len(coveredCpGs)} CpGs')
    logger.info(f"Reporting {len(coveredCpGs)} CpGs are covered by all tools and bgtruth")

    logger.info('Output data of coverage and meth-freq on joined CpG sites for correlation analysis')

    outfn_joined = os.path.join(out_dir, f"Meth_corr_plot_data_joined-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.csv")
    save_meth_corr_data(callresult_dict, bgTruth, coveredCpGs, outfn_joined)

    outfn_bgtruth = os.path.join(out_dir, f"Meth_corr_plot_data_bgtruth-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.csv")
    save_meth_corr_data(callresult_dict, bgTruth, set(list(bgTruth.keys())), outfn_bgtruth)

    # Report correlation matrix here
    df = pd.read_csv(outfn_joined, sep=sep)
    df = df.filter(regex='_freq$', axis=1)
    cordf = df.corr()
    logger.info(f'Correlation matrix is:\n{cordf}')
    corr_outfn = os.path.join(out_dir, f'Meth_corr_plot_data-{RunPrefix}-correlation-matrix.xlsx')
    cordf.to_excel(corr_outfn)

    logger.info(f'\n\n####################\n\n')

    # plot fig5a of correlation plot
    command = f"set -x; PYTHONPATH=$NanoCompareDir/src python $NanoCompareDir/src/plot_figure.py fig5a -i {outfn_joined} -o {out_dir}"
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")

    logger.info(f'\n\n####################\n\n')

    logger.info("### Site level correlation analysis DONE")
