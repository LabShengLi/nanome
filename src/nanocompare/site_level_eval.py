#!/home/liuya/anaconda3/envs/nmf/bin/python

"""
Generate site-level methy correlation plotting based on input data as a tsv

All usedful functions are located in nanocompare.meth_stats.meth_stats_common
"""
import argparse
import subprocess

from nanocompare.global_settings import get_tool_name, Top3ToolNameList, ToolNameList, location_filename_to_abbvname
from nanocompare.eval_common import *


def summary_cpgs_stats_results_table():
    """
    Study and summary each tool joined with bg-truth results, make table as dataframe
    :return:
    """
    logger.debug(f"Study set intersection of each tool with bgtruth")
    dataset = []
    bgtruthCpGs = set(list(bgTruth.keys()))
    joinedSet = None
    unionSet = set()
    for toolname in callname_list:
        ## CpG sites set with cov >= cutoff(3)
        logger.info(f'Study tool={toolname}')
        callSet = set(list(callresult_dict[toolname][1].keys()))
        if not joinedSet:
            joinedSet = set(callSet)
        else:
            joinedSet = joinedSet.intersection(callSet)
        unionSet = unionSet.union(callSet)
        toolOverlapBGTruthCpGs = bgtruthCpGs.intersection(set(list(callresult_dict[toolname][1].keys())))
        ret = {f'CpG sites in BG-Truth cov>={bgtruthCutt}': len(bgtruthCpGs), 'Total CpG sites by Nanopore tool': len(callresult_dict[toolname][0]), f'Total CpG sites by tool cov>={minToolCovCutt}': len(callresult_dict[toolname][1]), 'Joined CpG sites with BG-Truth': len(toolOverlapBGTruthCpGs)}

        cnt_calls = 0
        for cpg in callresult_dict[toolname][0]:
            cnt_calls += len(callresult_dict[toolname][0][cpg])
        ret.update({'Total calls by Nanopore reads': cnt_calls})

        # Add coverage of every regions by each tool here
        for bedfn in narrowCoordFileList[1:]:  # calculate how overlap with Singletons, Non-Singletons, etc.
            basefn = os.path.basename(bedfn)
            tagname = location_filename_to_abbvname[basefn]
            subset = filter_cpgkeys_using_bedfile(callSet, bedfn)
            ret.update({tagname: len(subset)})

        if args.beddir:  # add concordant and discordant region coverage if needed
            datasetBedDir = args.beddir
            concordantFileName = find_bed_filename(basedir=datasetBedDir, pattern=f'{args.dsname}*hg38_nonsingletons.concordant.bed')
            concordantSet = filter_cpgkeys_using_bedfile(callSet, concordantFileName)

            discordantFileName = find_bed_filename(basedir=datasetBedDir, pattern=f'{args.dsname}*hg38_nonsingletons.discordant.bed')
            discordantSet = filter_cpgkeys_using_bedfile(callSet, discordantFileName)

            ret.update({'Concordant': len(concordantSet), 'Discordant': len(discordantSet)})

        dataset.append(ret)

        logger.info(ret)

        logger.info(f'BG-Truth join with {toolname} get {len(toolOverlapBGTruthCpGs):,} CpGs')
        # outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-bgtruth-{name1}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
        # save_keys_to_bed(overlapCpGs, outfn)

    ret = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(joinedSet)}
    dataset.append(ret)

    ret = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(unionSet)}
    dataset.append(ret)

    # also report top3 joined and union set
    top3JointSet = None
    top3UnionSet = set()

    for callname in Top3ToolNameList:
        toolSet = top3_cpg_set_dict[callname]
        if not top3JointSet:
            top3JointSet = toolSet
        else:
            top3JointSet = top3JointSet.intersection(toolSet)
        top3UnionSet = top3UnionSet.union(toolSet)

    ret = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top3JointSet)}
    dataset.append(ret)

    ret = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top3UnionSet)}
    dataset.append(ret)

    df = pd.DataFrame(dataset, index=callname_list + ['Joined', 'Union', 'TOP3 Joined', 'TOP3 Union'])

    logger.info(df)

    # df = df.iloc[:, [4, 0, 1, 2, 5, 6, 7, 8, 9, 3]]

    outfn = os.path.join(out_dir, f'{RunPrefix}-summary-bgtruth-tools-bsCov{bgtruthCutt}-minCov{minToolCovCutt}.csv')
    df.to_csv(outfn, sep=args.sep)
    logger.info(f'save to {outfn}\n')


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Correlation data gen task')
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by spaces', required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name1>;<file-name1>", required=True)
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--beddir', type=str, help="base dir for bed files", default=None)  # need perform performance evaluation before, then get concordant, etc. bed files, like '/projects/li-lab/yang/results/2021-04-01'
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=',')
    parser.add_argument('--processors', type=int, help="running processors", default=8)
    parser.add_argument('--bgtruthcov-cutoff', type=int, help="cutoff of coverage in bg-truth", default=5)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff of coverage in nanopore calls", default=3)
    parser.add_argument('--baseFormat', type=int, help="base format after imported", default=1)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--enable-cache', action='store_true')
    parser.add_argument('--using-cache', action='store_true')

    return parser.parse_args()


def save_meth_corr_data(callresult_dict, bgTruth, reportCpGSet, outfn):
    """
    Save meth freq and cov results into csv file
    :param callresult_dict:
    :param bgTruth:
    :param reportCpGSet:
    :param outfn:
    :return:
    """
    outfile = open(outfn, 'w')

    header_list = ['chr', 'start', 'end', 'BGTruth_freq', 'BGTruth_cov', 'strand']

    for name in callname_list:
        header_list.extend([f'{name}_freq', f'{name}_cov'])
    outfile.write(sep.join(header_list))
    outfile.write("\n")

    for cpg in reportCpGSet:
        if baseFormat == 0:
            end = cpg[1] + 1
        elif baseFormat == 1:
            end = cpg[1]

        row_list = [cpg[0], str(cpg[1]), str(end), f'{bgTruth[cpg][0]:.3f}', str(bgTruth[cpg][1]), cpg[2]]

        for name in callname_list:
            if cpg in callresult_dict[name][1]:  # if cpg is in tool results
                row_list.extend([f'{callresult_dict[name][1][cpg][0]:.3f}', f'{callresult_dict[name][1][cpg][1]}'])
            else:  # if cpg key is not exist, we use NA as ''
                row_list.extend(['', ''])
        outfile.write(sep.join(row_list))
        outfile.write("\n")
    outfile.close()
    logger.info(f"save to {outfn}\n")


if __name__ == '__main__':
    set_log_debug_level()

    args = parse_arguments()

    enable_cache = args.enable_cache
    using_cache = args.using_cache

    if enable_cache:
        os.makedirs(cache_dir, exist_ok=True)

    RunPrefix = args.runid.replace('MethCorr-', '')

    # tool coverage cutoff 3
    minToolCovCutt = args.toolcov_cutoff

    # bgtruth coverage cutoff 5
    bgtruthCutt = args.bgtruthcov_cutoff

    # load into program format 0-base or 1-base
    baseFormat = args.baseFormat

    # output csv seperator: , or tab
    sep = args.sep

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir:{out_dir}')

    # Add logging files also to result output dir
    add_logging_file(os.path.join(out_dir, 'run-results.log'))

    logger.debug(args)

    callfn_dict = defaultdict()  # callname -> filename

    # callname -> [call0, call1], call0 is no-filter results, call1 is filter by cutoff, and convert to [meth-freq, meth-cov] results.
    callresult_dict = defaultdict()
    callname_list = []

    for callstr in args.calls:
        callencode, callfn = callstr.split(':')
        callname = get_tool_name(callencode)
        callfn_dict[callname] = callfn
        callname_list.append(callname)
        callresult_dict[callname] = [import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache, using_cache=using_cache)]
    logger.debug(callfn_dict)

    # we can import multiple replicates and join them
    encode, fnlist = args.bgtruth.split(':')
    fnlist = fnlist.split(';')
    logger.debug(f'BGTruth fnlist={fnlist}, encode={encode}')

    bgTruthList = []
    for fn in fnlist:
        # import if cov >= 1 firstly, then after join two replicates step, remove low coverage
        bgTruth1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True, using_cache=using_cache, enable_cache=enable_cache)
        bgTruthList.append(bgTruth1)

    # bgTruth = import_bgtruth(fn, encode, covCutoff=bgtruthCutt, baseFormat=baseFormat, enable_cache=enable_cache, using_cache=using_cache)

    bgTruth = combineBGTruthList(bgTruthList, covCutoff=bgtruthCutt)

    # Filter cov of nanopore tools
    for callname in callname_list:
        callresult_dict[callname].append(coverageFiltering(callresult_dict[callname][0], minCov=minToolCovCutt, toolname=callname))

    logger.debug(callfn_dict)
    logger.debug(callname_list)

    logger.info(f"Start gen venn data for each tool (cov>={minToolCovCutt}) with or without bgtruth (cov>={bgtruthCutt})")
    # Study five set venn data, no join with bgtruth, tool-cov > tool-cutoff=3
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

    logger.info(f"Start set intersection with all tools joined together with bgtruth")

    coveredCpGs = set(list(bgTruth.keys()))
    for name in callname_list:
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

    # plot fig5a of correlation plot
    command = f"set -x; PYTHONPATH=$NanoCompareDir/src python $NanoCompareDir/src/plot_figure.py fig5a -i {outfn_joined} -o {out_dir}"
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")

    summary_cpgs_stats_results_table()

    logger.info("### Site level correlation data for plotting generation DONE")
