#!/usr/bin/env python3

"""
Generate site-level methylation correlation results in Nanocompare paper.
"""
import argparse

from nanocompare.eval_common import *
from nanocompare.global_settings import Top3ToolNameList, location_filename_to_abbvname, get_tool_name


def summary_cpgs_stats_results_table():
    """
    Study and summary each tool joined with bg-truth results, make table as dataframe
    :return:
    """
    logger.debug(f"Report number of sites by methylation calling tools in each region, take times")
    dataset = []
    bgtruthCpGs = set(list(bgTruth.keys()))
    joinedSet = None
    unionSet = set()
    for toolname in loaded_callname_list:
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
            logger.info(f'We use Concordant and Discordant BED file at basedir={args.beddir}')
            datasetBedDir = args.beddir
            concordantFileName = find_bed_filename(basedir=datasetBedDir, pattern=f'{args.dsname}*hg38_nonsingletons*.concordant.bed')
            concordantSet = filter_cpgkeys_using_bedfile(callSet, concordantFileName)

            discordantFileName = find_bed_filename(basedir=datasetBedDir, pattern=f'{args.dsname}*hg38_nonsingletons*.discordant.bed')
            discordantSet = filter_cpgkeys_using_bedfile(callSet, discordantFileName)

            ret.update({'Concordant': len(concordantSet), 'Discordant': len(discordantSet)})

        dataset.append(ret)

        # logger.info(ret)

        logger.info(f'BG-Truth join with {toolname} get {len(toolOverlapBGTruthCpGs):,} CpGs')
        # outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-bgtruth-{name1}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
        # save_keys_to_bed(overlapCpGs, outfn)

    # add additional rows for Joined count
    ret = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(joinedSet)}
    dataset.append(ret)

    # add additional row for Unioned count
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

    df = pd.DataFrame(dataset, index=loaded_callname_list + ['Joined', 'Union', 'TOP3 Joined', 'TOP3 Union'])

    logger.info(df)

    outfn = os.path.join(out_dir, f'{RunPrefix}-summary-bgtruth-tools-bsCov{bgtruthCutt}-minCov{minToolCovCutt}.csv')
    df.to_csv(outfn, sep=args.sep)
    logger.info(f'save to {outfn}\n')


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

    # Output header
    for name in loaded_callname_list:
        header_list.extend([f'{name}_freq', f'{name}_cov'])
    outfile.write(sep.join(header_list))
    outfile.write("\n")

    for cpg in reportCpGSet:
        if baseFormat == 0:
            end = cpg[1] + 1
        elif baseFormat == 1:
            end = cpg[1]

        # Ouput ground truth
        row_list = [cpg[0], str(cpg[1]), str(end), f'{bgTruth[cpg][0]:.3f}', str(bgTruth[cpg][1]), cpg[2]]

        # Output each tool results
        for name in loaded_callname_list:
            if cpg in callresult_dict[name][1]:  # if cpg is in tool results
                row_list.extend([f'{callresult_dict[name][1][cpg][0]:.3f}', f'{callresult_dict[name][1][cpg][1]}'])
            else:  # if cpg key is not exist, we use NA as ''
                row_list.extend(['', ''])
        outfile.write(sep.join(row_list))
        outfile.write("\n")
    outfile.close()
    logger.info(f"save to {outfn}\n")


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Site level correlation analysis')
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by spaces', required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name1>;<file-name1>", required=True)
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--beddir', type=str, help="base dir for bed files", default=None)  # need perform performance evaluation before, then get concordant, etc. bed files, like '/projects/li-lab/yang/results/2021-04-01'
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=' ')
    parser.add_argument('--processors', type=int, help="running processors", default=8)
    parser.add_argument('--min-bgtruth-cov', type=int, help="cutoff of coverage in bg-truth", default=1)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff of coverage in nanopore calls", default=1)
    parser.add_argument('--baseFormat', type=int, help="base format after imported", default=1)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--enable-cache', action='store_true')
    parser.add_argument('--using-cache', action='store_true')

    return parser.parse_args()


def output_dict_to_bed(dictCalls, outfn, sep='\t'):
    """
    Assume dictCalls are key->value, key=(chr, 123, +), value=[(freq, cov), ...], note is 1-based format
    Output is format: 0-based format for analysis
    chr start end . . freq  cov

    :param dictCalls:
    :param outfn:
    :return:
    """
    with open(outfn, 'w') as outf:
        for key in dictCalls:
            strlist = [key[0], str(key[1] - 1), str(key[1]), '.', '.', key[2], str(dictCalls[key][0]), str(dictCalls[key][1])]
            outf.write(sep.join(strlist) + '\n')
    logger.debug(f'Output for TSS analysis: {outfn}')


if __name__ == '__main__':
    set_log_debug_level()

    args = parse_arguments()

    # cache function same with read level
    enable_cache = args.enable_cache
    using_cache = args.using_cache

    # if enable_cache:
    #     os.makedirs(cache_dir, exist_ok=True)

    # runid is always like 'MethCorr-K562_WGBS_2Reps', remove first word as RunPrefix like K562_WGBS_2Reps
    RunPrefix = args.runid.replace('TSS-', '')

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

    # Add logging files also to result output dir
    add_logging_file(os.path.join(out_dir, 'run-results.log'))

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
        bgTruth1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True, using_cache=using_cache, enable_cache=enable_cache)
        bgTruthList.append(bgTruth1)

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

        ontCall = import_call(callfn, callencode, baseFormat=baseFormat, enable_cache=enable_cache, using_cache=using_cache, include_score=False, deepmod_cluster_freq_cov_format=True)
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

    summary_cpgs_stats_results_table()

    logger.info("### Site level correlation analysis DONE")
