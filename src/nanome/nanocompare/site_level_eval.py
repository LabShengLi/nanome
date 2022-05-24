#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : site_level_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

"""
Generate site-level methylation correlation results in nanome paper.
"""

import argparse

import pybedtools
from scipy import stats

from nanome.common.eval_common import *
from nanome.common.global_settings import ToolNameList, NANOME_VERSION, load_genome_annotation_config, sing_tagname, \
    nonsing_tagname, save_done_file


def get_nsites_in_regions(callSet, bedfn, tagname):
    subset = filter_cpgkeys_using_bedfile(callSet, bedfn)
    ret = {tagname: len(subset)}
    return ret


def summary_cpgs_stats_results_table():
    """
    Study and summary each tool joined with bg-truth results, make table as dataframe
    :return:
    """

    logger.debug(f"Report number of sites by methylation calling tools in each region, take times...")

    eval_cov_summary_region = region_bed_list[1:]
    logger.debug(f"Evaluated on regions: {eval_cov_summary_region}")

    dataset = []
    bgtruthCpGs = set(list(bgTruth.keys())) if bgTruth else set()
    joinedSet = None
    unionSet = set()

    retList = []
    logger.debug("Start to study coverage now:")
    for toolname in loaded_callname_list:
        ## CpG sites set with cov >= cutoff(3)
        logger.info(f'Study CPG coverage for tool={toolname}')
        logger.info(f"Memory report: {get_current_memory_usage()}")
        callSet = list(set(callresult_dict_cov3[toolname].keys()))

        if not joinedSet:
            joinedSet = set(callSet)
        else:
            joinedSet = joinedSet.intersection(callSet)
        unionSet = unionSet.union(callSet)
        toolOverlapBGTruthCpGs = bgtruthCpGs.intersection(set(list(callresult_dict_cov3[toolname].keys())))
        row_dict = {f'CpG sites in BG-Truth cov>={bgtruthCutt}': len(bgtruthCpGs),
                    'Total CpG sites by Nanopore tool': call_cov1_cpg_sites[toolname],
                    f'Total CpG sites by tool cov>={minToolCovCutt}': len(callresult_dict_cov3[toolname]),
                    'Joined CpG sites with BG-Truth': len(toolOverlapBGTruthCpGs)}
        row_dict.update({'Total calls by Nanopore reads': call_cov1_calls[toolname]})

        callBed = calldict2bed(callSet)

        if not args.mpi:
            # Add coverage of every regions by each tool here
            bar = tqdm(eval_cov_summary_region)
            for (bedfn, tagname, region_bed) in bar:  # calculate how overlap with Singletons, Non-Singletons, etc.
                bar.set_description(f"CPG_cov-{args.dsname}-{toolname}-region-{tagname}")
                if not args.large_mem and region_bed is None:  # load in demand
                    region_bed = get_region_bed_tuple(bedfn,
                                                      enable_base_detection_bedfile=not args.disable_bed_check,
                                                      enable_cache=args.enable_cache, using_cache=args.using_cache,
                                                      cache_dir=ds_cache_dir
                                                      )[2]

                if region_bed is None:
                    logger.debug(f"region name={tagname} is not found")
                    continue
                intersect_bed = intersect_bed_regions(callBed, region_bed, bedfn)
                ret = {tagname: len(intersect_bed)}
                retList.append(ret)
        else:
            # multi-threading way
            # if some entry is strange 0, means memory is not enough, such as 50G for HL60
            global progress_bar_global_site
            progress_bar_global_site = tqdm(total=len(eval_cov_summary_region))
            progress_bar_global_site.set_description(f"MT-CPG_cov-{args.dsname}-{toolname}-all-regions")
            executor = ThreadPoolExecutor(max_workers=args.processors)
            all_tasks = []
            tag_list = []
            for (bedfn, tagname, region_bed) in eval_cov_summary_region:
                if not args.large_mem and region_bed is None:  # load in demand
                    region_bed = get_region_bed_tuple(bedfn,
                                                      enable_base_detection_bedfile=not args.disable_bed_check,
                                                      enable_cache=args.enable_cache, using_cache=args.using_cache,
                                                      cache_dir=ds_cache_dir
                                                      )[2]

                if region_bed is None:
                    logger.debug(f"region name={tagname} is not found")
                    continue
                future = executor.submit(get_num_intersect, callBed, region_bed, bedfn=bedfn, tagname=tagname)
                future.add_done_callback(update_progress_bar_site_level)
                all_tasks.append(future)
                tag_list.append(tagname)
            executor.shutdown()
            progress_bar_global_site.close()

            for future, tagname in zip(all_tasks, tag_list):
                ret = future.result()
                ## del bed_ret
                retList.append(ret)

        if concordant_bed is not None:
            intersect_bed = intersect_bed_regions(callBed, concordant_bed)
            ret = {'Concordant': len(intersect_bed)}
            retList.append(ret)

        if discordant_bed is not None:
            intersect_bed = intersect_bed_regions(callBed, discordant_bed)
            ret = {'Discordant': len(intersect_bed)}
            retList.append(ret)

        retDict = {}
        for e in retList:
            retDict.update(e)
        logger.debug(f"retDict={retDict}")

        ## Sanity check
        if sing_tagname in retDict and nonsing_tagname in retDict:
            sum_sing_nonsingle = retDict[sing_tagname] + retDict[nonsing_tagname]
        else:
            sum_sing_nonsingle = None
        if 'CG_20' in retDict and 'CG_40' in retDict and 'CG_60' in retDict and 'CG_80' in retDict and 'CG_100' in retDict:
            sum_cg = retDict['CG_20'] + retDict['CG_40'] + retDict['CG_60'] + retDict['CG_80'] + retDict['CG_100']
        else:
            sum_cg = None
        total_sites = row_dict[f'Total CpG sites by tool cov>={minToolCovCutt}']

        if sum_sing_nonsingle and sum_cg:
            logger.debug(
                f"\n\nSanity check: sum_sing_nonsingle={sum_sing_nonsingle:,}; sum_cg={sum_cg:,}; total={total_sites:,}")

        if sum_sing_nonsingle is not None and sum_sing_nonsingle != total_sites:
            logger.debug(
                f"Sanity check for {toolname}, total_sites={total_sites:,}, sum_sing_nonsingle={sum_sing_nonsingle:,}, some non-singletons are not captered by bed file")
            retDict[nonsing_tagname] = total_sites - retDict[sing_tagname]
            logger.debug(f"Updated, retDict={retDict}")

        ## add retDict into row_dict
        row_dict.update(retDict)
        dataset.append(row_dict)
        logger.debug(f'BG-Truth join with {toolname} get {len(toolOverlapBGTruthCpGs):,} CpGs')
        logger.debug(f"Memory report: {get_current_memory_usage()}")

    # add additional rows for Joined count
    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(joinedSet)}
    dataset.append(new_row_dict)

    # add additional row for Unioned count
    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(unionSet)}
    dataset.append(new_row_dict)

    # also report top3 joined and union set
    top3JointSet = None
    top3UnionSet = set()

    for callname in ToolNameList[:3]:
        if callname not in top3_cpg_set_dict:
            continue
        toolSet = top3_cpg_set_dict[callname]
        if not top3JointSet:
            top3JointSet = toolSet
        else:
            top3JointSet = top3JointSet.intersection(toolSet)
        top3UnionSet = top3UnionSet.union(toolSet)

    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top3JointSet) if top3JointSet else None}
    dataset.append(new_row_dict)

    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top3UnionSet) if top3UnionSet else None}
    dataset.append(new_row_dict)

    # also report top4 joined and union set
    top4JointSet = None
    top4UnionSet = set()

    for callname in ToolNameList[:4]:
        if callname not in callresult_dict_cov3:
            continue
        toolSet = set(callresult_dict_cov3[callname].keys())
        if not top4JointSet:
            top4JointSet = toolSet
        else:
            top4JointSet = top4JointSet.intersection(toolSet)
        top4UnionSet = top4UnionSet.union(toolSet)

    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top4JointSet) if top4JointSet else None}
    dataset.append(new_row_dict)

    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top4UnionSet) if top4UnionSet else None}
    dataset.append(new_row_dict)

    df = pd.DataFrame(dataset, index=loaded_callname_list +
                                     ['Joined', 'Union', 'TOP3 Joined', 'TOP3 Union', 'TOP4 Joined', 'TOP4 Union'])

    logger.debug(df)

    outfn = os.path.join(out_dir,
                         f'{RunPrefix}-summary-bgtruth-tools-bsCov{bgtruthCutt}-minCov{minToolCovCutt}.table.s10.xlsx')
    df.to_excel(outfn)
    logger.debug(f'save to {outfn}\n')
    logger.debug(f"Memory report: {get_current_memory_usage()}")


def correlation_report_on_regions(corr_infn, bed_tuple_list, dsname=None, runid=None, outdir=None, join_tag="join",
                                  mpi=True):
    """
    Calculate Pearson's correlation coefficient at different regions.
    :param corr_infn:
    :param beddir:
    :param dsname:
    :param outdir:
    :return:
    """
    global progress_bar_global_site
    progress_bar_global_site = tqdm(total=len(bed_tuple_list))
    progress_bar_global_site.set_description(f"MT-PCC-{join_tag}-{dsname}-regions")

    if mpi:
        executor = ThreadPoolExecutor(max_workers=args.processors)
    else:
        executor = ThreadPoolExecutor(max_workers=1)
    all_task = []
    ## input: df, bed_tuple
    ## return: list of dict [{tool1's pcc}, {toolk's pcc}]
    df = pd.read_csv(corr_infn)
    if df.isnull().values.any():
        df.fillna('.', inplace=True)

    for bed_tuple in bed_tuple_list:
        future = executor.submit(compute_pcc_at_region, df, bed_tuple)
        future.add_done_callback(update_progress_bar_site_level)
        all_task.append(future)
    executor.shutdown()
    progress_bar_global_site.close()

    ret_list = []  # list of dict for dataframe
    for future in all_task:
        ret_l1 = future.result()
        if ret_l1 is None:
            continue
        ret_list.extend(ret_l1)

    # logger.info(dataset)
    outdf = pd.DataFrame(ret_list)
    logger.debug(outdf)

    outfn = os.path.join(outdir, f'{runid}_{dsname}.corrdata.coe.pvalue.each.regions_{join_tag}.xlsx')
    outdf.to_excel(outfn)
    logger.debug(f'save to {outfn}')
    return outdf


def save_meth_corr_data(callresult_dict, bgTruth, reportCpGSet, outfn):
    """
    Save meth freq and cov results into csv file
    :param callresult_dict:
    :param bgTruth:
    :param reportCpGSet:
    :param outfn:
    :return:
    """
    # outfile = open(outfn, 'w')
    outfile = gzip.open(outfn, 'wt')

    header_list = ['Chr', 'Start', 'End', 'BGTruth_freq', 'BGTruth_cov', 'Strand']

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
            if cpg in callresult_dict[name]:  # if cpg is in tool results
                row_list.extend([f'{callresult_dict[name][cpg][0]:.3f}', f'{callresult_dict[name][cpg][1]}'])
            else:  # if cpg key is not exist, we use NA as ''
                row_list.extend(['', ''])
        outfile.write(sep.join(row_list))
        outfile.write("\n")
    outfile.close()
    logger.debug(f"save to {outfn}\n")


def update_progress_bar_site_level(*a):
    """
    Update progress for multiprocessing
    :param a:
    :return:
    """
    global progress_bar_global_site
    progress_bar_global_site.update()


def get_num_intersect(callBed, region_bed, bedfn="", tagname=None):
    ret_bed = intersect_bed_regions(callBed, region_bed, bedfn=bedfn)
    if len(ret_bed) < 1:
        logger.error(
            f"Found 0 intersect CPG for tool with {tagname} region, callBed={callBed}, region_bed={region_bed}")
    ret = {tagname: len(ret_bed)}
    return ret


def compute_pcc_at_region(df, bed_tuple):
    """
    Compute PCC of input DF with respect to bed region tuple
    Args:
        df:
        bed_tuple:

    Returns:

    """
    infn, tagname, coord_bed = bed_tuple
    logger.debug(f'tagname={tagname}, coord_fn={infn}')
    if not args.large_mem and tagname != genome_wide_tagname and coord_bed is None:  # load on demand
        eval_coord_bed = get_region_bed_tuple(
            infn, enable_base_detection_bedfile=enable_base_detection_bedfile,
            enable_cache=args.enable_cache, using_cache=args.using_cache,
            cache_dir=ds_cache_dir)[2]
    else:  # large memory, or genome wide - None
        eval_coord_bed = coord_bed

    if tagname != genome_wide_tagname and eval_coord_bed is None:
        logger.debug(f"Region name={tagname} is not found, not compute PCC")
        return None

    newdf = filter_corrdata_df_by_bedfile(df, eval_coord_bed, infn)
    if newdf is None:
        logger.debug(f"Found intersection 0 CPGs for tagname={tagname}, no report for PCC")
        return None

    # Computer COE and pvalue
    newdf = newdf.filter(regex='_freq$', axis=1)
    ret_list = []
    for i in range(1, len(newdf.columns)):
        toolname = str(newdf.columns[i]).replace('_freq', '')
        try:  # too few samples will fail
            # warnings.filterwarnings('ignore', category=PearsonRConstantInputWarning)
            mask_notna = newdf.iloc[:, i].notna().values
            coe, pval = stats.pearsonr(newdf.iloc[mask_notna, 0], newdf.iloc[mask_notna, i].astype(np.float64))
        except:  ## May be pearsonr function failed
            coe, pval = None, None

        # report to dataset
        ret = {
            'dsname': dsname,
            'Tool': toolname,
            'Location': tagname,
            '#Bases': newdf.iloc[:, i].notna().sum(),
            'COE': coe,
            'p-value': pval
        }
        ret_list.append(ret)
    logger.debug(f"tagname={tagname}, pcc_return={ret_list}")
    return ret_list


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='site_level_eval (NANOME)',
                                     description='Site-level correlation analysis in nanome paper')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{NANOME_VERSION}')
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str,
                        help="running prefix/output folder name, such as MethCorr-Dataset_WGBS_2Reps",
                        required=True)
    parser.add_argument('--calls', nargs='+',
                        help='all ONT call results <tool-name>:<file-encode>:<file-name> seperated by spaces, tool-name/file-encode can be Nanopolish, Megalodon, DeepSignal, Guppy, Tombo, METEORE, DeepMod, NANOME',
                        required=True)
    parser.add_argument('--bgtruth', type=str,
                        help="background truth file <encode-type>:<file-name1>;<file-name2>, encode-type can be 'encode' or 'bismark'",
                        default=None)
    parser.add_argument('--genome-annotation', type=str,
                        help='genome annotation dir, contain BED files',
                        default=None)
    parser.add_argument('--beddir', type=str,
                        help="base dir for concordant/discordant BED files generated by read-level analysis, make sure the dsname is same",
                        default=None)
    parser.add_argument('--min-bgtruth-cov', type=int, help="cutoff for coverage in bg-truth, default is >=5",
                        default=5)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff for coverage in nanopore tools, default is >=3",
                        default=3)
    parser.add_argument('--chrSet', nargs='+', help='chromosome list, default is human chr1-22, X and Y',
                        default=HUMAN_CHR_SET)
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=',')
    parser.add_argument('--processors', type=int, help="number of processors used, default is 1", default=1)

    parser.add_argument('-o', type=str, help=f"output base dir, default is {pic_base_dir}", default=pic_base_dir)
    parser.add_argument('--gen-venn', help="if generate CpGs files for venn data analysis", action='store_true')
    parser.add_argument('--sort', help="if sort outputs", action='store_true')
    parser.add_argument('--deduplicate', help="if deduplicate not unique records needed", action='store_true')
    parser.add_argument('--summary-coverage', help="if summarize coverage at each region",
                        action='store_true')
    parser.add_argument('--region-coe-report', help="if report PCC value at each region",
                        action='store_true')
    parser.add_argument('--report-no-join', help="if output no-join report also", action='store_true')
    parser.add_argument('--enable-cache', help="if enable cache functions", action='store_true')
    parser.add_argument('--using-cache', help="if use cache files", action='store_true')
    parser.add_argument('--plot', help="if plot the correlation matrix figure", action='store_true')
    parser.add_argument('--bedtools-tmp', type=str, help=f'bedtools temp dir, default is {global_temp_dir}',
                        default=global_temp_dir)
    parser.add_argument('--cache-dir', type=str,
                        help=f'cache dir used for loading calls/bs-seq(speed up running), default is {global_cache_dir}',
                        default=global_cache_dir)
    parser.add_argument('--large-mem', help="if using large memory (>100GB) for speed up", action='store_true')
    parser.add_argument('--disable-bed-check',
                        help="if disable auto-checking the 0/1 base format for genome annotations",
                        action='store_true')
    parser.add_argument('--mpi',
                        help="if using multi-processing/threading for evaluation, it can speed-up but need more memory",
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
    bed_temp_dir = os.path.join(args.bedtools_tmp, f"{dsname}_corr")
    os.makedirs(bed_temp_dir, exist_ok=True)
    pybedtools.helpers.set_tempdir(bed_temp_dir)

    ## Set cache dir for each dataset
    if args.enable_cache or args.using_cache:
        ds_cache_dir = os.path.join(args.cache_dir, dsname)
        # os.makedirs(ds_cache_dir, exist_ok=True)
    else:
        ds_cache_dir = None

    # cache function same with read level
    enable_cache = args.enable_cache
    using_cache = args.using_cache

    # runid is always like 'MethCorr-K562_WGBS_2Reps', remove first word as RunPrefix like K562_WGBS_2Reps
    RunPrefix = args.runid.replace('MethCorr-', '')

    # tool coverage cutoff 1, or 3, 5
    minToolCovCutt = args.toolcov_cutoff

    # bgtruth coverage cutoff 1, or 5, 10  --min-bgtruth-cov
    bgtruthCutt = args.min_bgtruth_cov

    # We import and report use 1-base start format
    baseFormat = 1

    # output csv seperator: , or tab
    sep = args.sep

    out_dir = os.path.join(args.o, args.runid)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir:{out_dir}')

    # Add logging files also to result output dir
    add_logging_file(os.path.join(out_dir, 'run-results.log'))
    logger.debug(args)

    if args.config:
        load_genome_annotation_config(verbose=True)
    logger.debug(f'\n\n####################\n\n')

    if args.bgtruth:
        # we import multiple (1 or 2) replicates and join them
        encode, fnlist = args.bgtruth.split(':')
        fnlist = fnlist.split(';')

        if len(fnlist) > 2:
            raise Exception(f'Currently only support the number of bgtruth upto two, but found more: {fnlist}')

        logger.debug(f'BGTruth fnlist={fnlist}, encode={encode}')

        bgTruthList = []
        for fn in fnlist:
            if len(fn) == 0:  # incase of input like 'bismark:/a/b/c;'
                continue
            # import if cov >= 1 firstly, then after join two replicates step, remove low coverage
            bgTruth1 = import_bgtruth(fn, encode, covCutoff=1, baseFormat=baseFormat, includeCov=True,
                                      filterChr=args.chrSet,
                                      using_cache=using_cache, enable_cache=enable_cache, cache_dir=ds_cache_dir)
            bgTruthList.append(bgTruth1)

        # Combine one/two replicates, using cutoff=1 or 5
        bgTruth = combineBGTruthList(bgTruthList, covCutoff=bgtruthCutt)

        logger.info(f'Combined BS-seq data (cov>={bgtruthCutt}), all methylation level sites={len(bgTruth):,}')
        logger.debug(f"Memory report: {get_current_memory_usage()}")
        logger.debug(f'\n\n####################\n\n')
    else:
        bgTruth = None

    callfn_dict = defaultdict()  # callname -> filename

    # callname -> [call0, call1], call0 is no-filter results, call1 is filter by cutoff, and convert to [meth-freq, meth-cov] results.
    callresult_dict_cov1 = defaultdict()
    call_cov1_cpg_sites = defaultdict(int)
    call_cov1_calls = defaultdict(int)

    callresult_dict_cov3 = defaultdict()
    loaded_callname_list = []

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

        if len(callfn) == 0:
            continue

        callfn_dict[toolname] = callfn

        loaded_callname_list.append(toolname)

        # For site level evaluation, only need (freq, cov) results, no score needed. Especially for DeepMod, we must import as freq and cov format from DeepMod.Cluster encode
        # Do not filter bgtruth, because we use later for overlapping (without bg-truth)
        callresult_dict_cov1[toolname] = import_call(
            callfn, callencode, baseFormat=baseFormat, filterChr=args.chrSet,
            enable_cache=enable_cache, using_cache=using_cache,
            include_score=False, siteLevel=True, cache_dir=ds_cache_dir,
            toolname=toolname, score_cutoff=score_cutoff)

        # Stats the total cpgs and calls for each calls
        cnt_calls = 0
        for cpg in callresult_dict_cov1[toolname]:
            cnt_calls += len(callresult_dict_cov1[toolname][cpg])
        call_cov1_calls[toolname] = cnt_calls
        call_cov1_cpg_sites[toolname] = len(callresult_dict_cov1[toolname])

    logger.info(f"Import calls from tools done for toollist={list(call_cov1_cpg_sites.keys())}")
    logger.info(f"Memory report: {get_current_memory_usage()}")

    logger.debug(f'Start apply cutoff={minToolCovCutt} to methylation calls, take time')
    # Cutoff of read cov >= 1 or 3, 5 for nanopore tools
    for callname in loaded_callname_list:
        callresult_dict_cov3[callname] = readLevelToSiteLevelWithCov(callresult_dict_cov1[callname],
                                                                     minCov=minToolCovCutt, toolname=callname)
    ## Destroy cov1 for memory saving
    del callresult_dict_cov1
    logger.debug(f"Memory report: {get_current_memory_usage()}")

    logger.debug(f'\n\n####################\n\n')

    top3_cpg_set_dict = defaultdict()
    for callname in ToolNameList[:3]:
        if callname in callresult_dict_cov3:
            top3_cpg_set_dict[callname] = set(callresult_dict_cov3[callname].keys())

    if args.gen_venn:
        logger.info('CPG overlapping analysis')
        logger.debug(
            f"Start gen venn data for each tool (cov>={minToolCovCutt}) and BS-seq (cov>={args.min_bgtruth_cov})")

        # Generate all tools and bsseq covered cpgs files for set evaluation
        logger.debug("We generate sets file for each tool and bg-truth")
        venn_outdir = os.path.join(out_dir, 'venn_data')
        os.makedirs(venn_outdir, exist_ok=True)

        if bgTruth:
            bg_cpgs = bgTruth.keys()
            outfn = os.path.join(venn_outdir,
                                 f'{args.dsname}.bgtruth.cpg.sites.cov{args.min_bgtruth_cov}.setsfile.txt.gz')
            ontcalls_to_setsfile_for_venn_analysis(bg_cpgs, outfn)
            outfn_sort = outfn.replace('.setsfile.txt.gz', '.setsfile.sort.txt.gz')
            sort_set_txt_file(outfn, outfn_sort)
            os.remove(outfn)

        for callname in callresult_dict_cov3.keys():
            call_keys = callresult_dict_cov3[callname].keys()
            outfn = os.path.join(venn_outdir,
                                 f'{args.dsname}.{callname}.cpg.sites.cov{args.toolcov_cutoff}.setsfile.txt.gz')
            ontcalls_to_setsfile_for_venn_analysis(call_keys, outfn)
        logger.debug(f"Memory report: {get_current_memory_usage()}")
        logger.debug(f'\n\n####################\n\n')

        if args.sort:
            ## sort all setsfile.txt.gz
            logger.debug(f"Start to sort setsfile for venn data")
            flist = glob.glob(os.path.join(venn_outdir, '*.setsfile.txt.gz'))
            with Pool(args.processors) as p:
                input_params_list = []
                for infn in flist:
                    input_params_list.append((infn, infn.replace('setsfile.txt.gz', 'setsfile.sort.txt.gz'), args.deduplicate))
                logger.debug(f"input_params_list={input_params_list}")
                p.starmap(sort_set_txt_file, input_params_list)
            for infn in flist:
                os.remove(infn)
            logger.debug(f"Sort setsfile done")

    if bgTruth:  # Having bgtruth params, then report PCC performance
        logger.debug(f"Start getting intersection (all joined) sites by tools and bgtruth")
        coveredCpGs = set(list(bgTruth.keys()))  # joined sets, start with bs-seq
        coveredCpGs001 = set(list(bgTruth.keys()))  # no change later

        sitesDataset = defaultdict(list)

        for name in loaded_callname_list:
            coveredCpGs = coveredCpGs.intersection(set(list(callresult_dict_cov3[name].keys())))
            logger.debug(f'Join {name} get {len(coveredCpGs):,} CpGs')

            joinBSWithEachToolSet = coveredCpGs001.intersection(set(list(callresult_dict_cov3[name].keys())))
            sitesDataset['Dataset'].append(args.dsname)
            sitesDataset['Method'].append(name)
            sitesDataset[f'Sites-cov{args.toolcov_cutoff}'].append(len(callresult_dict_cov3[name]))
            sitesDataset[f'BS-seq-cov{args.min_bgtruth_cov}-all'].append(len(coveredCpGs001))
            sitesDataset[f'Join-with-BSseq-cov{args.min_bgtruth_cov}-all'].append(len(joinBSWithEachToolSet))
        df = pd.DataFrame.from_dict(sitesDataset)
        outfn = os.path.join(out_dir,
                             f'{RunPrefix}_{args.dsname}.tools.cov{args.toolcov_cutoff}.join.with.bsseq.cov{args.min_bgtruth_cov}.site.level.report.csv')
        df.to_csv(outfn)

        # Output sites report for all tools and BS-seq
        logger.info(
            f"Joined {len(coveredCpGs):,} CpGs are covered by all tools (cov >= {args.toolcov_cutoff}) and BS-seq (cov >= {args.min_bgtruth_cov})")
        logger.debug('Output data of meth-freq and coverage on joined CpG sites as datasets for correlation analysis')
        outfn_joined = os.path.join(out_dir,
                                    f"Meth_corr_plot_data_joined-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.csv.gz")
        save_meth_corr_data(callresult_dict_cov3, bgTruth, coveredCpGs, outfn_joined)
        if args.sort:
            outfn_joined_sorted = os.path.join(out_dir,
                                               f"Meth_corr_plot_data_joined-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.sorted.csv.gz")
            sort_bed_file(infn=outfn_joined, outfn=outfn_joined_sorted, has_header=True)
            os.remove(outfn_joined)
        else:
            outfn_joined_sorted = outfn_joined

        logger.debug(
            'Output data of meth-freq and coverage on bgTruth related CpG sites as datasets for correlation analysis')
        outfn_bgtruth = os.path.join(out_dir,
                                     f"Meth_corr_plot_data_bgtruth-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.csv.gz")
        save_meth_corr_data(callresult_dict_cov3, bgTruth, set(list(bgTruth.keys())), outfn_bgtruth)
        if args.sort:
            outfn_bgtruth_sorted = os.path.join(out_dir,
                                                f"Meth_corr_plot_data_bgtruth-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.sorted.csv.gz")
            sort_bed_file(infn=outfn_bgtruth, outfn=outfn_bgtruth_sorted, has_header=True)
            os.remove(outfn_bgtruth)
        else:
            outfn_bgtruth_sorted = outfn_bgtruth

        # Report correlation matrix for joined results
        df = pd.read_csv(outfn_joined_sorted, sep=sep)
        df = df.filter(regex='_freq$', axis=1)
        cordf = df.corr()
        # Count CpGs
        num_join_with_bsseq = [len(df.iloc[:, 0])] * len(df.columns)
        cordf = pd.concat(
            [cordf, pd.Series(num_join_with_bsseq, index=cordf.index).rename('CpGs_all_tools_with_BSseq')],
            axis=1)

        logger.debug(f'Correlation matrix (joined CpGs) is:\n{cordf}')
        corr_outfn = os.path.join(out_dir,
                                  f'Meth_corr_plot_data_joined-{RunPrefix}-correlation-matrix-toolcov{minToolCovCutt}-bsseqcov{bgtruthCutt}.xlsx')
        cordf.to_excel(corr_outfn)

        if args.report_no_join:
            # Report correlation matrix for no-join results
            df = pd.read_csv(outfn_bgtruth_sorted, sep=sep)
            df = df.filter(regex='_freq$', axis=1)
            cordf = df.corr()
            # Count CpGs
            num_join_with_bsseq = [len(df.iloc[:, 0])]
            for k in range(1, len(df.columns)):
                num_join_with_bsseq.append(df.iloc[:, k].notna().sum())
            cordf = pd.concat([cordf, pd.Series(num_join_with_bsseq, index=cordf.index).rename('CpGs_with_BSseq')],
                              axis=1)

            logger.debug(f'Correlation matrix (No-joined CpGs) is:\n{cordf}')
            corr_outfn = os.path.join(out_dir,
                                      f'Meth_corr_plot_data_no_join-{RunPrefix}-correlation-matrix-toolcov{minToolCovCutt}-bsseqcov{bgtruthCutt}.xlsx')
            cordf.to_excel(corr_outfn)

        logger.debug(f'\n\n####################\n\n')

    if args.region_coe_report or args.summary_coverage:
        ## load region bed list
        logger.debug("Create region bed list firstly, take times......")

        # Evaluated all region filename lists, bed objects
        # assume all files are located in args.genome_annotation dir
        annot_dir = args.genome_annotation if args.genome_annotation is not None else '.'

        # regions_full_filepath = [os.path.join(annot_dir, cofn) for cofn in narrowCoordNameList[1:]] + \
        #                         [os.path.join(annot_dir, cofn) for cofn in cg_density_coord_name_list] + \
        #                         [os.path.join(annot_dir, cofn) for cofn in rep_coord_name_list]
        # region file path from genome-wide, singletons, to genic/intergenic, cg-density, and repetitive, the concordant and discoradnt wil be discovered later
        regions_full_filepath = [None] + [os.path.join(annot_dir, cofn) for cofn in region_filename_dict.keys()]

        if args.large_mem:  # load all in memory
            region_bed_list = get_region_bed_pairs_list_mp(
                regions_full_filepath,
                processors=args.processors,
                enable_base_detection_bedfile=not args.disable_bed_check,
                enable_cache=args.enable_cache,
                using_cache=args.using_cache,
                cache_dir=ds_cache_dir)
            logger.info(f"Memory report: {get_current_memory_usage()}")
        else:  # load bed coord later
            region_bed_list = [(infn, get_region_tagname(infn), None,)
                               for infn in regions_full_filepath]

        if args.beddir:  # add concordant and discordant region coverage if needed
            logger.debug(f'We are finding Concordant and Discordant BED file at basedir={args.beddir}')
            concordantFileName = find_bed_filename(basedir=args.beddir,
                                                   pattern=f'*{args.dsname}*.concordant.bed.gz')
            concordant_bed = get_region_bed(concordantFileName) if concordantFileName is not None else None

            discordantFileName = find_bed_filename(basedir=args.beddir,
                                                   pattern=f'*{args.dsname}*.discordant.bed.gz')
            discordant_bed = get_region_bed(discordantFileName) if discordantFileName is not None else None
        else:
            concordant_bed = None
            discordant_bed = None
        ## Add concordant/discordant if possible
        if concordant_bed is not None:
            region_bed_list += [(concordantFileName, 'Concordant', concordant_bed,)]
        if discordant_bed is not None:
            region_bed_list += [(discordantFileName, 'Discordant', discordant_bed,)]

        logger.debug(f"Evaluated on regions: {region_bed_list}")

    if args.region_coe_report and bgTruth:
        eval_genomic_context_tuple = region_bed_list

        # file like: Meth_corr_plot_data_joined-TestData_RRBS_2Reps-bsCov1-minToolCov1-baseFormat1.sorted.csv.gz
        fnlist = glob.glob(os.path.join(out_dir,
                                        f'Meth_corr_plot_data_joined-*bsCov{args.min_bgtruth_cov}-minToolCov{args.toolcov_cutoff}*.csv.gz'))
        if len(fnlist) < 1:
            raise Exception(f'Found no file for fnlist={fnlist}, for dir={out_dir}')
        logger.debug(f'Find file: {fnlist}')

        basefn = os.path.basename(fnlist[0])
        tagname = basefn.replace('Meth_corr_plot_data_joined-', '')
        dsname = tagname[:tagname.find('_')]

        logger.info(f"Start report PCC in difference genomic regions based on file={fnlist[0]}, dsname={dsname}")
        correlation_report_on_regions(
            fnlist[0], bed_tuple_list=eval_genomic_context_tuple, dsname=dsname,
            runid=args.runid,
            outdir=out_dir)
        logger.debug(f"Memory report: {get_current_memory_usage()}")

        if args.report_no_join:
            # file like: Meth_corr_plot_data_bgtruth-HL60_RRBS_2Reps_NANOME-bsCov5-minToolCov3-baseFormat1.sorted.csv.gz
            fnlist = glob.glob(os.path.join(out_dir,
                                            f'Meth_corr_plot_data_bgtruth-*bsCov{args.min_bgtruth_cov}-minToolCov{args.toolcov_cutoff}*.csv.gz'))
            if len(fnlist) < 1:
                raise Exception(f'Found no file for fnlist={fnlist}, for dir={out_dir}')
            logger.debug(f'Find file: {fnlist}')

            basefn = os.path.basename(fnlist[0])
            tagname = basefn.replace('Meth_corr_plot_data_bgtruth-', '')
            dsname = tagname[:tagname.find('_')]

            logger.info(
                f"Start report no-joined sites PCC in difference genomic regions based on file={fnlist[0]}, dsname={dsname}")
            ## Due to large memory for no-joined, using non-mpi method
            correlation_report_on_regions(
                fnlist[0], bed_tuple_list=eval_genomic_context_tuple, dsname=dsname,
                runid=args.runid,
                outdir=out_dir,
                join_tag="no_join", mpi=False)
            logger.debug(f"Memory report: {get_current_memory_usage()}")

    if args.summary_coverage:
        logger.info("Start summarize CPG coverage at each regions")
        summary_cpgs_stats_results_table()
        logger.info(f"Memory report: {get_current_memory_usage()}")

    if args.plot:
        # plot fig5a of correlation plot
        command = f"set -x; plot_figure.py fig5a -i {outfn_joined_sorted} -o {out_dir}"
        subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")

    save_done_file(out_dir)
    logger.info(f"Memory report: {get_current_memory_usage()}")
    logger.info("### Site level correlation analysis DONE")
