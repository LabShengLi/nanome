#!/usr/bin/env python3

"""
Generate site-level methylation correlation results in nanome paper.
"""
import argparse
import subprocess

import pybedtools

from nanocompare.eval_common import *
from nanocompare.global_settings import get_tool_name, Top3ToolNameList, ToolNameList, save_done_file


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
    dataset = []
    joinedSet = None
    unionSet = set()

    ## create region bed files
    logger.info("Create region bed list firstly, take times.")
    region_fn_list = narrowCoordFileList[
                     1:] + cg_density_file_list + rep_file_list
    region_fn_list = narrowCoordFileList[
                     1:3] # only singleton/non-singleton
    # regionbed_list = []
    # for region_fn in region_fn_list:
    #     logger.info(f"Load file {region_fn}")
    #     basefn = os.path.basename(region_fn)
    #     tagname = location_filename_to_abbvname[basefn]
    #     region_bed = get_region_bed(region_fn)
    #     regionbed_list.append((region_fn, tagname, region_bed))

    # TODO: check bed script intersect only
    # outfn = os.path.join(pic_base_dir, basefn)
    # region_bed.saveas(outfn)
    # logger.info(f"save to {outfn}")
    region_bed_list = get_region_bed_pairs_list_mp(region_fn_list)

    retList = []
    logger.info("Start to study coverage now:")
    for toolname in loaded_callname_list:
        ## CpG sites set with cov >= cutoff(3)
        logger.info(f'\n\nStudy tool={toolname}')
        callSet = list(set(callresult_dict_cov3[toolname].keys()))

        if not joinedSet:
            joinedSet = set(callSet)
        else:
            joinedSet = joinedSet.intersection(callSet)
        unionSet = unionSet.union(callSet)
        row_dict = {
                    'Total CpG sites by Nanopore tool': call_cov1_cpg_sites[toolname],
                    f'Total CpG sites by tool cov>={minToolCovCutt}': len(callresult_dict_cov3[toolname])}
        row_dict.update({'Total calls by Nanopore reads': call_cov1_calls[toolname]})

        callBed = calldict2bed(callSet)

        # tagname -> set of intersected with bed region name
        intersect_set_dict={}
        # Add coverage of every regions by each tool here
        for (bedfn, tagname, region_bed) in tqdm(
                region_bed_list):  # calculate how overlap with Singletons, Non-Singletons, etc.
            # get_nsites_in_regions(callSet, bedfn, tagname)
            intersect_bed = intersect_bed_regions(callBed, region_bed, bedfn)
            intersect_set_dict[tagname] = set(bedtxt2dict(intersect_bed).keys())
            ret = {tagname: len(intersect_bed)}
            retList.append(ret)
        retDict = {}
        for e in retList:
            retDict.update(e)
        logger.info(f"retDict={retDict}")

        ## Sanity check
        sum_sing_nonsingle = retDict['Singletons'] + retDict['Non-singletons']
        total_sites = row_dict[f'Total CpG sites by tool cov>={minToolCovCutt}']
        logger.info(
            f"\n\nSanity check: sum_sing_nonsingle={sum_sing_nonsingle:,}; total_sites={total_sites:,}")

        if sum_sing_nonsingle != total_sites:
            logger.error(
                f"Sanity check for {toolname}, total_sites={total_sites:,}, sum_sing_nonsingle={sum_sing_nonsingle:,}, some non-singletons are not captered by bed file")
            retDict['Non-singletons'] = total_sites - retDict['Singletons']
            logger.error(f"Updated, retDict={retDict}")

            logger.info("Start to investigate problems......")
            toolSet = set(callSet)
            logger.info(f"toolSet={len(toolSet)}")

            subtractSingletonSet = toolSet - intersect_set_dict['Singletons']
            logger.info(f"subtractSingletonSet={len(subtractSingletonSet)}")

            subtractNonsingletonSet = subtractSingletonSet - intersect_set_dict['Non-singletons']
            logger.info(f"subtractNonsingletonSet={len(subtractNonsingletonSet)}")

            logger.info(subtractNonsingletonSet)

            bed_not_covered = calldict2bed(subtractNonsingletonSet)
            outfn = os.path.join(out_dir, f'{args.runid}-not_covered_sites.bed')
            bed_not_covered.saveas(outfn)

            logger.info("### DONE")
        row_dict.update(retDict)
        dataset.append(row_dict)

    df = pd.DataFrame(dataset, index=loaded_callname_list)
    logger.info(df)
    outfn = os.path.join(out_dir,
                         f'{RunPrefix}-summary-bgtruth-tools-bsCov{bgtruthCutt}-minCov{minToolCovCutt}.table.s10.xlsx')
    df.to_excel(outfn)
    logger.info(f'save to {outfn}\n')


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Site level correlation analysis')
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by spaces',
                        required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name1>;<file-name1>",
                        required=True)
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--beddir', type=str, help="base dir for bed files",
                        default=None)  # need perform performance evaluation before, then get concordant, etc. bed files, like '/projects/li-lab/yang/results/2021-04-01'
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=',')
    parser.add_argument('--processors', type=int, help="running processors", default=16)
    parser.add_argument('--min-bgtruth-cov', type=int, help="cutoff for coverage in bg-truth", default=5)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff for coverage in nanopore tools", default=3)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--gen-venn', help="generate venn data", action='store_true')
    parser.add_argument('--summary-coverage', help="generate summary table for coverage at each region",
                        action='store_true')
    parser.add_argument('--enable-cache', action='store_true')
    parser.add_argument('--using-cache', action='store_true')
    parser.add_argument('--plot', help="plot the correlation figure", action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    set_log_debug_level()

    ## Set tmp dir for bedtools
    bedtool_tmp_dir = "/fastscratch/liuya/nanocompare/bedtools_tmp"
    os.makedirs(bedtool_tmp_dir, exist_ok=True)
    pybedtools.helpers.set_tempdir(bedtool_tmp_dir)

    args = parse_arguments()

    # cache function same with read level
    enable_cache = args.enable_cache
    using_cache = args.using_cache

    # if enable_cache:
    #     os.makedirs(cache_dir, exist_ok=True)

    # runid is always like 'MethCorr-K562_WGBS_2Reps', remove first word as RunPrefix like K562_WGBS_2Reps
    RunPrefix = args.runid.replace('Sanity_Nanopolish-', '')

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
    # add_logging_file(os.path.join(out_dir, 'run-results.log'))

    logger.debug(args)

    logger.info(f'\n\n####################\n\n')

    callfn_dict = defaultdict()  # callname -> filename

    # callname -> [call0, call1], call0 is no-filter results, call1 is filter by cutoff, and convert to [meth-freq, meth-cov] results.
    callresult_dict_cov1 = defaultdict()
    call_cov1_cpg_sites = defaultdict(int)
    call_cov1_calls = defaultdict(int)
    callresult_dict_cov3 = defaultdict()
    loaded_callname_list = []

    for callstr in args.calls:
        callencode, callfn = callstr.split(':')

        if len(callfn) == 0:
            continue

        callname = get_tool_name(callencode)

        if callname not in ['Nanopolish', 'Megalodon']:
            continue
        callfn_dict[callname] = callfn
        loaded_callname_list.append(callname)

        callresult_dict_cov1[callname] = import_call(callfn, callencode, baseFormat=baseFormat,
                                                     enable_cache=enable_cache, using_cache=using_cache,
                                                     include_score=False, siteLevel=True)

        # Stats the total cpgs and calls for each calls
        cnt_calls = 0
        for cpg in callresult_dict_cov1[callname]:
            cnt_calls += len(callresult_dict_cov1[callname][cpg])
        call_cov1_calls[callname] = cnt_calls
        call_cov1_cpg_sites[callname] = len(callresult_dict_cov1[callname])

    logger.debug(loaded_callname_list)

    logger.info(f'Start apply cutoff={minToolCovCutt} to methylation calls, take time')

    # Cutoff of read cov >= 1 or 3, 5 for nanopore tools
    for callname in loaded_callname_list:
        callresult_dict_cov3[callname] = readLevelToSiteLevelWithCov(callresult_dict_cov1[callname],
                                                                     minCov=minToolCovCutt, toolname=callname)
        ## Destroy cov1 for memory saving
        del callresult_dict_cov1[callname]

    logger.info(f'\n\n####################\n\n')

    logger.info(f"Start getting intersection (all joined) sites by tools and bgtruth")

    if args.summary_coverage:
        logger.info("Start summarize coverage at each regions")
        summary_cpgs_stats_results_table()

    save_done_file(out_dir)
    logger.info("### Site level correlation analysis DONE")
