#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : site_level_eval.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

"""
Generate site-level methylation correlation results in nanome paper.
"""

import argparse

import pybedtools

from nanocompare.eval_common import *
from nanocompare.global_settings import get_tool_name, Top3ToolNameList, ToolNameList, save_done_file, \
    narrowCoordNameList, cg_density_coord_name_list, rep_coord_name_list, nanome_version


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
    bgtruthCpGs = set(list(bgTruth.keys()))
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

        # Add coverage of every regions by each tool here
        for (bedfn, tagname, region_bed) in tqdm(
                region_bed_list):  # calculate how overlap with Singletons, Non-Singletons, etc.
            if not args.large_mem:  # load in demand
                region_bed = get_region_bed_tuple(bedfn, enable_base_detection_bedfile=not args.disable_bed_check)[2]

            if region_bed is None:
                logger.warning(f"region name={tagname} is not found")
                continue
            intersect_bed = intersect_bed_regions(callBed, region_bed, bedfn)
            ret = {tagname: len(intersect_bed)}
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
        sum_sing_nonsingle = retDict['Singletons'] + retDict['Non-singletons']
        if 'CG_20' in retDict and 'CG_40' in retDict and 'CG_60' in retDict and 'CG_80' in retDict and 'CG_100' in retDict:
            sum_cg = retDict['CG_20'] + retDict['CG_40'] + retDict['CG_60'] + retDict['CG_80'] + retDict['CG_100']
        else:
            sum_cg = 0
        total_sites = row_dict[f'Total CpG sites by tool cov>={minToolCovCutt}']
        logger.debug(
            f"\n\nSanity check: sum_sing_nonsingle={sum_sing_nonsingle:,}; sum_cg={sum_cg:,}; total={total_sites:,}")

        if sum_sing_nonsingle != total_sites:
            logger.error(
                f"Sanity check for {toolname}, total_sites={total_sites:,}, sum_sing_nonsingle={sum_sing_nonsingle:,}, some non-singletons are not captered by bed file")
            retDict['Non-singletons'] = total_sites - retDict['Singletons']
            logger.error(f"Updated, retDict={retDict}")
        row_dict.update(retDict)
        dataset.append(row_dict)
        logger.debug(f'BG-Truth join with {toolname} get {len(toolOverlapBGTruthCpGs):,} CpGs')
        # outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-bgtruth-{name1}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
        # save_keys_to_bed(overlapCpGs, outfn)

    # add additional rows for Joined count
    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(joinedSet)}
    dataset.append(new_row_dict)

    # add additional row for Unioned count
    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(unionSet)}
    dataset.append(new_row_dict)

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

    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top3JointSet)}
    dataset.append(new_row_dict)

    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top3UnionSet)}
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

    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top4JointSet)}
    dataset.append(new_row_dict)

    new_row_dict = {f'Total CpG sites by tool cov>={minToolCovCutt}': len(top4UnionSet)}
    dataset.append(new_row_dict)

    df = pd.DataFrame(dataset, index=loaded_callname_list +
                                     ['Joined', 'Union', 'TOP3 Joined', 'TOP3 Union', 'TOP4 Joined', 'TOP4 Union'])

    logger.debug(df)

    outfn = os.path.join(out_dir,
                         f'{RunPrefix}-summary-bgtruth-tools-bsCov{bgtruthCutt}-minCov{minToolCovCutt}.table.s10.xlsx')
    df.to_excel(outfn)
    logger.debug(f'save to {outfn}\n')


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
            if cpg in callresult_dict[name]:  # if cpg is in tool results
                row_list.extend([f'{callresult_dict[name][cpg][0]:.3f}', f'{callresult_dict[name][cpg][1]}'])
            else:  # if cpg key is not exist, we use NA as ''
                row_list.extend(['', ''])
        outfile.write(sep.join(row_list))
        outfile.write("\n")
    outfile.close()
    logger.debug(f"save to {outfn}\n")


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(prog='site_level_eval (NANOME)',
                                     description='Site-level correlation analysis in nanome paper')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v{nanome_version}')
    parser.add_argument('--dsname', type=str, help="dataset name", required=True)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name> seperated by spaces',
                        required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name1>;<file-name1>",
                        required=True)
    parser.add_argument('--beddir', type=str, help="base dir for concordant/discordant bed files",
                        default=None)
    parser.add_argument('--min-bgtruth-cov', type=int, help="cutoff for coverage in bg-truth", default=5)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff for coverage in nanopore tools", default=3)
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=',')
    parser.add_argument('--processors', type=int, help="running processors", default=1)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    parser.add_argument('--gen-venn', help="generate CpGs for venn data analysis", action='store_true')
    parser.add_argument('--summary-coverage', help="generate summary table for coverage at each region",
                        action='store_true')
    parser.add_argument('--region-coe-report', help="report table for PCC value at each region",
                        action='store_true')
    parser.add_argument('--enable-cache', help="if enable cache functions", action='store_true')
    parser.add_argument('--using-cache', help="if use cache files", action='store_true')
    parser.add_argument('--plot', help="plot the correlation matrix figure", action='store_true')
    parser.add_argument('--bedtools-tmp', type=str, help='bedtools temp dir', default=temp_dir)
    parser.add_argument('--genome-annotation', type=str, help='genome annotation dir',
                        default=os.path.join(data_base_dir, 'genome-annotation'))
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
    os.makedirs(args.bedtools_tmp, exist_ok=True)
    pybedtools.helpers.set_tempdir(args.bedtools_tmp)

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
    logger.debug(f'\n\n####################\n\n')

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
                                  using_cache=using_cache, enable_cache=enable_cache)
        bgTruthList.append(bgTruth1)

    # Combine one/two replicates, using cutoff=1 or 5
    bgTruth = combineBGTruthList(bgTruthList, covCutoff=bgtruthCutt)

    logger.info(f'Combined BS-seq data (cov>={bgtruthCutt}), all methylation level sites={len(bgTruth):,}')

    logger.debug(f'\n\n####################\n\n')

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
        callfn_dict[callname] = callfn

        # We do now allow import DeepMod.C for site level evaluation, in current version
        if callencode == 'DeepMod.C':
            raise Exception(
                f'{callencode} is not allowed for site level evaluation, please use DeepMod.Cluster file here')

        loaded_callname_list.append(callname)

        # For site level evaluation, only need (freq, cov) results, no score needed. Especially for DeepMod, we must import as freq and cov format from DeepMod.Cluster encode
        # Do not filter bgtruth, because we use later for overlapping (without bg-truth)
        callresult_dict_cov1[callname] = import_call(callfn, callencode, baseFormat=baseFormat,
                                                     enable_cache=enable_cache, using_cache=using_cache,
                                                     include_score=False, siteLevel=True)

        # Stats the total cpgs and calls for each calls
        cnt_calls = 0
        for cpg in callresult_dict_cov1[callname]:
            cnt_calls += len(callresult_dict_cov1[callname][cpg])
        call_cov1_calls[callname] = cnt_calls
        call_cov1_cpg_sites[callname] = len(callresult_dict_cov1[callname])

    logger.info(f"Import calls from tools done for toollist={list(call_cov1_cpg_sites.keys())}")
    logger.info(f"Memory report: {get_current_memory_usage()}")

    logger.debug(f'Start apply cutoff={minToolCovCutt} to methylation calls, take time')
    # Cutoff of read cov >= 1 or 3, 5 for nanopore tools
    for callname in loaded_callname_list:
        callresult_dict_cov3[callname] = readLevelToSiteLevelWithCov(callresult_dict_cov1[callname],
                                                                     minCov=minToolCovCutt, toolname=callname)
        ## Destroy cov1 for memory saving
        del callresult_dict_cov1[callname]

    logger.debug(f'\n\n####################\n\n')

    top3_cpg_set_dict = defaultdict()
    for callname in Top3ToolNameList:
        top3_cpg_set_dict[callname] = set(callresult_dict_cov3[callname].keys())

    if args.gen_venn:
        logger.info('CPG overlapping analysis')
        logger.debug(
            f"Start gen venn data for each tool (cov>={minToolCovCutt}) and BS-seq (cov>={args.min_bgtruth_cov})")

        # Generate all tools and bsseq covered cpgs files for set evaluation
        logger.debug("We generate sets file for each tool and bg-truth")
        venn_outdir = os.path.join(out_dir, 'venn_data')
        os.makedirs(venn_outdir, exist_ok=True)

        bg_cpgs = bgTruth.keys()
        outfn = os.path.join(venn_outdir, f'{args.dsname}.bgtruth.cpg.sites.cov{args.min_bgtruth_cov}.setsfile.txt.gz')
        ontcalls_to_setsfile_for_venn_analysis(bg_cpgs, outfn)
        for callname in ToolNameList:
            if callname not in callresult_dict_cov3:
                continue
            call_keys = callresult_dict_cov3[callname].keys()
            outfn = os.path.join(venn_outdir,
                                 f'{args.dsname}.{callname}.cpg.sites.cov{args.toolcov_cutoff}.setsfile.txt.gz')
            ontcalls_to_setsfile_for_venn_analysis(call_keys, outfn)
        logger.debug(f'\n\n####################\n\n')

    logger.debug(f"Start getting intersection (all joined) sites by tools and bgtruth")
    coveredCpGs = set(list(bgTruth.keys()))
    coveredCpGs001 = set(list(bgTruth.keys()))

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
                         f'{args.dsname}.tools.cov{args.toolcov_cutoff}.join.with.bsseq.cov{args.min_bgtruth_cov}.site.level.report.csv')
    df.to_csv(outfn)

    # Output sites report for all tools and BS-seq
    logger.info(f"Joined {len(coveredCpGs):,} CpGs are covered by all tools (cov >= {args.toolcov_cutoff}) and BS-seq (cov >= {args.min_bgtruth_cov})")
    logger.debug('Output data of meth-freq and coverage on joined CpG sites as datasets for correlation analysis')
    outfn_joined = os.path.join(out_dir,
                                f"Meth_corr_plot_data_joined-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.csv.gz")
    save_meth_corr_data(callresult_dict_cov3, bgTruth, coveredCpGs, outfn_joined)
    outfn_joined_sorted = os.path.join(out_dir,
                                       f"Meth_corr_plot_data_joined-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.sorted.csv.gz")
    sort_bed_file(infn=outfn_joined, outfn=outfn_joined_sorted, has_header=True)
    os.remove(outfn_joined)

    logger.debug(
        'Output data of meth-freq and coverage on bgTruth related CpG sites as datasets for correlation analysis')
    outfn_bgtruth = os.path.join(out_dir,
                                 f"Meth_corr_plot_data_bgtruth-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.csv.gz")
    save_meth_corr_data(callresult_dict_cov3, bgTruth, set(list(bgTruth.keys())), outfn_bgtruth)
    outfn_bgtruth_sroted = os.path.join(out_dir,
                                        f"Meth_corr_plot_data_bgtruth-{RunPrefix}-bsCov{bgtruthCutt}-minToolCov{minToolCovCutt}-baseFormat{baseFormat}.sorted.csv.gz")
    sort_bed_file(infn=outfn_bgtruth, outfn=outfn_bgtruth_sroted, has_header=True)
    os.remove(outfn_bgtruth)

    # Report correlation matrix here
    df = pd.read_csv(outfn_joined_sorted, sep=sep)
    df = df.filter(regex='_freq$', axis=1)
    cordf = df.corr()
    logger.debug(f'Correlation matrix is:\n{cordf}')
    corr_outfn = os.path.join(out_dir,
                              f'Meth_corr_plot_data-{RunPrefix}-correlation-matrix-toolcov{minToolCovCutt}-bsseqcov{bgtruthCutt}.xlsx')
    cordf.to_excel(corr_outfn)

    logger.debug(f'\n\n####################\n\n')

    if args.region_coe_report or args.summary_coverage:
        ## load region bed list
        logger.debug("Create region bed list firstly, take times......")

        # Evaluated all region filename lists, bed objects
        # assume all files are located in args.genome_annotation dir
        regions_full_filepath = [os.path.join(args.genome_annotation, cofn) for cofn in narrowCoordNameList[1:]] + \
                                [os.path.join(args.genome_annotation, cofn) for cofn in cg_density_coord_name_list] + \
                                [os.path.join(args.genome_annotation, cofn) for cofn in rep_coord_name_list]
        # logger.debug(f"Evaluated regions: {regions_full_filepath}")

        if args.large_mem:  # load all in memory
            region_bed_list = get_region_bed_pairs_list_mp(regions_full_filepath, processors=args.processors, enable_base_detection_bedfile=not args.disable_bed_check)
        else:  # load bed coord later
            region_bed_list = [(infn, map_region_fn_to_name(infn), None,)
                               for infn in regions_full_filepath]

        if args.beddir:  # add concordant and discordant region coverage if needed
            logger.debug(f'We use Concordant and Discordant BED file at basedir={args.beddir}')
            concordantFileName = find_bed_filename(basedir=args.beddir,
                                                   pattern=f'{args.dsname}*hg38_nonsingletons*.concordant.bed.gz')
            concordant_bed = get_region_bed(concordantFileName)
            discordantFileName = find_bed_filename(basedir=args.beddir,
                                                   pattern=f'{args.dsname}*hg38_nonsingletons*.discordant.bed.gz')
            discordant_bed = get_region_bed(discordantFileName)
        else:
            concordant_bed = None
            discordant_bed = None

    if args.region_coe_report:
        eval_genomic_context_tuple = [('x.x.Genome-wide', 'Genome-wide', None)] + region_bed_list
        if concordant_bed is not None:
            eval_genomic_context_tuple += [(concordantFileName, 'Concordant', concordant_bed,)]
        if discordant_bed is not None:
            eval_genomic_context_tuple += [(discordantFileName, 'Discordant', discordant_bed,)]

        logger.debug(f"Evaluated on regions: {eval_genomic_context_tuple}")

        # file like: Meth_corr_plot_data_joined-TestData_RRBS_2Reps-bsCov1-minToolCov1-baseFormat1.sorted.csv.gz
        fnlist = glob.glob(os.path.join(out_dir, f'Meth_corr_plot_data_joined-*bsCov{args.min_bgtruth_cov}-minToolCov{args.toolcov_cutoff}*.csv.gz'))
        if len(fnlist) < 1:
            raise Exception(f'Found no file for fnlist={fnlist}, for dir={out_dir}')
        logger.debug(f'Find file: {fnlist}')

        basefn = os.path.basename(fnlist[0])
        tagname = basefn.replace('Meth_corr_plot_data_joined-', '')
        dsname = tagname[:tagname.find('_')]

        logger.info(f"Start report PCC in difference genomic regions based on file={fnlist[0]}, dsname={dsname}")
        correlation_report_on_regions(fnlist[0], bed_tuple_list=eval_genomic_context_tuple, dsname=dsname,
                                      outdir=out_dir, large_mem=args.large_mem,
                                      enable_base_detection_bedfile=not args.disable_bed_check)

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
