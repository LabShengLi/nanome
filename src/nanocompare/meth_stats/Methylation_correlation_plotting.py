#!/home/liuya/anaconda3/envs/nmf/bin/python

"""
Generate methy correlation plotting input data as a tsv

Sample usage:
    python /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/meth_stats/Methylation_correlation_plotting.py --calls DeepSignal:/fastscratch/liuya/nanocompare/K562-Runs/K562-DeepSignal-N50/K562-DeepSignal-N50-meth-call/K562.deepsignal.call_mods.combine.tsv Tombo:/fastscratch/liuya/nanocompare/K562-Runs/K562-Tombo-N1/K562-Tombo-N1-meth-call/K562.tombo.perReadsStatsOnlyCG.combine.tsv Nanopolish:/fastscratch/liuya/nanocompare/K562-Runs/K562-Nanopolish-N50/K562-Nanopolish-N50-meth-call/K562.nanopolish.methylation_calls.combine.tsv DeepMod:/fastscratch/liuya/nanocompare/deepmod-read-level1.tsv --bgtruth bed:/projects/li-lab/yang/workspace/nano-compare/data/bgtruth-data/K562_joined.bed.gz --runid K562_WGBS_Joined_NewRuns

All usedful functions are located in nanocompare.meth_stats.meth_stats_common
"""
import argparse

from nanocompare.global_config import *
from nanocompare.meth_stats.meth_stats_common import *


# nanocompare_prj = "/projects/li-lab/yang/workspace/nano-compare/src"
# sys.path.append(nanocompare_prj)


def scatter_plot_cov():
    logger.debug(f"Study scatter plot of coverage for Tombo with other tool with no cutoff")

    scatter_analysis_cov(Tombo_calls0, Nanopolish_calls0, outdir=out_dir, RunPrefix=RunPrefix, tool1_name='Tombo', tool2_name='Nanopolish')
    scatter_analysis_cov(Tombo_calls0, DeepSignal_calls0, outdir=out_dir, RunPrefix=RunPrefix, tool1_name='Tombo', tool2_name='Deepsignal')

    scatter_analysis_cov(DeepMod_calls0, Nanopolish_calls0, outdir=out_dir, RunPrefix=RunPrefix, tool1_name='DeepMod', tool2_name='Nanopolish')
    scatter_analysis_cov(DeepMod_calls0, DeepSignal_calls0, outdir=out_dir, RunPrefix=RunPrefix, tool1_name='DeepMod', tool2_name='Deepsignal')


def study_intersection():
    logger.debug(f"Study set intersection of each tool with bgtruth")
    dataset = []
    bgtruthCpGs = set(list(bgTruth.keys()))
    for call1, name1, call1_before_cutoff in zip(all_calls, name_calls, all_calls_before_cutoff):
        if name1 == "DeepMod_cluster":
            continue
        if call1 is None:
            # dataset.append({'bgtruth': len(bgtruthCpGs), 'tool': np.nan, f'tool-covcut{minToolCovCutt}': np.nan, 'joined': np.nan})
            continue
        overlapCpGs = bgtruthCpGs.intersection(set(list(call1.keys())))
        ret = {f'CpG sites in BG-Truth cov>={bgtruthCutt}': len(bgtruthCpGs), 'Total CpG sites by Nanopore tool': len(set(list(call1_before_cutoff))), f'Total CpG sites by cov>={minToolCovCutt}': len(set(list(call1.keys()))), 'Joined CpG sites with BG-Truth': len(overlapCpGs)}

        cnt_calls = 0

        if name1 != "DeepMod_cluster":
            for cpg in call1_before_cutoff:
                cnt_calls += len(call1_before_cutoff[cpg])
        else:
            for cpg in call1_before_cutoff:
                cnt_calls += call1_before_cutoff[cpg][1]
        ret.update({'Total calls by Nanopore reads': cnt_calls})
        dataset.append(ret)

        logger.info(f'BG-Truth join with {name1} get {len(overlapCpGs):,} CpGs')
        # outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-bgtruth-{name1}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
        # save_keys_to_bed(overlapCpGs, outfn)

    df = pd.DataFrame(dataset, index=name_calls[:-1])
    df = df.iloc[:, [4, 0, 1, 2, 3]]
    outfn = os.path.join(out_dir, f'{RunPrefix}-summary-bgtruth-tools-bsCov{bgtruthCutt}-minCov{minToolCovCutt}.csv')
    df.to_csv(outfn)

    return

    logger.debug(f"Study set intersection of deepsignal, nanopolish with bgtruth")
    dataset = []
    overlapCpGs = set(list(bgTruth.keys()))
    for call1, name1 in zip(all_calls, name_calls):
        if call1 is None:
            continue
        if name1 not in ['DeepSignal', 'Nanopolish']:
            continue
        overlapCpGs = overlapCpGs.intersection(set(list(call1.keys())))
    # outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-deepsignal-nanopolish-bgtruth-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
    # save_keys_to_bed(overlapCpGs, outfn)
    logger.info(f'Reporting deepsignal, nanopolish with bgtruth, joined results = {len(overlapCpGs)}')

    logger.debug(f"Study set intersection of deepsignal, nanopolish, deepmod and tombo")
    dataset = []
    overlapCpGs = None
    for call1, name1 in zip(all_calls, name_calls):
        if call1 is None:
            continue
        if name1 not in ['DeepSignal', 'Nanopolish', 'DeepMod', 'Tombo']:
            continue
        if overlapCpGs is None:
            overlapCpGs = set(list(call1.keys()))
        else:
            overlapCpGs = overlapCpGs.intersection(set(list(call1.keys())))
    # outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-deepsignal-nanopolish-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
    # save_keys_to_bed(overlapCpGs, outfn)
    logger.info(f'Reporting deepsignal, nanopolish,deepmod,and tombo joined results = {len(overlapCpGs)}')

    logger.debug(f"Next, study set intersection of deepsignal, nanopolish but not in Tombo")
    newOverlapCpGs = overlapCpGs - Tombo_calls.keys()

    # outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-deepsignal-nanopolish-notombo-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
    # save_keys_to_bed(newOverlapCpGs, outfn)
    logger.info(f'Reporting deepsignal, nanopolish but not in Tombo results = {len(newOverlapCpGs)}')


def study_high_cov_tool1_low_cov_tool2():
    logger.debug(f'Study and find high cov in DeepSignal and Nanopolish, but low in Tombo (before cutoff)')

    df = get_high_cov_call1_low_cov_call2_df(Nanopolish_calls0, Tombo_calls0, low_cutoff=minToolCovCutt)
    # logger.debug(df)
    outfn = os.path.join(out_dir, f'{RunPrefix}-high-cov-nanopolish-low-cov-tombo-baseCount{baseFormat}.bed')
    df.to_csv(outfn, sep='\t', index=False)

    outfn = os.path.join(out_dir, f'{RunPrefix}-high-cov-nanopolish-low-cov-tombo-baseCount{baseFormat}.pkl')
    df.to_pickle(outfn)

    logger.debug(f'Study DeepSignal cov>={minToolCovCutt} , but low in Tombo cov < {minToolCovCutt}')

    df = get_high_cov_call1_low_cov_call2_df(DeepSignal_calls0, Tombo_calls0, low_cutoff=minToolCovCutt, call1_name='deepsignal')
    logger.debug(df)
    outfn = os.path.join(out_dir, f'{RunPrefix}-high-cov-deepsignal-low-cov-tombo-baseCount{baseFormat}.bed')
    df.to_csv(outfn, sep='\t', index=False)

    outfn = os.path.join(out_dir, f'{RunPrefix}-high-cov-deepsignal-low-cov-tombo-baseCount{baseFormat}.pkl')
    df.to_pickle(outfn)


def output_bed_of_tools_bgtruth():
    logger.debug(f"Output bed files of each tool, and bgtruth")

    # Note bgTruth in format of {'chr\t123\t123\n':[0.56, 15], etc.}
    outfn = os.path.join(out_dir, f'{RunPrefix}-meth-cov-bgtruth-baseCount{baseFormat}.bed')
    save_call_or_bgtruth_to_bed(bgTruth, outfn, callBaseFormat=baseFormat)

    # Note each tool call filtered results is in format of {'chr\t123\t123\n':[0.56, 15], etc.}
    for call1, name1, call2 in zip(all_calls, name_calls, all_calls_before_cutoff):
        outfn = os.path.join(out_dir, f'{RunPrefix}-meth-cov{minToolCovCutt}-tool-{name1}-baseCount{baseFormat}.bed')
        save_call_or_bgtruth_to_bed(call1, outfn, callBaseFormat=baseFormat)

        outfn = os.path.join(out_dir, f'{RunPrefix}-meth-covcutoff{minToolCovCutt}-tool-{name1}-baseCount{baseFormat}.myCpG.txt')
        save_call_to_methykit_txt(call1, outfn, callBaseFormat=baseFormat)

        if 'DeepMod_cluster' == name1:
            continue
        outfn = os.path.join(out_dir, f'{RunPrefix}-meth-nocovcutoff-tool-{name1}-baseCount{baseFormat}.myCpG.txt')
        save_call_to_methykit_txt(call2, outfn, callBaseFormat=baseFormat, is_cov=False)


def summary_cpgs_joined_results_table():
    """
    Study and summary each tool joined with bg-truth results, make table as dataframe
    :return:
    """
    logger.debug(f"Study set intersection of each tool with bgtruth")
    dataset = []
    bgtruthCpGs = set(list(bgTruth.keys()))
    for name in callname_list:
        overlapCpGs = bgtruthCpGs.intersection(set(list(callresult_dict[name][1].keys())))
        ret = {f'CpG sites in BG-Truth cov>={bgtruthCutt}': len(bgtruthCpGs), 'Total CpG sites by Nanopore tool': len(callresult_dict[name][0]), f'Total CpG sites by cov>={minToolCovCutt}': len(callresult_dict[name][1]), 'Joined CpG sites with BG-Truth': len(overlapCpGs)}

        cnt_calls = 0

        if name != "DeepMod_cluster":
            for cpg in callresult_dict[name][0]:
                cnt_calls += len(callresult_dict[name][0][cpg])
        else:
            raise Exception("DeepMod_cluster need to be deal with")
        ret.update({'Total calls by Nanopore reads': cnt_calls})
        dataset.append(ret)

        logger.info(f'BG-Truth join with {name} get {len(overlapCpGs):,} CpGs')
        # outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-bgtruth-{name1}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
        # save_keys_to_bed(overlapCpGs, outfn)

    df = pd.DataFrame(dataset, index=callname_list)
    df = df.iloc[:, [4, 0, 1, 2, 3]]
    outfn = os.path.join(out_dir, f'{RunPrefix}-summary-bgtruth-tools-bsCov{bgtruthCutt}-minCov{minToolCovCutt}.csv')
    df.to_csv(outfn, sep=args.sep)
    logger.info(f'save to {outfn}\n')


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Correlation data gen task')
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name>', required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name>", required=True)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--sep', type=str, help="seperator for output csv file", default=',')
    parser.add_argument('--processors', type=int, help="running processors", default=8)
    parser.add_argument('--bgtruthcov-cutoff', type=int, help="cutoff of coverage in bg-truth", default=5)
    parser.add_argument('--toolcov-cutoff', type=int, help="cutoff of coverage in nanopore calls", default=4)
    parser.add_argument('--baseFormat', type=int, help="cutoff of coverage in nanopore calls", default=0)
    parser.add_argument('-o', type=str, help="output dir", default=pic_base_dir)
    return parser.parse_args()


if __name__ == '__main__':
    set_log_debug_level()

    args = parse_arguments()
    logger.debug(args)

    RunPrefix = args.runid  # "K562_WGBS_rep_ENCFF721JMB"
    # tool coverage cutoff 4
    minToolCovCutt = args.toolcov_cutoff

    # bgtruth coverage cutoff 10
    bgtruthCutt = args.bgtruthcov_cutoff

    # load into program format 0-base or 1-base
    baseFormat = args.baseFormat

    sep = args.sep

    out_dir = os.path.join(args.o, RunPrefix)
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Output to dir:{out_dir}')

    callfn_dict = defaultdict()  # callname -> filename
    callresult_dict = defaultdict()
    callname_list = []

    # with Pool(processes=args.processors) as pool:
    for callstr in args.calls:
        callname, callfn = callstr.split(':')
        callfn_dict[callname] = callfn
        callname_list.append(callname)
        callresult_dict[callname] = [import_call(callfn, callname, baseFormat=baseFormat)]
        # callresult_dict[callname] = pool.apply_async(import_call, (callfn, callname,))

    encode, fn = args.bgtruth.split(':')

    logger.debug(f'fn={fn}, encode={encode}')
    bgTruth = import_bgtruth(fn, encode, cov=bgtruthCutt, baseFormat=baseFormat)
    # bgTruth = pool.apply_async(import_bgtruth, (fn, encode,))
    # pool.close()
    # pool.join()
    # logger.info("Join processes of import")

    for callname in callname_list:
        # callresult_dict[callname] = [callresult_dict[callname].get()]
        callresult_dict[callname].append(coverageFiltering(callresult_dict[callname][0], minCov=minToolCovCutt, toolname=callname))
    # bgTruth = bgTruth.get()

    logger.debug(callfn_dict)
    logger.debug(callname_list)

    logger.info(f"Start set intersection with all tools joined together with bgtruth")

    coveredCpGs = set(list(bgTruth.keys()))
    for name in callname_list:
        coveredCpGs = coveredCpGs.intersection(set(list(callresult_dict[name][1].keys())))
        logger.info(f'Join {name} get {len(coveredCpGs)} CpGs')
    logger.info(f"Reporting {len(coveredCpGs)} CpGs are covered by all tools and bgtruth")

    logger.info('Output data of coverage and meth-freq on joined CpG sites for correlation analysis')
    outfn = os.path.join(out_dir, f"Meth_corr_plot_data-{RunPrefix}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.csv")

    outfile = open(outfn, 'w')

    # outfile.write("chr\tstart\tend\tBSseq_freq\tBSseq_cov\tstrand")

    header_list = ['chr', 'start', 'end', 'BGTruth_freq', 'BGTruth_cov', 'strand']

    for name in callname_list:
        header_list.extend([f'{name}_freq', f'{name}_cov'])
        # outfile.write(f'{sep}'.join(f'{sep}{name}_freq', f'{name}_cov'))
        # outfile.write(f"\t{name}_freq\t{name}_cov")
    outfile.write(sep.join(header_list))
    outfile.write("\n")

    for cpg in coveredCpGs:
        row_list = [cpg[0], str(cpg[1]), str(cpg[1] + 1), f'{bgTruth[cpg][0]:.3f}', str(bgTruth[cpg][1]), cpg[2]]
        # outfile.write(f'{sep}'.join(cpg[0], cpg[1], cpg[1] + 1, bgTruth[cpg][0], bgTruth[cpg][1], cpg[2]))
        # outfile.write(f"{cpg[0]}\t{cpg[1]}\t{cpg[1] + 1}\t{bgTruth[cpg][0]}\t{bgTruth[cpg][1]}\t{cpg[2]}")

        for name in callname_list:
            row_list.extend([f'{callresult_dict[name][1][cpg][0]:.3f}', f'{callresult_dict[name][1][cpg][1]}'])
        # outfile.write(f'{sep}'.join(f'{sep}{callresult_dict[name][1][cpg][0]}', f'{callresult_dict[name][1][cpg][1]}'))
        # outfile.write(f"\t{callresult_dict[name][1][cpg][0]}\t{callresult_dict[name][1][cpg][1]}")
        outfile.write(sep.join(row_list))
        outfile.write("\n")
    outfile.close()
    logger.info(f"save to {outfn}\n")

    summary_cpgs_joined_results_table()

    logger.info("Methylation correlation plotting data generation program finished.")

    sys.exit(0)

    is_scatter_plot = False
    output_bed = False
    study_high_low_cov = False

    if study_high_low_cov:
        study_high_cov_tool1_low_cov_tool2()

    if output_bed:
        output_bed_of_tools_bgtruth()

    if is_scatter_plot:
        scatter_plot_cov()

    study_intersection()
