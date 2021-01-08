#!/home/liuya/anaconda3/envs/nmf/bin/python

"""
Generate methy correlation plotting input data as a tsv

Sample usage:
    <pyfilename> $DeepSignal_calls $Tombo_calls $Nanopolish_calls \
			$DeepMod_calls $DeepMod_cluster_calls $bgTruth $parser $RunPrefix

All usedful functions are located in nanocompare.meth_stats.meth_stats_common
"""
import argparse
import sys
from multiprocessing import Pool

nanocompare_prj = "/projects/li-lab/yang/workspace/nano-compare/src"
sys.path.append(nanocompare_prj)

from nanocompare.meth_stats.meth_stats_common import *

from global_config import *


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


def parse_arguments():
    """
    usage: volume_calculation.py [-h] [-n N] [-t T] [--show]
                             [--input INPUT [INPUT ...]] [--output OUTPUT]
                             [--dcm] [--single-scan] [--lu-score LU_SCORE]
                             [--le-score LE_SCORE]
                             cmd

    Volume calculation for lung and lesion

    positional arguments:
      cmd                   name of command: compute, combine, or gen-pixel-info

    optional arguments:
      -h, --help            show this help message and exit
      -n N                  the total number of tasks (1-27)
      -t T                  the current task id (1-N)
      --show                show prediction images if using this switch
      --input INPUT [INPUT ...]
                            the input dir that contains scanid of pic/dcm files
      --output OUTPUT       the input pic dir
      --dcm                 folders are scanid that containing DCM files if using
                            this switch
      --single-scan         folders are directly the scanid folder if using this
                            switch
      --lu-score LU_SCORE   the lung field detection score
      --le-score LE_SCORE   the lesion field detection score
    :return:
    """
    parser = argparse.ArgumentParser(description='Correlation data gen task')
    # parser.add_argument('-n', type=int, help="the total number of tasks (1-27)", default=1)
    parser.add_argument('--calls', nargs='+', help='all ONT call results <tool-name>:<file-name>', required=True)
    parser.add_argument('--bgtruth', type=str, help="background truth file <encode-type>:<file-name>", required=True)
    parser.add_argument('--runid', type=str, help="running prefix", required=True)
    parser.add_argument('--processors', type=int, help="running processors", default=8)
    # parser.add_argument('-o', type=str, help="running prefix", default=pic_base_dir)

    # parser.add_argument('-n', type=int, help="the total number of tasks (1-27)", default=1)
    # parser.add_argument('-t', type=int, help="the current task id (1-N)", default=1)
    # parser.add_argument('-i', type=str, help="input file", default=None)
    # parser.add_argument('-o', type=str, help="output dir or file", default=pic_base_dir)
    # parser.add_argument('--o2', type=str, help="second output dir or file", default=None)
    # parser.add_argument('--ibam', type=str, help="input bam/sam file", default=None)
    # parser.add_argument('--basecallDir', type=str, help="basecallDir dir name", default=None)
    # parser.add_argument('--methcallDir', type=str, help="methcallDir dir name", default=None)
    # parser.add_argument('--processors', type=int, help="Number of processors", default=8)
    # parser.add_argument('--mpi', action='store_true')

    return parser.parse_args()


def import_call(fn, callname):
    """
    Import fn as callname and return
        call0   -   original results
        call1   -   covrage cutoff results
    :param fn:
    :param callname:
    :return:
    """
    if callname == 'DeepSignal':
        logger.debug(f"Start load DeepSignal")
        calls0 = importPredictions_DeepSignal(fn, baseFormat=baseFormat)

    elif callname == 'Tombo':
        logger.debug(f"Start load Tombo")
        calls0 = importPredictions_Tombo(fn, baseFormat=baseFormat)
    elif callname == 'Nanopolish':
        logger.debug(f"Start load Nanopolish")
        calls0 = importPredictions_Nanopolish(fn, baseFormat=baseFormat, logLikehoodCutt=2.0)
    elif callname == 'DeepMod':
        logger.debug(f"Start load DeepMod")
        calls0 = importPredictions_DeepMod_Read_Level(fn, baseFormat=baseFormat)
    else:
        raise Exception(f'Not support {callname} for file {fn} now')

    calls1 = coverageFiltering(calls0, minCov=minToolCovCutt, toolname=callname)
    logger.debug(f'Import {name} finished!')
    return [calls0, calls1]


def import_bgtruth(fn, encode):
    logger.debug(f"Start load bgTruth")

    if encode == "encode":
        bgTruth = importGroundTruth_BedMethyl_from_Encode(fn, covCutt=bgtruthCutt, baseCount=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/K562/ENCSR765JPC_WGBS_hg38/ENCFF721JMB.bed", chrFilter="chr20")
    elif encode == "oxBS_sudo":
        bgTruth = importGroundTruth_oxBS(fn, covCutt=bgtruthCutt, baseCount=baseFormat)
    elif encode == "bed":  # like K562 bg-truth
        bgTruth = importGroundTruth_coverage_output_from_Bismark(fn, covCutt=bgtruthCutt, baseFormat=baseFormat, includeCov=True)
    elif encode == "bismark_bedgraph":
        bgTruth = importGroundTruth_coverage_output_from_Bismark_BedGraph(fn, baseCount=baseFormat)
    elif encode == "bismark":  # for genome-wide Bismark results, such as HL60, etc.
        bgTruth = importGroundTruth_genome_wide_output_from_Bismark(fn, covCutt=bgtruthCutt, baseCount=baseFormat)
    else:
        raise Exception("Methylation_correlation_plotting.py ERROR: Unknown bacground truth parser configuration. Aborting. FYI: currently supported are: encode, oxBS_sudo, bismark")

    logger.debug(f'Import BG-Truth finished!')
    return bgTruth


def summary_table():
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
    df.to_csv(outfn)
    logger.info(f'save to {outfn}')


if __name__ == '__main__':
    set_log_debug_level()

    ### python Methylation_correlation_plotting.py  --calls a:3234 b:33433 --bgtruth c:3242424
    args = parse_arguments()
    logger.debug(args)

    RunPrefix = args.runid  # "K562_WGBS_rep_ENCFF721JMB"
    out_dir = os.path.join(pic_base_dir, RunPrefix)
    os.makedirs(out_dir, exist_ok=True)

    # tool coverage cutoff 4
    minToolCovCutt = 4

    # bgtruth coverage cutoff 10
    bgtruthCutt = 10

    # load into program format 0-base or 1-base
    baseFormat = 0

    callfn_dict = defaultdict()  # callname -> filename
    callresult_dict = defaultdict()
    callname_list = []

    with Pool(processes=args.processors) as pool:
        for callstr in args.calls:
            callname, callfn = callstr.split(':')
            callfn_dict[callname] = callfn
            callname_list.append(callname)
            # callresult_dict[callname] = import_call(callfn, callname)
            callresult_dict[callname] = pool.apply_async(import_call, (callfn, callname,))

        logger.info(args.bgtruth.split(':'))

        encode = args.bgtruth.split(':')[0]
        fn = args.bgtruth.split(':')[1]

        logger.info(f'fn={fn}, encode={encode}')
        # bgTruth = import_bgtruth(fn, encode)
        bgTruth = pool.apply_async(import_bgtruth, (fn, encode,))
        pool.close()
        pool.join()
    logger.info("Join processes of import")

    for name in callname_list:
        callresult_dict[callname] = callresult_dict[callname].get()
    bgTruth = bgTruth.get()

    # logger.info(callfn_dict)
    # logger.info(callname_list)

    logger.debug(f"Start set intersection with all tools joined together with bgtruth")

    coveredCpGs = set(list(bgTruth.keys()))
    for name in callname_list:
        coveredCpGs = coveredCpGs.intersection(set(list(callresult_dict[name][1].keys())))
        logger.info(f'Join {name} get {len(coveredCpGs)} CpGs')
    logger.info(f"Reporting {len(coveredCpGs)} CpGs are covered by all tools and bgtruth")

    logger.debug('Output data of coverage and meth-freq on joined CpG sites for correlation analysis')
    outfn = os.path.join(out_dir, f"Meth_corr_plot_data-{RunPrefix}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.tsv")
    logger.info(f"Start output results to {outfn}")

    outfile = open(outfn, 'w')
    outfile.write("chr\tstart\tend\tBSseq_freq\tBSseq_cov\tstrand")

    for name in callname_list:
        outfile.write(f"\t{name}_freq\t{name}_cov")
    outfile.write("\n")

    for cpg in coveredCpGs:
        outfile.write(f"{cpg[0]}\t{cpg[1]}\t{cpg[1] + 1}\t{bgTruth[cpg][0]}\t{bgTruth[cpg][1]}\t{cpg[2]}")

        for name in callname_list:
            outfile.write(f"\t{callresult_dict[name][1][cpg][0]}\t{callresult_dict[name][1][cpg][1]}")
        outfile.write("\n")
    outfile.close()
    logger.info(f"save to {outfn}")

    summary_table()

    logger.info("Methylation correlation plotting data generation program finished.")

    sys.exit(0)

    is_scatter_plot = False
    output_bed = False
    study_high_low_cov = False

    logger.info(f'bgtruth cov={bgtruthCutt}, tool cov={minToolCovCutt}, baseFormat={baseFormat}')

    out_dir = os.path.join(pic_base_dir, RunPrefix)
    os.makedirs(out_dir, exist_ok=True)

    logger.debug(list(enumerate(argv)))

    logger.debug(f"Start load DeepSignal")
    if argv[1] == 'NO':
        DeepSignal_calls = None
    else:
        DeepSignal_calls0 = importPredictions_DeepSignal(argv[1], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepSignal_runs/K562/K562_DeepSignal.MethCalls.Joined.tsv")
        DeepSignal_calls = coverageFiltering(DeepSignal_calls0, minCov=minToolCovCutt)

    logger.debug(f"Start load Tombo")
    if argv[2] == 'NO':
        Tombo_calls = None
    else:
        Tombo_calls0 = importPredictions_Tombo(argv[2], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/K562/K562_Tombo.batch_all.batch_0.perReadsStats.bed")
        Tombo_calls = coverageFiltering(Tombo_calls0, minCov=minToolCovCutt)

    logger.debug(f"Start load Nanopolish")
    if argv[3] == 'NO':
        Nanopolish_calls = None
    else:
        # Nanopolish_calls = importPredictions_Nanopolish_v2(argv[3])  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.tsv", IncludeNonSingletons = True)
        Nanopolish_calls0 = importPredictions_Nanopolish(argv[3], baseFormat=baseFormat, logLikehoodCutt=2.0)
        Nanopolish_calls = coverageFiltering(Nanopolish_calls0, minCov=minToolCovCutt)

    logger.debug(f"Start load DeepMod")
    if argv[4] == 'NO':
        DeepMod_calls = None
    else:
        DeepMod_calls0 = importPredictions_DeepMod_Read_Level(argv[4], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/K562/K562.C.combined.bed")
        DeepMod_calls = coverageFiltering(DeepMod_calls0, minCov=minToolCovCutt)

    logger.debug(f"Start load DeepMod_clustered")
    if argv[5] == 'NO':
        DeepMod_calls_clustered = None
    else:
        DeepMod_calls_clustered0 = importPredictions_DeepMod_clustered(argv[5], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/K562/K562.C_clusterCpG.combined.bed")
        DeepMod_calls_clustered = coverageFiltering(DeepMod_calls_clustered0, minCov=minToolCovCutt, byLength=False)

    name_calls = ['DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'DeepMod_cluster']
    all_calls = [DeepSignal_calls, Tombo_calls, Nanopolish_calls, DeepMod_calls, DeepMod_calls_clustered]
    callfn_dict = {name_calls[k]: all_calls[k] for k in range(len(all_calls))}

    all_calls_before_cutoff = [DeepSignal_calls0, Tombo_calls0, Nanopolish_calls0, DeepMod_calls0, DeepMod_calls_clustered0]

    if study_high_low_cov:
        study_high_cov_tool1_low_cov_tool2()

    if output_bed:
        output_bed_of_tools_bgtruth()

    if is_scatter_plot:
        scatter_plot_cov()

    study_intersection()
