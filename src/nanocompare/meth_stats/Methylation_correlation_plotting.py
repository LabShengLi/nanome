#!/home/liuya/anaconda3/envs/nmf/bin/python

"""
Generate methy correlation plotting input data as a tsv

Sample usage:
    <pyfilename> $DeepSignal_calls $Tombo_calls $Nanopolish_calls \
			$DeepMod_calls $DeepMod_cluster_calls $bgTruth $parser $RunPrefix

All usedful functions are located in nanocompare.meth_stats.meth_stats_common
"""
import sys

nanocompare_prj = "/projects/li-lab/yang/workspace/nano-compare/src"
sys.path.append(nanocompare_prj)

from sys import argv

from nanocompare.meth_stats.meth_stats_common import *

from global_config import *

if __name__ == '__main__':
    set_log_debug_level()

    # tool coverage cutoff 4
    minToolCovCutt = 4

    # bgtruth coverage cutoff 10
    bgtruthCutt = 5

    # load into program format 0-base or 1-base
    baseFormat = 0

    logger.info(f'bgtruth cov={bgtruthCutt}, tool cov={minToolCovCutt}, baseFormat={baseFormat}')

    RunPrefix = argv[8]  # "K562_WGBS_rep_ENCFF721JMB"
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
        DeepMod_calls0 = importPredictions_DeepMod(argv[4], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/K562/K562.C.combined.bed")
        DeepMod_calls = coverageFiltering(DeepMod_calls0, minCov=minToolCovCutt)

    logger.debug(f"Start load DeepMod_clustered")
    if argv[5] == 'NO':
        DeepMod_calls_clustered = None
    else:
        DeepMod_calls_clustered0 = importPredictions_DeepMod_clustered(argv[5], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/K562/K562.C_clusterCpG.combined.bed")
        DeepMod_calls_clustered = coverageFiltering(DeepMod_calls_clustered0, minCov=minToolCovCutt, byLength=False)

    logger.debug(f"Start load bgTruth")

    if argv[7] == "encode":
        bgTruth = importGroundTruth_BedMethyl_from_Encode(argv[6], covCutt=bgtruthCutt, baseCount=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/K562/ENCSR765JPC_WGBS_hg38/ENCFF721JMB.bed", chrFilter="chr20")
    elif argv[7] == "oxBS_sudo":
        bgTruth = importGroundTruth_oxBS(argv[6], covCutt=bgtruthCutt, baseCount=baseFormat)
    elif argv[7] == "bed":  # like K562 bg-truth
        bgTruth = importGroundTruth_coverage_output_from_Bismark(argv[6], covCutt=bgtruthCutt, baseFormat=baseFormat, includeCov=True)
    elif argv[7] == "bismark_bedgraph":
        bgTruth = importGroundTruth_coverage_output_from_Bismark_BedGraph(argv[6], baseCount=baseFormat)
    elif argv[7] == "bismark":  # for genome-wide Bismark results, such as HL60, etc.
        bgTruth = importGroundTruth_genome_wide_output_from_Bismark(argv[6], covCutt=bgtruthCutt, baseCount=baseFormat)
    else:
        raise Exception("Methylation_correlation_plotting.py ERROR: Unknown bacground truth parser configuration. Aborting. FYI: currently supported are: encode, oxBS_sudo, bismark")
        sys.exit(-1)

    name_calls = ['DeepSignal', 'Nanopolish', 'DeepMod', 'Tombo', 'DeepMod_cluster']
    all_calls = [DeepSignal_calls, Nanopolish_calls, DeepMod_calls, Tombo_calls, DeepMod_calls_clustered]

    all_calls_before_cutoff = [DeepSignal_calls0, Nanopolish_calls0, DeepMod_calls0, Tombo_calls0, DeepMod_calls_clustered0]

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

    logger.debug(f"Study scatter plot of coverage for Tombo with other tool with no cutoff")

    scatter_analysis_cov(Tombo_calls0, Nanopolish_calls0, outdir=out_dir, RunPrefix=RunPrefix, tool1_name='Tombo', tool2_name='Nanopolish')
    scatter_analysis_cov(Tombo_calls0, DeepSignal_calls0, outdir=out_dir, RunPrefix=RunPrefix, tool1_name='Tombo', tool2_name='Deepsignal')

    scatter_analysis_cov(DeepMod_calls0, Nanopolish_calls0, outdir=out_dir, RunPrefix=RunPrefix, tool1_name='DeepMod', tool2_name='Nanopolish')
    scatter_analysis_cov(DeepMod_calls0, DeepSignal_calls0, outdir=out_dir, RunPrefix=RunPrefix, tool1_name='DeepMod', tool2_name='Deepsignal')

    logger.debug(f"Study set intersection of each tool with bgtruth")
    dataset = []
    bgtruthCpGs = set(list(bgTruth.keys()))
    for call1, name1, call1_before_cutoff in zip(all_calls, name_calls, all_calls_before_cutoff):
        if call1 is None:
            dataset.append({'bgtruth': len(bgtruthCpGs), 'tool': np.nan, f'tool-covcut{minToolCovCutt}': np.nan, 'joined': np.nan})
            continue
        overlapCpGs = bgtruthCpGs.intersection(set(list(call1.keys())))
        ret = {'CpG sites in BG-Truth': len(bgtruthCpGs), 'Total CpG sites by Nanopore tool': len(set(list(call1_before_cutoff))), f'Total CpG sites by cov>={minToolCovCutt}': len(set(list(call1.keys()))), 'Joined CpG sites with BG-Truth': len(overlapCpGs)}

        cnt_calls = 0

        if name1 != "DeepMod_cluster":
            for cpg in call1_before_cutoff:
                cnt_calls += len(call1_before_cutoff[cpg])
        else:
            for cpg in call1_before_cutoff:
                cnt_calls += call1_before_cutoff[cpg][1]
        ret.update({'Total reads by Nanopore tool': cnt_calls})
        dataset.append(ret)

        logger.info(f'BG-Truth join with {name1} get {len(overlapCpGs):,} CpGs')
        outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-bgtruth-{name1}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
        save_keys_to_bed(overlapCpGs, outfn)

    df = pd.DataFrame(dataset, index=name_calls)
    outfn = os.path.join(out_dir, f'{RunPrefix}-summary-bgtruth-tools-bsCov{bgtruthCutt}-minCov{minToolCovCutt}.csv')
    df.to_csv(outfn)

    logger.debug(f"Study set intersection of deepsignal, nanopolish with bgtruth")
    dataset = []
    overlapCpGs = set(list(bgTruth.keys()))
    for call1, name1 in zip(all_calls, name_calls):
        if call1 is None:
            continue
        if name1 not in ['DeepSignal', 'Nanopolish']:
            continue
        overlapCpGs = overlapCpGs.intersection(set(list(call1.keys())))
    outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-deepsignal-nanopolish-bgtruth-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
    save_keys_to_bed(overlapCpGs, outfn)
    logger.info(f'Reporting deepsignal, nanopolish with bgtruth, joined results = {len(overlapCpGs)}')

    logger.debug(f"Study set intersection of deepsignal, nanopolish")
    dataset = []
    overlapCpGs = None
    for call1, name1 in zip(all_calls, name_calls):
        if call1 is None:
            continue
        if name1 not in ['DeepSignal', 'Nanopolish']:
            continue
        if overlapCpGs is None:
            overlapCpGs = set(list(call1.keys()))
        else:
            overlapCpGs = overlapCpGs.intersection(set(list(call1.keys())))
    outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-deepsignal-nanopolish-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
    save_keys_to_bed(overlapCpGs, outfn)
    logger.info(f'Reporting deepsignal, nanopolish joined results = {len(overlapCpGs)}')

    logger.debug(f"Next, study set intersection of deepsignal, nanopolish but not in Tombo")
    newOverlapCpGs = overlapCpGs - Tombo_calls.keys()

    outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-deepsignal-nanopolish-notombo-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
    save_keys_to_bed(newOverlapCpGs, outfn)
    logger.info(f'Reporting deepsignal, nanopolish but not in Tombo results = {len(newOverlapCpGs)}')

    logger.debug(f"Start set intersection with all joined together (4+1 tools with bgtruth)")
    coveredCpGs = set(list(bgTruth.keys()))

    name_calls_new_order = ['DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'DeepMod_cluster']

    # for call1, name1 in zip(all_calls, name_calls):

    call_dict = {name_calls[k]: all_calls[k] for k in range(len(all_calls))}
    for name1 in name_calls_new_order:
        if call_dict[name1] is None:
            continue
        coveredCpGs = coveredCpGs.intersection(set(list(call_dict[name1].keys())))
        logger.info(f'Join {name1} get {len(coveredCpGs)} CpGs')
    logger.info(f"{len(coveredCpGs)} CpGs are covered by all tools and bgtruth")

    outfn = os.path.join(out_dir, f"Meth_corr_plot_data-{RunPrefix}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.tsv")
    logger.info(f"Start output results to {outfn}")

    outfile = open(outfn, 'w')
    outfile.write("chr\tstart\tend\tBSseq_freq\tBSseq_cov\tstrand\tDeepSignal_freq\tDeepSignal_cov\tTombo_freq\tTombo_cov\tNanopolish_freq\tNanopolish_cov\tDeepMod_freq\tDeepMod_cov\tDeepMod_clust_freq\tDeepMod_clust_cov\n")

    for cpg in coveredCpGs:
        # coords = cpg.strip().split("\t")
        outfile.write(f"{cpg[0]}\t{cpg[1]}\t{cpg[1] + 1}\t{bgTruth[cpg][0]}\t{bgTruth[cpg][1]}\t{cpg[2]}")

        for call1 in all_calls[:]:
            if call1 is None:
                outfile.write("\t.\t.")
                continue
            outfile.write(f"\t{call1[cpg][0]}\t{call1[cpg][1]}")
        outfile.write("\n")
    outfile.close()
    logger.info(f"save to {outfn}")
    logger.info("Methylation correlation plotting data generation program finished.")
