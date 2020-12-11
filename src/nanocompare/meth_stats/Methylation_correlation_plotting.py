#!/home/liuya/anaconda3/envs/nmf/bin/python

"""
Generate methy correlation plotting input data as a tsv

Sample usage:
    <pyfilename> $DeepSignal_calls $Tombo_calls $Nanopolish_calls \
			$DeepMod_calls $DeepMod_cluster_calls $bgTruth $parser $RunPrefix

All usedful functions are located in nanocompare.meth_stats.meth_stats_common
"""
import sys

from sys import argv

nanocompare_prj = "/projects/li-lab/yang/workspace/nano-compare/src"
sys.path.append(nanocompare_prj)

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

    logger.info(f'bgtruth cov={bgtruthCutt}, tool cov={minToolCovCutt}')

    RunPrefix = argv[8]  # "K562_WGBS_rep_ENCFF721JMB"
    out_dir = os.path.join(pic_base_dir, RunPrefix)
    os.makedirs(out_dir, exist_ok=True)

    logger.debug(list(enumerate(argv)))

    logger.debug(f"Start load DeepSignal")
    if argv[1] == 'NO':
        DeepSignal_calls = None
    else:
        DeepSignal_calls = importPredictions_DeepSignal(argv[1], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepSignal_runs/K562/K562_DeepSignal.MethCalls.Joined.tsv")
        DeepSignal_calls = coverageFiltering(DeepSignal_calls, minCov=minToolCovCutt)

    logger.debug(f"Start load Tombo")
    if argv[2] == 'NO':
        Tombo_calls = None
    else:
        Tombo_calls = importPredictions_Tombo(argv[2], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/K562/K562_Tombo.batch_all.batch_0.perReadsStats.bed")
        Tombo_calls = coverageFiltering(Tombo_calls, minCov=minToolCovCutt)

    logger.debug(f"Start load Nanopolish")
    if argv[3] == 'NO':
        Nanopolish_calls = None
    else:
        # Nanopolish_calls = importPredictions_Nanopolish_v2(argv[3])  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.tsv", IncludeNonSingletons = True)
        Nanopolish_calls = importPredictions_Nanopolish(argv[3], baseFormat=baseFormat)
        Nanopolish_calls = coverageFiltering(Nanopolish_calls, minCov=minToolCovCutt)

    logger.debug(f"Start load DeepMod")
    if argv[4] == 'NO':
        DeepMod_calls = None
    else:
        DeepMod_calls = importPredictions_DeepMod(argv[4], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/K562/K562.C.combined.bed")
        DeepMod_calls = coverageFiltering(DeepMod_calls, minCov=minToolCovCutt)

    logger.debug(f"Start load DeepMod_clusteredResultParsing")
    if argv[5] == 'NO':
        DeepMod_calls_clustered = None
    else:
        DeepMod_calls_clustered = importPredictions_DeepMod_clustered(argv[5], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/K562/K562.C_clusterCpG.combined.bed")
        DeepMod_calls_clustered = coverageFiltering(DeepMod_calls_clustered, minCov=minToolCovCutt, byLength=False)

    logger.debug(f"Start load bgTruth")

    if argv[7] == "encode":
        bgTruth = importGroundTruth_BedMethyl_from_Encode(argv[6], covCutt=bgtruthCutt, baseCount=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/K562/ENCSR765JPC_WGBS_hg38/ENCFF721JMB.bed", chrFilter="chr20")
    elif argv[7] == "oxBS_sudo":
        bgTruth = importGroundTruth_oxBS(argv[6], covCutt=bgtruthCutt, baseCount=baseFormat)
    elif argv[7] == "bismark":
        bgTruth = importGroundTruth_coverage_output_from_Bismark(argv[6], covCutt=bgtruthCutt, baseFormat=baseFormat, includeCov=True)
    elif argv[7] == "bismark_bedgraph":
        bgTruth = importGroundTruth_coverage_output_from_Bismark_BedGraph(argv[6], baseCount=baseFormat)
    else:
        raise Exception("Methylation_correlation_plotting.py ERROR: Unknown bacground truth parser configuration. Aborting. FYI: currently supported are: encode, oxBS_sudo, bismark")
        sys.exit(-1)

    name_calls = ['DeepSignal', 'Tombo', 'Nanopolish', 'DeepMod', 'DeepMod_cluster']
    all_calls = [DeepSignal_calls, Tombo_calls, Nanopolish_calls, DeepMod_calls, DeepMod_calls_clustered]

    logger.debug(f"Output bed files of each tool, and bgtruth")

    # Note bgTruth in format of {'chr\t123\t123\n':[0.56, 15], etc.}
    outfn = os.path.join(out_dir, f'{RunPrefix}-meth-cov-bgtruth-baseCount{baseFormat}.bed')
    save_call_or_bgtruth_to_bed(bgTruth, outfn)

    # Note each tool call filtered results is in format of {'chr\t123\t123\n':[0.56, 15], etc.}
    for call1, name1 in zip(all_calls, name_calls):
        outfn = os.path.join(out_dir, f'{RunPrefix}-meth-cov-{name1}-baseCount{baseFormat}.bed')
        save_call_or_bgtruth_to_bed(call1, outfn)

    logger.debug(f"Study set intersection of each tool with bgtruth")
    dataset = []
    bgtruthCpGs = set(list(bgTruth.keys()))
    for call1, name1 in zip(all_calls, name_calls):
        if call1 is None:
            dataset.append({'bgtruth': len(bgtruthCpGs), 'tool': np.nan, 'joined': np.nan})
            continue
        overlapCpGs = bgtruthCpGs.intersection(set(list(call1.keys())))
        dataset.append({'bgtruth': len(bgtruthCpGs), 'tool': len(set(list(call1.keys()))), 'joined': len(overlapCpGs)})
        logger.info(f'BG-Truth join with {name1} get {len(overlapCpGs)} CpGs')
        outfn = os.path.join(out_dir, f'{RunPrefix}-joined-cpgs-bgtruth-{name1}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-baseCount{baseFormat}.bed')
        save_keys_to_bed(overlapCpGs, outfn)

    df = pd.DataFrame(dataset, index=name_calls)
    outfn = os.path.join(out_dir, f'{RunPrefix}-summary-bgtruth-tools-bsCov{bgtruthCutt}-minCov{minToolCovCutt}.csv')
    df.to_csv(outfn)

    logger.debug(f"Start set intersection with all joined together (4+1 tools with bgtruth)")
    coveredCpGs = set(list(bgTruth.keys()))
    for call1, name1 in zip(all_calls, name_calls):
        if call1 is None:
            continue
        coveredCpGs = coveredCpGs.intersection(set(list(call1.keys())))
        logger.info(f'Join {name1} get {len(coveredCpGs)} CpGs')

    logger.info(f"{len(coveredCpGs)} CpGs are covered by all tools and bgtruth")

    outfn = os.path.join(out_dir, f"Meth_corr_plot_data-{RunPrefix}-bsCov{bgtruthCutt}-minCov{minToolCovCutt}-time-{current_time_str()}.tsv")
    logger.info(f"Start output results to {outfn}")

    outfile = open(outfn, 'w')
    outfile.write("chr\tstart\tend\tstrand\tDeepSignal_freq\tDeepSignal_cov\tTombo_freq\tTombo_cov\tNanopolish_freq\tNanopolish_cov\tDeepMod_freq\tDeepMod_cov\tDeepMod_clust_freq\tDeepMod_clust_cov\tBSseq\n")

    for cpg in coveredCpGs:
        coords = cpg.strip().split("\t")
        outfile.write("{}\t{}\t{}\t{}\t".format(coords[0], coords[1], coords[2], coords[3]))

        for call1 in all_calls:
            if call1 is None:
                outfile.write("\t\t")
                continue
            outfile.write(f"{call1[cpg][0]}\t{call1[cpg][1]}\t")

        outfile.write(f"{bgTruth[cpg][0]}\n")

    outfile.close()
    logger.info(f"save to {outfn}")
