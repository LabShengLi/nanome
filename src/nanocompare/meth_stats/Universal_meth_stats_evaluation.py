"""

Evaluation based on methylation calls of four tools, compute the performance results(F1, accuracy, etc.)

/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepSignal_runs/NA19240/NA19240_DeepSignal.MethCalls.Joined.tsv /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/NA19240/NA19240_Tombo.batch_all.batch_0.perReadsStats.bed /projects/li-lab/NanoporeData/WR_ONT_analyses/NA19240_nanopolish/NA19240.methylation_calls.tsv /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/NA19240/NA19240.C.combined.bed /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/joined_reps/RRBS/extractBismark/NA19240_joined_RRBS.Read_R1.Rep_1_trimmed_bismark_bt2.bismark.cov.gz NA19240_RRBS_joined/NA19240_RRBS_joined bismark 10 Tombo


"""

"""
   This script will generate all performance results, bed files of singleton, non-singleton, based on results by DL call tool related to BGTruth results
"""

# this script was based on this jupyter notebook: http://helix067:9912/notebooks/li-lab/NanoporeData/WR_ONT_analyses/ai/GitHub/nanoCompare/Scripts/UniversalMethStatsEvaluation.ipynb
# copied Wed Jul 10 14:27:33 EDT 2019
# example command run: python /projects/li-lab/NanoporeData/WR_ONT_analyses/ai/GitHub/nanoCompare/Scripts/UniversalMethStatsEvaluation.standalone_01.py $DeepSignal_calls $Tombo_calls $Nanopolish_calls $bgTruth $RunPrefix $parser

# where variable names came from:
# DeepSignal_calls = importPredictions_DeepSignal(argv[1])#"/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepSignal_runs/K562/K562_DeepSignal.MethCalls.Joined.chr20.tsv")
# Tombo_calls = importPredictions_Tombo(argv[2])#"/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/K562/K562_Tombo.perReadsStats.chr20.bed")
# Nanopolish_calls = importPredictions_Nanopolish_2(argv[3])#"/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.chr_chr20.tsv", IncludeNonSingletons = True)
# bgTruth = importGroundTruth_BedMethyl_from_Encode(argv[4])#"/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/K562/ENCSR765JPC_WGBS_hg38/ENCFF721JMB.bed", chrFilter="chr20")
# RunPrefix = argv[5]#"K562_WGBS_rep_ENCFF721JMB/K562_WGBS_rep_ENCFF721JMB"
# parser = "encode"

####################

import sys

nanocompare_prj = "/projects/li-lab/yang/workspace/nano-compare/src"
sys.path.append(nanocompare_prj)

from nanocompare.nanocompare_global_settings import singletonsFile, narrowCoord, nonsingletonsFile
from nanocompare.meth_stats.meth_stats_common import *

from global_config import *

from sys import argv

if __name__ == '__main__':

    set_log_debug_level()

    baseFormat = 1

    minCov = int(argv[8])
    dsname = argv[9]
    RunPrefix = argv[6]  # "K562_WGBS_rep_ENCFF721JMB/K562_WGBS_rep_ENCFF721JMB"

    out_dir = os.path.join(pic_base_dir, RunPrefix)
    os.makedirs(out_dir, exist_ok=True)

    tool = argv[9]

    logger.info(f"Start DeepSignal")
    DeepSignal_calls = importPredictions_DeepSignal(argv[1], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepSignal_runs/K562/K562_DeepSignal.MethCalls.Joined.chr20.tsv")

    logger.info(f"Start Tombo")
    Tombo_calls = importPredictions_Tombo(argv[2], baseFormat=baseFormat, cutoff=2.0)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/K562/K562_Tombo.perReadsStats.chr20.bed")

    logger.info(f"Start Nanopolish")
    Nanopolish_calls = importPredictions_Nanopolish(argv[3], baseFormat=baseFormat, logLikehoodCutt=2.0)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.chr_chr20.tsv", IncludeNonSingletons = True)

    logger.info(f"Start DeepMod")
    DeepMod_calls = importPredictions_DeepMod(argv[4], baseFormat=baseFormat)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.chr_chr20.tsv", IncludeNonSingletons = True)

    logger.info(f"Start bgTruth")
    if argv[7] == "encode":
        bgTruth = importGroundTruth_BedMethyl_from_Encode(argv[5], covCutt=minCov)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/K562/ENCSR765JPC_WGBS_hg38/ENCFF721JMB.bed", chrFilter="chr20")
    elif argv[7] == "oxBS_sudo":
        bgTruth = importGroundTruth_oxBS(argv[5], covCutt=minCov)
    elif argv[7] == "bismark":
        bgTruth = importGroundTruth_coverage_output_from_Bismark(argv[5], baseFormat=baseFormat, covCutt=minCov)
    else:
        logger.error("UniversalMethStatsEvaluation.standalone_01.py ERROR: Unknown bacground truth parser configuration. Aborting. FYI: currently supported are: encode, oxBS_sudo, bismark")
        exit()

    ####################################### Preprocessing neccessary files:

    # singletonsFile = "hg38_singletons.bed"
    # nonsingletonsFile = "hg38_nonsingletons.bed"
    #
    # narrowCoord = [False, singletonsFile, nonsingletonsFile, "ONT.hg38.cpgIslandExt.bed", "ONT.hg38.cpgShoresExt.bed", "ONT.hg38.cpgShelvesExt.bed", "ONT.hg38.exonFeature.bed", "ONT.hg38.geneFeature.bed", "ONT.hg38.intergenic.bed", "ONT.hg38.intronFeature.bed", "ONT.hg38.promoterFeature.flank_100.bed", "ONT.hg38.promoterFeature.flank_1000.bed",
    #         "ONT.hg38.promoterFeature.flank_200.bed", "ONT.hg38.promoterFeature.flank_2000.bed", "ONT.hg38.promoterFeature.flank_500.bed", "ONT.hg38.promoterFeature.flank_750.bed"]

    relateCoord = narrowCoord
    ## add missing region files:
    singletonsFilePrefix = singletonsFile.replace(".bed", '')
    # relateCoord.append("{}/{}.{}.mixed.bed".format(out_dir, RunPrefix, singletonsFilePrefix))
    relateCoord.append("{}/{}.{}.absolute.bed".format(out_dir, RunPrefix, singletonsFilePrefix))

    nonsingletonsFilePrefix = nonsingletonsFile.replace(".bed", '')
    # relateCoord.append("{}/{}.{}.other.bed".format(out_dir, RunPrefix, nonsingletonsFilePrefix))
    # relateCoord.append("{}/{}.{}.fullyMixed.bed".format(out_dir, RunPrefix, nonsingletonsFilePrefix))
    relateCoord.append("{}/{}.{}.discordant.bed".format(out_dir, RunPrefix, nonsingletonsFilePrefix))
    relateCoord.append("{}/{}.{}.concordant.bed".format(out_dir, RunPrefix, nonsingletonsFilePrefix))

    singletonsPostprocessing(bgTruth, singletonsFile, RunPrefix, outdir=out_dir)
    nonSingletonsPostprocessing(bgTruth, nonsingletonsFile, RunPrefix, outdir=out_dir)

    ####################################### Computing stats per program:

    logger.info("\n\n############\n\n")

    DeepSignal_Tombo = combine2programsCalls(DeepSignal_calls, Tombo_calls)
    DeepSignal_Tombo_Nanopolish = combine2programsCalls(DeepSignal_Tombo, Nanopolish_calls)
    DeepSignal_Tombo_Nanopolish_DeepMod = combine2programsCalls(DeepSignal_Tombo_Nanopolish, DeepMod_calls)
    DeepSignal_Tombo_Nanopolish_DeepMod_Bacground = combine2programsCalls(DeepSignal_Tombo_Nanopolish_DeepMod, bgTruth, outfileName="{}/{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed".format(out_dir, RunPrefix))

    logger.info(f"Data points for joined 4 tools with bg-truth stats: {len(DeepSignal_Tombo_Nanopolish_DeepMod_Bacground)}")

    DeepSignal_Tombo_corr = combine2programsCalls_4Corr(DeepSignal_calls, Tombo_calls)
    DeepSignal_Tombo_Nanopolish_corr = combine2programsCalls_4Corr(DeepSignal_Tombo_corr, Nanopolish_calls)
    DeepSignal_Tombo_Nanopolish_DeepMod_corr = combine2programsCalls_4Corr(DeepSignal_Tombo_Nanopolish_corr, DeepMod_calls)
    # DeepSignal_Tombo_Nanopolish_Bacground = combine2programsCalls_4Corr(DeepSignal_Tombo_Nanopolish_corr, oxBS, outfileName = "DeepSignal_Tombo_Nanopolish_Bacground.4Corr.bed")
    DeepSignal_Tombo_Nanopolish_DeepMod_Bacground_corr = combine2programsCalls(DeepSignal_Tombo_Nanopolish_DeepMod_corr, bgTruth, outfileName="{}/{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.4Corr.bed".format(out_dir, RunPrefix))

    logger.info(f"Data points for correlation: {len(DeepSignal_Tombo_Nanopolish_DeepMod_Bacground_corr)}")

    logger.info("\n\n############\n\n")

    # RunPrefix = f'{RunPrefix}/{RunPrefix}'

    secondFilterBed = "{}/{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed".format(out_dir, RunPrefix)  # If using joined together CpGs sites, or False

    # secondFilterBed = False

    secondFilterBed_4Corr = "{}/{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.4Corr.bed".format(out_dir, RunPrefix)

    toolCalls = {'DeepSignal': DeepSignal_calls, 'Tombo': Tombo_calls, 'Nanopolish': Nanopolish_calls, 'DeepMod': DeepMod_calls}

    logger.debug(relateCoord)  # all coordinate generated

    perf_dir = os.path.join(out_dir, 'performance-results-joined')
    os.makedirs(perf_dir, exist_ok=True)

    perf_dir_nojoined = os.path.join(out_dir, 'performance-results-nojoined')
    os.makedirs(perf_dir_nojoined, exist_ok=True)

    for tool in toolCalls:
        tmpPrefix = f'{RunPrefix}.{tool}'
        logger.info(tmpPrefix)

        # Note: relateCoord - all singleton (absolute and mixed) and non-singleton generated bed. ranges
        #       secondFilterBed - joined sites of four tools and bg-truth. points

        # step 1: with joined results of all tools
        df = report_per_read_performance(toolCalls[tool], bgTruth, tmpPrefix, narrowedCoordinatesList=relateCoord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)
        df = df[["prefix", "coord", "accuracy", "roc_auc", "F1_5C", "F1_5mC", "precision_5C", "recall_5C", "precision_5mC", "recall_5mC", "corrMix", "Corr_mixedSupport", "corrAll", "Corr_allSupport", "Csites_called", "mCsites_called", "referenceCpGs", "Csites", "mCsites"]]

        df['Tool'] = tool
        df['Dataset'] = dsname

        outfn = os.path.join(perf_dir, f"{RunPrefix}.{tool}.report.tsv")
        df.to_csv(outfn, sep='\t')
        logger.info(f"save to {outfn}")
        # logger.debug(df.head())

        # step 2: no joined results eval
        df = report_per_read_performance(toolCalls[tool], bgTruth, tmpPrefix, narrowedCoordinatesList=relateCoord, secondFilterBed=False, secondFilterBed_4Corr=secondFilterBed_4Corr)
        df = df[["prefix", "coord", "accuracy", "roc_auc", "F1_5C", "F1_5mC", "precision_5C", "recall_5C", "precision_5mC", "recall_5mC", "corrMix", "Corr_mixedSupport", "corrAll", "Corr_allSupport", "Csites_called", "mCsites_called", "referenceCpGs", "Csites", "mCsites"]]

        df['Tool'] = tool
        df['Dataset'] = dsname

        outfn = os.path.join(perf_dir_nojoined, f"{RunPrefix}.{tool}.report.tsv")
        df.to_csv(outfn, sep='\t')
        logger.info(f"save to {outfn}")
        # logger.debug(df.head())

    # sys.exit(0)
    #
    # # tool="DeepSignal"
    # if tool == "DeepSignal":
    #     prefix = "{}.{}_calls.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground".format(RunPrefix, tool)
    #     secondFilterBed = "{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed".format(RunPrefix)  # filter coordinates for stats like AUC, F1 etc.
    #     secondFilterBed_4Corr = "{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.4Corr.bed".format(RunPrefix)  # filter coordinates for correlation
    #     df = combine_ONT_and_BS(DeepSignal_calls, bgTruth, prefix, narrowedCoordinates=narrowCoord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)
    #     df = df[["prefix", "coord", "accuracy", "roc_auc", "precision_5C", "recall_5C", "F1_5C", "Csites", "precision_5mC", "recall_5mC", "F1_5mC", "mCsites", "referenceCpGs", "corrMix", "Corr_mixedSupport", "corrAll", "Corr_allSupport"]]
    #     df.to_csv("{}.report.tsv".format(prefix), sep='\t')
    #     logger.info(df.head())
    #
    # # tool="Tombo"
    # if tool == "Tombo":
    #     prefix = "{}.{}_calls.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground".format(RunPrefix, tool)
    #     secondFilterBed = "{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed".format(RunPrefix)  # filter coordinates for stats like AUC, F1 etc.
    #     secondFilterBed_4Corr = "{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.4Corr.bed".format(RunPrefix)  # filter coordinates for correlation
    #     df = combine_ONT_and_BS(Tombo_calls, bgTruth, prefix, narrowedCoordinates=narrowCoord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)
    #     df = df[["prefix", "coord", "accuracy", "roc_auc", "precision_5C", "recall_5C", "F1_5C", "Csites", "precision_5mC", "recall_5mC", "F1_5mC", "mCsites", "referenceCpGs", "corrMix", "Corr_mixedSupport", "corrAll", "Corr_allSupport"]]
    #     df.to_csv("{}.report.tsv".format(prefix), sep='\t')
    #     print(df.head())
    #
    # # tool="Nanopolish"
    # if tool == "Nanopolish":
    #     prefix = "{}.{}_calls.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground".format(RunPrefix, tool)
    #     secondFilterBed = "{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed".format(RunPrefix)  # filter coordinates for stats like AUC, F1 etc.
    #     secondFilterBed_4Corr = "{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.4Corr.bed".format(RunPrefix)  # filter coordinates for correlation
    #     df = combine_ONT_and_BS(Nanopolish_calls, bgTruth, prefix, narrowedCoordinates=narrowCoord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)
    #     df = df[["prefix", "coord", "accuracy", "roc_auc", "precision_5C", "recall_5C", "F1_5C", "Csites", "precision_5mC", "recall_5mC", "F1_5mC", "mCsites", "referenceCpGs", "corrMix", "Corr_mixedSupport", "corrAll", "Corr_allSupport"]]
    #     df.to_csv("{}.report.tsv".format(prefix), sep='\t')
    #     print(df.head())
    #
    # # tool="DeepMod"
    # if tool == "DeepMod":
    #     prefix = "{}.{}_calls.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground".format(RunPrefix, tool)
    #     secondFilterBed = "{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed".format(RunPrefix)  # filter coordinates for stats like AUC, F1 etc.
    #     secondFilterBed_4Corr = "{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.4Corr.bed".format(RunPrefix)  # filter coordinates for correlation
    #     df = combine_ONT_and_BS(DeepMod_calls, bgTruth, prefix, narrowedCoordinates=narrowCoord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)
    #     df = df[["prefix", "coord", "accuracy", "roc_auc", "precision_5C", "recall_5C", "F1_5C", "Csites", "precision_5mC", "recall_5mC", "F1_5mC", "mCsites", "referenceCpGs", "corrMix", "Corr_mixedSupport", "corrAll", "Corr_allSupport"]]
    #     df.to_csv("{}.report.tsv".format(prefix), sep='\t')
    #     print(df.head())
