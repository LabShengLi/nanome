import csv
import os
import pickle
import sys

from lilab.tcga.global_tcga import set_log_debug_level, logger, results_base_dir, pkl_base_dir, ensure_dir, pic_base_dir
from nanocompare.meth_stats.Universal_meth_stats_evaluation import importPredictions_DeepMod, importPredictions_DeepSignal, importPredictions_Tombo, importPredictions_Nanopolish_v2, importGroundTruth_BedMethyl_from_Encode, importGroundTruth_oxBS, importGroundTruth_coverage_output_from_Bismark, nonSingletonsPostprocessing2, singletonsPostprocessing2, \
    combine2programsCalls, combine2programsCalls_4Corr, report_per_read_performance, save_ontcalls_to_pkl, load_ontcalls_pkl, importGroundTruth_coverage_output_from_Bismark_BedGraph
from nanocompare.global_settings import tools_abbr, narrowCoord, singletonsFile, nonsingletonsFile
import study.venn as venn
import matplotlib.pyplot as plt


def save_only_tool_ontcalls(tsvFilename):
    """
    Try save ONT calls of each method into pkl, save to pkl/nanocompare/ONTCalls/dsname/

    Note: this will only scan all tool's ONT calls, not deal with bg truth

    :param tsvFilename:
    :return:
    """
    logger.debug(f"load config file:{tsvFilename}")

    infile = open(tsvFilename, 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')

    for row in csvfile:
        if row['status'] == "submit":
            logger.debug(f"read row={row}")

            dsname = row['Dataset']
            runPrefix = row['RunPrefix']
            encoder = row['parser']
            bgtruth_fn = row['bgTruth']
            minCov = int(row['minCov'])

            outdirONT = os.path.join(pkl_base_dir, 'nanocompare', 'ONTCalls', dsname)
            ensure_dir(outdirONT)

            DeepSignal_calls = importPredictions_DeepSignal(row['DeepSignal_calls'])
            outfn = os.path.join(outdirONT, f"Ontcalls.DeepSignal.{dsname}.pkl")
            save_ontcalls_to_pkl(DeepSignal_calls, outfn)

            Tombo_calls = importPredictions_Tombo(row['Tombo_calls'])
            outfn = os.path.join(outdirONT, f"Ontcalls.Tombo.{dsname}.pkl")
            save_ontcalls_to_pkl(Tombo_calls, outfn)

            Nanopolish_calls = importPredictions_Nanopolish_v2(row['Nanopolish_calls'])
            outfn = os.path.join(outdirONT, f"Ontcalls.Nanopolish.{dsname}.pkl")
            save_ontcalls_to_pkl(Nanopolish_calls, outfn)

            DeepMod_calls = importPredictions_DeepMod(row['DeepMod_calls'])
            outfn = os.path.join(outdirONT, f"Ontcalls.DeepMod.{dsname}.pkl")
            save_ontcalls_to_pkl(DeepMod_calls, outfn)
    logger.info("####### save_only_tool_ontcalls finished!")


def save_bgtruth_and_singleton_bed(tsvFilename):
    """
    read config file, load ontcalls

    save bgtruth, joined bed files as pkl data and .bed file to pkl/nanocompare/runprefix/

    :param tsvFilename:
    :return:
    """
    logger.debug(f"load config file:{tsvFilename}")

    infile = open(tsvFilename, 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')

    for row in csvfile:
        if row['status'] == "submit":
            logger.debug(f"read row={row}")

            dsname = row['Dataset']
            runPrefix = row['RunPrefix']
            encoder = row['parser']
            bgtruth_fn = row['bgTruth']
            minCov = int(row['minCov'])

            outdir = os.path.join(pkl_base_dir, 'nanocompare', runPrefix)
            ensure_dir(outdir)

            DeepSignal_calls = load_ontcall_by_tool(dsname, 'DeepSignal')
            # DeepSignal_calls = importPredictions_DeepSignal(row['DeepSignal_calls'])
            # outfn = os.path.join(outdirONT, "DeepSignal_ontcalls.pkl")
            # save_ontcalls_to_pkl(DeepSignal_calls, outfn)

            Tombo_calls = load_ontcall_by_tool(dsname, 'Tombo')
            # Tombo_calls = importPredictions_Tombo(row['Tombo_calls'])
            # outfn = os.path.join(outdirONT, "Tombo_ontcalls.pkl")
            # save_ontcalls_to_pkl(Tombo_calls, outfn)

            Nanopolish_calls = load_ontcall_by_tool(dsname, 'Nanopolish')
            # Nanopolish_calls = importPredictions_Nanopolish_2(row['Nanopolish_calls'])
            # outfn = os.path.join(outdirONT, "Nanopolish_ontcalls.pkl")
            # save_ontcalls_to_pkl(Nanopolish_calls, outfn)

            DeepMod_calls = load_ontcall_by_tool(dsname, 'DeepMod')
            # DeepMod_calls = importPredictions_DeepMod(row['DeepMod_calls'])
            # outfn = os.path.join(outdirONT, "DeepMod_ontcalls.pkl")
            # save_ontcalls_to_pkl(DeepMod_calls, outfn)

            if encoder == "encode":
                bgTruth = importGroundTruth_BedMethyl_from_Encode(bgtruth_fn, covCutt=minCov)
            elif encoder == "oxBS_sudo":
                bgTruth = importGroundTruth_oxBS(bgtruth_fn, covCutt=minCov)
            elif encoder == "bismark":
                bgTruth = importGroundTruth_coverage_output_from_Bismark(bgtruth_fn, covCutt=minCov)
            elif encoder == "bedgraph":
                bgTruth = importGroundTruth_coverage_output_from_Bismark_BedGraph(bgtruth_fn, covCutt=minCov)
                raise Exception(f"Not used encoder={encoder} before, please check code first.")
            else:
                logger.error("UniversalMethStatsEvaluation.standalone_01.py ERROR: Unknown bacground truth parser configuration. Aborting. FYI: currently supported are: encode, oxBS_sudo, bismark")
                raise Exception(f"unknown encoder={encoder}")

            outfn = os.path.join(outdir, f"Ontcalls.bgTruth.{runPrefix}.{encoder}.pkl")
            save_ontcalls_to_pkl(bgTruth, outfn)

            # singletonsFile = "hg38_singletons.bed"
            # nonsingletonsFile = "hg38_nonsingletons.bed"

            # generate fully-methylated coordinate include absolute, concordant, discordant, etc.
            singletonsPostprocessing2(bgTruth, singletonsFile, runPrefix, outdir)
            nonSingletonsPostprocessing2(bgTruth, nonsingletonsFile, runPrefix, outdir)

            logger.info("\n\n############\n\n")

            DeepSignal_Tombo = combine2programsCalls(DeepSignal_calls, Tombo_calls)
            DeepSignal_Tombo_Nanopolish = combine2programsCalls(DeepSignal_Tombo, Nanopolish_calls)
            DeepSignal_Tombo_Nanopolish_DeepMod = combine2programsCalls(DeepSignal_Tombo_Nanopolish, DeepMod_calls)

            outfn = os.path.join(outdir, f"{runPrefix}.joined.bed")
            DeepSignal_Tombo_Nanopolish_DeepMod_Bacground = combine2programsCalls(DeepSignal_Tombo_Nanopolish_DeepMod, bgTruth, outfileName=outfn)

            logger.info(f"Data points for all stats: {len(DeepSignal_Tombo_Nanopolish_DeepMod_Bacground)}")

            # outfn = os.path.join(outdir, f"Joined_ontcalls.pkl")
            # save_ontcalls_to_pkl(DeepSignal_Tombo_Nanopolish_DeepMod_Bacground, outfn)

            DeepSignal_Tombo_corr = combine2programsCalls_4Corr(DeepSignal_calls, Tombo_calls)
            DeepSignal_Tombo_Nanopolish_corr = combine2programsCalls_4Corr(DeepSignal_Tombo_corr, Nanopolish_calls)
            DeepSignal_Tombo_Nanopolish_DeepMod_corr = combine2programsCalls_4Corr(DeepSignal_Tombo_Nanopolish_corr, DeepMod_calls)
            # DeepSignal_Tombo_Nanopolish_Bacground = combine2programsCalls_4Corr(DeepSignal_Tombo_Nanopolish_corr, oxBS, outfileName = "DeepSignal_Tombo_Nanopolish_Bacground.4Corr.bed")
            outfn = os.path.join(outdir, f"{runPrefix}.4Corr.joined.bed")

            DeepSignal_Tombo_Nanopolish_DeepMod_Bacground_corr = combine2programsCalls(DeepSignal_Tombo_Nanopolish_DeepMod_corr, bgTruth, outfileName=outfn)

            logger.info(f"Data points for correlation: {len(DeepSignal_Tombo_Nanopolish_DeepMod_Bacground_corr)}")

            # outfn = os.path.join(outdir, f"Joined4Cor_ontcalls.pkl")
            # save_ontcalls_to_pkl(DeepSignal_Tombo_Nanopolish_DeepMod_Bacground_corr, outfn)

            logger.info("\n\n############\n\n")

    pass
    logger.info("####### save_bgtruth_and_singleton_bed finished!")


def ven4tools(ontCalls, outSuffix=""):
    """
    venn diagram of 4 tools ontCall
    :param ontCalls:
    :return:
    """

    j01 = combine2programsCalls(ontCalls[tools_abbr[0]], ontCalls[tools_abbr[1]])
    j02 = combine2programsCalls(ontCalls[tools_abbr[0]], ontCalls[tools_abbr[2]])
    j03 = combine2programsCalls(ontCalls[tools_abbr[0]], ontCalls[tools_abbr[3]])
    j12 = combine2programsCalls(ontCalls[tools_abbr[1]], ontCalls[tools_abbr[2]])
    j13 = combine2programsCalls(ontCalls[tools_abbr[1]], ontCalls[tools_abbr[3]])
    j23 = combine2programsCalls(ontCalls[tools_abbr[2]], ontCalls[tools_abbr[3]])

    j012 = combine2programsCalls(j01, ontCalls[tools_abbr[2]])
    j013 = combine2programsCalls(j01, ontCalls[tools_abbr[3]])
    j123 = combine2programsCalls(j12, ontCalls[tools_abbr[3]])
    j0123 = combine2programsCalls(j012, ontCalls[tools_abbr[3]])

    labels = {'0001': len(ontCalls[tools_abbr[3]]), '0010': len(ontCalls[tools_abbr[2]]), '0100': len(ontCalls[tools_abbr[1]]), '1000': len(ontCalls[tools_abbr[0]]), '1111': len(j0123)}

    fig, ax = venn.venn4(labels, names=tools_abbr)

    import lilab.tcga.utils as utils
    outfn = os.path.join(pic_base_dir, f'venn_4tools_{outSuffix}_{utils.current_time_str()}.png')
    fig.savefig(outfn, dpi=600, bbox_inches='tight')
    plt.close()
    logger.debug(f"save to {outfn}")

    pass


def ven_diagram_prepare(tsvFilename):
    """
    prepare data for venn diagram

    Note: running save_ontcalls_bgtruth_to_pkls before running report_performance_results

    :param tsvFilename:
    :return:
    """
    logger.debug(f"ven_diagram_prepare, load config file:{tsvFilename}")

    infile = open(tsvFilename, 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')

    for row in csvfile:
        if row['status'] == "submit":
            logger.debug(f"read row={row}\n")

            runPrefix = row['RunPrefix']
            dsname = row['Dataset']
            encoder = row['parser']

            outdir = os.path.join(pkl_base_dir, 'nanocompare', runPrefix, 'venn_diag')
            ensure_dir(outdir)

            bgCall = load_bgtruth(runPrefix, encoder)

            ontCalls = dict()
            for tool in tools_abbr:
                ontCalls[tool] = load_ontcall_by_tool(dsname, tool)
                # # TODO: caculate related to truth results based on
                #
                # secondFilterBed = os.path.join(pkl_base_dir, 'nanocompare', runPrefix, f"{runPrefix}.joined.bed")  # filter coordinates for stats like AUC, F1 etc.
                # secondFilterBed_4Corr = os.path.join(pkl_base_dir, 'nanocompare', runPrefix, f"{runPrefix}.4Corr.joined.bed")  # filter coordinates for correlation
                #
                # # df = combine_ONT_and_BS(DeepSignal_calls, bgTruth, prefix, narrowedCoordinates=narrowCoord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)
                #
                # tmpPrefix = f'{runPrefix}.{tool}'
                # df = report_df_ont_with_bgtruth(toolCalls, bgtruthCalls, tmpPrefix, narrowedCoordinatesList=repCord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)
                #
                # df = df[["prefix", "coord", "accuracy", "roc_auc", "F1_5C", "F1_5mC", "precision_5C", "recall_5C", "precision_5mC", "recall_5mC", "referenceCpGs", "corrMix", "Corr_mixedSupport", "corrAll", "Corr_allSupport", "Csites", "mCsites", "Csites1", "mCsites1"]]
                #
                # df['Tool'] = tool
                # df['Dataset'] = dsname
                #
                # outfn = os.path.join(outdir, f"{runPrefix}.{tool}.report.tsv")
                # df.to_csv(outfn, sep='\t')
                # logger.info(f"save to {outfn}")

            ven4tools(ontCalls, outSuffix=f'{dsname}')
            ret1 = combine2programsCalls(ontCalls['DeepSignal'], ontCalls['DeepMod'])
            ret2 = combine2programsCalls(ontCalls['DeepSignal'], ontCalls['Nanopolish'])

            logger.info(f'DeepSignal and DeepMod join:{len(ret1)}, {len(ret2)}')


def report_performance_results(tsvFilename):
    """
    report performance results based on tsv config

    Note: running save_ontcalls_bgtruth_to_pkls before running report_performance_results

    :param tsvFilename:
    :return:
    """
    logger.debug(f"report_performance_results, load config file:{tsvFilename}")

    infile = open(tsvFilename, 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')

    for row in csvfile:
        if row['status'] == "submit":
            logger.debug(f"read row={row}\n")

            runPrefix = row['RunPrefix']
            dsname = row['Dataset']
            encoder = row['parser']

            repCord = narrowCoord

            ## add recognized region bed files
            repCord1 = []
            singletonsFilePrefix = singletonsFile.replace(".bed", '')
            # repCord1.append("{}.{}.mixed.bed".format(runPrefix, singletonsFilePrefix))
            repCord1.append("{}.{}.absolute.bed".format(runPrefix, singletonsFilePrefix))

            nonsingletonsFilePrefix = nonsingletonsFile.replace(".bed", '')
            # repCord1.append("{}.{}.other.bed".format(runPrefix, nonsingletonsFilePrefix))
            # repCord1.append("{}.{}.fullyMixed.bed".format(runPrefix, nonsingletonsFilePrefix))
            repCord1.append("{}.{}.discordant.bed".format(runPrefix, nonsingletonsFilePrefix))
            repCord1.append("{}.{}.concordant.bed".format(runPrefix, nonsingletonsFilePrefix))
            # repCord1.append(f"{runPrefix}.joined.bed")
            # repCord1.append(f"{runPrefix}.4Corr.joined.bed")

            repCord1 = [os.path.join(pkl_base_dir, 'nanocompare', runPrefix, cofn) for cofn in repCord1]
            repCord = repCord + repCord1
            # logger.debug(f"repCord={repCord}")

            outdir = os.path.join(pkl_base_dir, 'nanocompare', runPrefix, 'performance_results')
            ensure_dir(outdir)

            bgtruthCalls = load_bgtruth(runPrefix, encoder)

            for tool in tools_abbr:
                toolCalls = load_ontcall_by_tool(dsname, tool)
                # TODO: caculate related to truth results based on

                secondFilterBed = os.path.join(pkl_base_dir, 'nanocompare', runPrefix, f"{runPrefix}.joined.bed")  # filter coordinates for stats like AUC, F1 etc.
                secondFilterBed_4Corr = os.path.join(pkl_base_dir, 'nanocompare', runPrefix, f"{runPrefix}.4Corr.joined.bed")  # filter coordinates for correlation

                # df = combine_ONT_and_BS(DeepSignal_calls, bgTruth, prefix, narrowedCoordinates=narrowCoord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)

                tmpPrefix = f'{runPrefix}.{tool}'
                df = report_per_read_performance(toolCalls, bgtruthCalls, tmpPrefix, narrowedCoordinatesList=repCord, secondFilterBed=secondFilterBed, secondFilterBed_4Corr=secondFilterBed_4Corr)

                df = df[["prefix", "coord", "accuracy", "roc_auc", "F1_5C", "F1_5mC", "precision_5C", "recall_5C", "precision_5mC", "recall_5mC", "referenceCpGs", "corrMix", "Corr_mixedSupport", "corrAll", "Corr_allSupport", "Csites", "mCsites", "Csites1", "mCsites1"]]

                df['Tool'] = tool
                df['Dataset'] = dsname

                outfn = os.path.join(outdir, f"{runPrefix}.{tool}.report.tsv")
                df.to_csv(outfn, sep='\t')
                logger.info(f"save to {outfn}")

                logger.info(df.head())


def load_ontcall_by_tool(dsname, tool):
    """
    Load ont Calls from pkl
    :param runPrefix:
    :param tool:
    :return:
    """
    infn = os.path.join(pkl_base_dir, "nanocompare", "ONTCalls", dsname, f"Ontcalls.{tool}.{dsname}.pkl")
    ret = load_ontcalls_pkl(infn)
    logger.debug(f'load ontcall by {tool} from {dsname} ok. CpGs={len(ret)}')
    return ret
    # return ontCalls


def load_bgtruth(runPrefix, encoder):
    """
    Load ont Calls from pkl
    :param runPrefix:
    :param tool:
    :return:
    """
    infn = os.path.join(pkl_base_dir, "nanocompare", runPrefix, f"Ontcalls.bgTruth.{runPrefix}.{encoder}.pkl")
    ret = load_ontcalls_pkl(infn)
    logger.debug(f'load bgtruth encoded by {encoder} from {runPrefix} ok.  CpGs={len(ret)}')

    return ret
    # with open(infn, 'rb') as handle:
    #     ontCalls = pickle.load(handle)
    # return ontCalls


def test1():
    dsname = 'NA19240'
    tool = 'DeepSignal'

    for tool in tools_abbr:
        logger.info(f'dsname={dsname}, tool={tool}')
        ONTCall = load_ontcall_by_tool(dsname=dsname, tool=tool)
        logger.info(f'tool={tool}, ONTCall={len(ONTCall)}')


if __name__ == '__main__':
    set_log_debug_level()
    #
    # test1()
    # sys.exit(0)

    # tsvFilename = '/projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComarePerformance_deprecated.tsv'

    # tsvFilename = '/projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComarePerformance_paper.tsv'

    tsvFilename = '/projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComarePerformance_paper2.tsv'

    # save_only_tool_ontcalls(tsvFilename)

    # save_bgtruth_and_singleton_bed(tsvFilename)
    # report_performance_results(tsvFilename)

    # ven_diagram_prepare(tsvFilename)
