import os
import pickle

from sklearn.metrics import roc_curve, auc

from lilab.tcga.global_tcga import logger, pic_base_dir, pkl_base_dir
from lilab.tcga.picture import plot_roc_curve
from nanocompare.meth_stats.Universal_meth_stats_evaluation import importPredictions_DeepMod, importGroundTruth_coverage_output_from_Bismark, nonSingletonsPostprocessing, singletonsPostprocessing, combine2programsCalls, combine2programsCalls_4Corr, \
    importPredictions_Nanopolish_3, importPredictions_DeepSignal3, importPredictions_Tombo3, computePerReadStats_v2_for_roc_auc
from nanocompare.nanocompare_global_settings import nanocompare_basedir


def main():
    """
    read four tools results and bg truth, then join common CpGs, save results for ROC curves
    :return:
    """
    RunPrefix = "NA19240"

    DeepsignalInfn = "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepSignal_runs/NA19240/NA19240_DeepSignal.MethCalls.Joined.tsv"
    DeepmodInfn = "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepMod_runs/NA19240/NA19240.C.combined.bed"
    NanopolishInfn = "/projects/li-lab/NanoporeData/WR_ONT_analyses/NA19240_nanopolish/NA19240.methylation_calls.tsv"
    TomboInfn = "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/NA19240_perChr/bk/NA19240_allChr.bed.CpGs.bed"
    BgtruthInfn = "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/EncodeMethyl/joined_reps/RRBS/extractBismark/NA19240_joined_RRBS.Read_R1.Rep_1_trimmed_bismark_bt2.bismark.cov.gz"

    logger.info(f"Start DeepSignal")
    DeepSignal_calls = importPredictions_DeepSignal3(DeepsignalInfn)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_DeepSignal_runs/K562/K562_DeepSignal.MethCalls.Joined.chr20.tsv")

    outfn = os.path.join(pic_base_dir, f"DeepSignal_Calls_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(DeepSignal_calls, handle)

    # with open(outfn, 'rb') as handle:
    #     b = pickle.load(handle)
    # logger.info(f"open b={b}")
    #
    # return

    logger.info(f"Start Tombo")
    Tombo_calls = importPredictions_Tombo3(TomboInfn)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/K562/K562_Tombo.perReadsStats.chr20.bed")

    outfn = os.path.join(pic_base_dir, f"Tombo_calls_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(Tombo_calls, handle)

    logger.info(f"Start Nanopolish3")
    Nanopolish_calls = importPredictions_Nanopolish_3(NanopolishInfn)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.chr_chr20.tsv", IncludeNonSingletons = True)

    outfn = os.path.join(pic_base_dir, f"Nanopolish_calls_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(Nanopolish_calls, handle)

    logger.info(f"Start DeepMod")
    DeepMod_calls = importPredictions_DeepMod(DeepmodInfn)  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/Leukemia_ONT/K562.nanopolish/K562.methylation_calls.chr_chr20.tsv", IncludeNonSingletons = True)

    outfn = os.path.join(pic_base_dir, f"DeepMod_calls_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(DeepMod_calls, handle)

    minCov = 5
    logger.info(f"Start bgTruth")
    bgTruth = importGroundTruth_coverage_output_from_Bismark(BgtruthInfn, covCutt=minCov)

    outfn = os.path.join(pic_base_dir, f"bgTruth_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(bgTruth, handle)

    ####################################### Preprocessing neccessary files:

    singletonsFile = "hg38_singletons.bed"
    nonsingletonsFile = "hg38_nonsingletons.bed"

    narrowCoord = [False, singletonsFile, nonsingletonsFile, "ONT.hg38.cpgIslandExt.bed", "ONT.hg38.cpgShoresExt.bed", "ONT.hg38.cpgShelvesExt.bed", "ONT.hg38.exonFeature.bed", "ONT.hg38.geneFeature.bed", "ONT.hg38.intergenic.bed", "ONT.hg38.intronFeature.bed", "ONT.hg38.promoterFeature.flank_100.bed", "ONT.hg38.promoterFeature.flank_1000.bed",
            "ONT.hg38.promoterFeature.flank_200.bed", "ONT.hg38.promoterFeature.flank_2000.bed", "ONT.hg38.promoterFeature.flank_500.bed", "ONT.hg38.promoterFeature.flank_750.bed"]

    ## add missing region files:
    singletonsFilePrefix = singletonsFile.replace(".bed", '')
    narrowCoord.append("{}.{}.mixed.bed".format(RunPrefix, singletonsFilePrefix))
    narrowCoord.append("{}.{}.absolute.bed".format(RunPrefix, singletonsFilePrefix))
    nonsingletonsFilePrefix = nonsingletonsFile.replace(".bed", '')
    narrowCoord.append("{}.{}.other.bed".format(RunPrefix, nonsingletonsFilePrefix))
    narrowCoord.append("{}.{}.fullyMixed.bed".format(RunPrefix, nonsingletonsFilePrefix))
    narrowCoord.append("{}.{}.discordant.bed".format(RunPrefix, nonsingletonsFilePrefix))
    narrowCoord.append("{}.{}.concordant.bed".format(RunPrefix, nonsingletonsFilePrefix))

    nonSingletonsPostprocessing(bgTruth, nonsingletonsFile, RunPrefix)
    singletonsPostprocessing(bgTruth, singletonsFile, RunPrefix)

    ####################################### Computing stats per program:

    logger.info("\n\n############\n\n")

    DeepSignal_Tombo = combine2programsCalls(DeepSignal_calls, Tombo_calls)
    DeepSignal_Tombo_Nanopolish = combine2programsCalls(DeepSignal_Tombo, Nanopolish_calls)
    DeepSignal_Tombo_Nanopolish_DeepMod = combine2programsCalls(DeepSignal_Tombo_Nanopolish, DeepMod_calls)
    DeepSignal_Tombo_Nanopolish_DeepMod_Bacground = combine2programsCalls(DeepSignal_Tombo_Nanopolish_DeepMod, bgTruth, outfileName="{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed".format(RunPrefix))

    logger.info(f"Data points for all stats: {len(DeepSignal_Tombo_Nanopolish_DeepMod_Bacground)}")

    DeepSignal_Tombo_corr = combine2programsCalls_4Corr(DeepSignal_calls, Tombo_calls)
    DeepSignal_Tombo_Nanopolish_corr = combine2programsCalls_4Corr(DeepSignal_Tombo_corr, Nanopolish_calls)
    DeepSignal_Tombo_Nanopolish_DeepMod_corr = combine2programsCalls_4Corr(DeepSignal_Tombo_Nanopolish_corr, DeepMod_calls)
    # DeepSignal_Tombo_Nanopolish_Bacground = combine2programsCalls_4Corr(DeepSignal_Tombo_Nanopolish_corr, oxBS, outfileName = "DeepSignal_Tombo_Nanopolish_Bacground.4Corr.bed")
    DeepSignal_Tombo_Nanopolish_DeepMod_Bacground_corr = combine2programsCalls(DeepSignal_Tombo_Nanopolish_DeepMod_corr, bgTruth, outfileName="{}.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.4Corr.bed".format(RunPrefix))

    logger.info(f"Data points for correlation: {len(DeepSignal_Tombo_Nanopolish_DeepMod_Bacground_corr)}")

    logger.info("\n\n############\n\n")

    DeepSignalJoined = {key: DeepSignal_calls[key] for key in DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.keys()}

    DeepModJoined = {key: DeepMod_calls[key] for key in DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.keys()}

    NanopolishJoined = {key: Nanopolish_calls[key] for key in DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.keys()}

    TomboJoined = {key: Tombo_calls[key] for key in DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.keys()}

    BgtruthJoined = {key: bgTruth[key] for key in DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.keys()}

    outfn = os.path.join(pic_base_dir, f"DeepSignalJoined_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(DeepSignalJoined, handle)

    outfn = os.path.join(pic_base_dir, f"DeepModJoined_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(DeepModJoined, handle)

    outfn = os.path.join(pic_base_dir, f"NanopolishJoined_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(NanopolishJoined, handle)

    outfn = os.path.join(pic_base_dir, f"TomboJoined_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(TomboJoined, handle)

    outfn = os.path.join(pic_base_dir, f"BgtruthJoined_{RunPrefix}.pkl")
    with open(outfn, 'wb') as handle:
        pickle.dump(BgtruthJoined, handle)

    pass


def load_joined(tool="Nanopolish", dsname="NA19240"):
    infn = os.path.join(pkl_base_dir, "nanocompare", "NA19240_AUC_cut5", f"{tool}Joined_{dsname}.pkl")
    with open(infn, 'rb') as handle:
        cpgDict = pickle.load(handle)
    # logger.info(f"open {tool} for {dsname} cpgDict={cpgDict}")
    return cpgDict
    pass


def load_joined3(tool="Nanopolish", dsname="NA19240"):
    infn = os.path.join(pkl_base_dir, 'nanocompare', 'NA19240_AUC_cut5', f"{tool}Joined_{dsname}.pkl")
    with open(infn, 'rb') as handle:
        cpgDict = pickle.load(handle)
    # logger.info(f"open {tool} for {dsname} cpgDict={cpgDict}")
    return cpgDict
    pass


if __name__ == '__main__':

    # main()
    # sys.exit(0)

    na19240_flist = ['NA19240.hg38_nonsingletons.concordant.bed', 'NA19240.hg38_nonsingletons.discordant.bed']
    na19240_flist = [os.path.join("/projects/liuya/results/pkl/nanocompare/NA19240_AUC_cut5/", co) for co in na19240_flist]

    singletonsFile = "hg38_singletons.bed"
    nonsingletonsFile = "hg38_nonsingletons.bed"

    narrowCoord = [singletonsFile, nonsingletonsFile, "ONT.hg38.cpgIslandExt.bed", "ONT.hg38.exonFeature.bed", "ONT.hg38.intergenic.bed", "ONT.hg38.promoterFeature.flank_500.bed", "ONT.hg38.intronFeature.bed"]  # "ONT.hg38.intronFeature.bed",

    narrowCoord = [os.path.join(nanocompare_basedir, "reports", co) for co in narrowCoord]

    narrowCoord = [False] + na19240_flist + narrowCoord

    narrowCoord_abbr = ['Genome-wide', 'Concordant', 'Discordant', 'Singletons', 'Non-singletons', 'CpG Island', 'Exon', 'Intergenic', 'Promoters', 'Intron']

    cpgDict1 = load_joined3()
    cpgDict2 = load_joined3(tool="DeepMod")
    cpgDict3 = load_joined3(tool="DeepSignal")
    cpgDict4 = load_joined3(tool="Tombo")
    #
    cpgDict5 = load_joined(tool="Bgtruth")

    cpgDictList = [cpgDict3, cpgDict4, cpgDict2, cpgDict1]

    toolsName = ['DeepSignal', 'Tombo', 'DeepMod', 'Nanopolish']

    for t, cordBed in enumerate(narrowCoord):
        cordAbbr = narrowCoord_abbr[t]
        logger.debug(f"cordBed={cordBed}, cordAbbr={cordAbbr}")

        fpr = dict()
        tpr = dict()
        roc_auc = dict()

        i = 0
        for cpgDict in cpgDictList:

            y_score, y = computePerReadStats_v2_for_roc_auc(cpgDict, cpgDict5, title="NA19240", bedFile=cordBed)
            fpr[i], tpr[i], _ = roc_curve(y, y_score)
            roc_auc[i] = auc(fpr[i], tpr[i])
            i = i + 1

        plot_roc_curve(fpr, tpr, roc_auc=roc_auc, labels=toolsName, coord=cordAbbr)
