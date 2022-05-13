"""
Sanity check for 5hmC detected by nanopore tools but not by BS-seq,
generate the agreement and discrepancy area between tool and BS-seq.

We will use these two areas to check the 5hmC level: we found that the discrepancy have a much higher 5hmC level.
"""
import gzip

import pandas as pd
from scipy.stats import pearsonr, wilcoxon

from nanome.common.eval_common import import_bgtruth, import_call, readLevelToSiteLevelWithCov, combineBGTruthList
from nanome.common.global_config import *


def output_cpg_set_0base(cpg_set, outfn, bgtruthDict=None):
    if bgtruthDict is None:
        raise Exception("Please specify bgtruthDict")
    outf = gzip.open(outfn, 'wt')
    for key in cpg_set:
        # key is (chr, start, strand)
        if isinstance(bgtruthDict[key], list) or isinstance(bgtruthDict[key], tuple):
            cov = bgtruthDict[key][1]
        else:
            cov = None
        if cov is None:
            outf.write(f"{key[0]}\t{key[1] - 1}\t{key[1]}\t.\t.\t{key[2]}\n")
        else:
            outf.write(f"{key[0]}\t{key[1] - 1}\t{key[1]}\t.\t.\t{key[2]}\t{cov}\n")
    outf.close()
    logger.info(f"save to {outfn}")


def correlation_calls(call1, call2, keySet):
    v1 = []
    v2 = []
    for key in keySet:
        if isinstance(call1[key], list) or isinstance(call1[key], tuple):
            vv1 = call1[key][0]
        else:
            vv1 = call1[key]

        if isinstance(call2[key], list) or isinstance(call2[key], tuple):
            vv2 = call2[key][0]
        else:
            vv2 = call2[key]
        v1.append(vv1)
        v2.append(vv2)
    coe, pvalue = pearsonr(v1, v2)
    return coe, pvalue


def wilcoxon_test_calls(call1, call2, keySet):
    v1 = []
    v2 = []
    for key in keySet:
        if isinstance(call1[key], list) or isinstance(call1[key], tuple):
            vv1 = call1[key][0]
        else:
            vv1 = call1[key]

        if isinstance(call2[key], list) or isinstance(call2[key], tuple):
            vv2 = call2[key][0]
        else:
            vv2 = call2[key]
        v1.append(vv1)
        v2.append(vv2)
    stats, pvalue = wilcoxon(v1, v2)
    return stats, pvalue
    pass


if __name__ == '__main__':
    set_log_debug_level()
    dsname = "APL"
    # dsname = "NA12878"

    logger.info(f"Processing {dsname}:")
    if dsname == "APL":
        mlmlFilename = "/pod/2/li-lab/Nanopore_compare/data/APL_5hmC_BSseq/APL.mlml.addstrand.bscov.oxbscov.selected.bed.gz"
        mlmlDict = import_bgtruth(mlmlFilename, encode="5hmc_ziwei", covCutoff=1, includeCov=True)
        logger.info(f"mlmlDict={len(mlmlDict):,}")

        bgtruthFilename = "/projects/li-lab/Nanopore_compare/data/APL/APL-bs_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.convert.add.strand.tsv.gz"
        # bgtruthDict is key -> (meth_freq, cov)
        bgtruthDict = import_bgtruth(bgtruthFilename, encode="bismark", covCutoff=1, includeCov=True)
        logger.info(f"bgtruthDict={len(bgtruthDict):,}")

        newBgtruthDict = {}
        joinKeys = set(bgtruthDict.keys()).intersection(set(mlmlDict.keys()))
        for key in joinKeys:
            newBgtruthDict[key] = bgtruthDict[key]
        logger.info(f"newBgtruthDict={len(newBgtruthDict):,}")

        # megalodonDict: key-> (meth_freq, cov)
        megalodonFilename = "/pod/2/li-lab/Nanopore_compare/data/APL/APL.megalodon.per_read.sorted.bed.gz"
        megalodonDict = import_call(megalodonFilename, encode="Megalodon.ZW")
        megalodonDict = readLevelToSiteLevelWithCov(megalodonDict, minCov=3, toolname="Megalodon")

        # nanopolishDict: key-> (meth_freq, cov)
        nanopolishFilename = "/pod/2/li-lab/Nanopore_compare/data/APL/APL.nanopolish.methylation_calls.combine.tsv.gz"
        nanopolishDict = import_call(nanopolishFilename, encode="Nanopolish")
        nanopolishDict = readLevelToSiteLevelWithCov(nanopolishDict, minCov=3, toolname="Nanopolish")

        # DeepSignal
        deepsignalFilename = "/pod/2/li-lab/Nanopore_compare/data/APL/APL.deepsignal.call_mods.combine.tsv.gz"
        deepsignalDict = import_call(deepsignalFilename, encode="DeepSignal")
        deepsignalDict = readLevelToSiteLevelWithCov(deepsignalDict, minCov=3, toolname="DeepSignal")

        # Guppy
        guppyFilename = "/pod/2/li-lab/Nanopore_compare/data/APL/APL.Guppy.fast5mod.persite.raw.tsv.gz"
        guppyDict = import_call(guppyFilename, encode="Guppy")
        guppyDict = readLevelToSiteLevelWithCov(guppyDict, minCov=3, toolname="Guppy")

        logger.info(
            f"newBgtruthDict={len(newBgtruthDict)}, megalodonDict={len(megalodonDict)}, nanopolishDict={len(nanopolishDict)}, deepsignalDict={len(deepsignalDict)}, guppyDict={len(guppyDict)}")
    elif dsname == "NA12878":
        enable_cache = True
        using_cache = True

        bgtruth1Filename = "/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz"
        # bgtruthDict is key -> (meth_freq, cov)
        bgtruth1Dict = import_bgtruth(bgtruth1Filename, encode="encode", covCutoff=1, includeCov=True,
                                      enable_cache=enable_cache, using_cache=using_cache)
        logger.info(f"bgtruth1Dict={len(bgtruth1Dict):,}")

        bgtruth2Filename = "/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz"
        # bgtruthDict is key -> (meth_freq, cov)
        bgtruth2Dict = import_bgtruth(bgtruth2Filename, encode="encode", covCutoff=1, includeCov=True,
                                      enable_cache=enable_cache, using_cache=using_cache)
        logger.info(f"bgtruth2Dict={len(bgtruth2Dict):,}")

        bgtruthDict = combineBGTruthList([bgtruth1Dict, bgtruth2Dict], covCutoff=1)
        logger.info(f"combined bgtruthDict={len(bgtruthDict):,}")

        newBgtruthDict = bgtruthDict

        # newBgtruthDict = {}
        # joinKeys = set(bgtruthDict.keys()).intersection(set(mlmlDict.keys()))
        # for key in joinKeys:
        #     newBgtruthDict[key] = bgtruthDict[key]
        # logger.info(f"newBgtruthDict={len(newBgtruthDict):,}")

        # megalodonDict: key-> (meth_freq, cov)
        megalodonFilename = "/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.megalodon.per_read.combine_allchrs.bed.gz"
        megalodonDict = import_call(megalodonFilename, encode="Megalodon", enable_cache=enable_cache,
                                    using_cache=using_cache)
        megalodonDict = readLevelToSiteLevelWithCov(megalodonDict, minCov=3, toolname="Megalodon")

        # nanopolishDict: key-> (meth_freq, cov)
        nanopolishFilename = "/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.nanopolish.methylation_calls.combine_allchrs.tsv.gz"
        nanopolishDict = import_call(nanopolishFilename, encode="Nanopolish", enable_cache=enable_cache,
                                     using_cache=using_cache)
        nanopolishDict = readLevelToSiteLevelWithCov(nanopolishDict, minCov=3, toolname="Nanopolish")

        # DeepSignal
        deepsignalFilename = "/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.deepsignal.call_mods.combine_allchrs.tsv.gz"
        deepsignalDict = import_call(deepsignalFilename, encode="DeepSignal", enable_cache=enable_cache,
                                     using_cache=using_cache)
        deepsignalDict = readLevelToSiteLevelWithCov(deepsignalDict, minCov=3, toolname="DeepSignal")

        # Guppy
        guppyFilename = "/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.guppy.fast5mod_site_level.combine_allchrs.tsv.gz"
        guppyDict = import_call(guppyFilename, encode="Guppy", enable_cache=enable_cache,
                                using_cache=using_cache)
        guppyDict = readLevelToSiteLevelWithCov(guppyDict, minCov=3, toolname="Guppy")

        logger.info(
            f"megalodonDict={len(megalodonDict)}, nanopolishDict={len(nanopolishDict)}, deepsignalDict={len(deepsignalDict)}, guppyDict={len(guppyDict)}")

    for callDict, callName in zip([megalodonDict, nanopolishDict, deepsignalDict, guppyDict],
                                  ['Megalodon', 'Nanopolish', 'DeepSignal', 'Guppy']):
        logger.debug(f"Analyze tool={callName}, CpGs={len(callDict):,}")
        joinSet = set(newBgtruthDict.keys()).intersection(set(callDict.keys()))

        coe, pvalue = correlation_calls(newBgtruthDict, callDict, joinSet)
        logger.info(f"BS-seq vs. {callName} COE={coe:.3f}, pvalue={pvalue:e}, CpGs={len(joinSet):,}")

        ret_cons_let_05 = set()
        ret_cons_let_10 = set()
        ret_diff_get_20 = set()
        ret_diff_get_40 = set()

        agre_let_05_5mc_bgtruth_list = []
        agre_let_05_5mc_tool_list = []
        agre_let_10_5mc_bgtruth_list = []
        agre_let_10_5mc_tool_list = []

        diff_get_20_5mc_bgtruth_list = []
        diff_get_20_5mc_tool_list = []
        diff_get_40_5mc_bgtruth_list = []
        diff_get_40_5mc_tool_list = []

        for key in joinSet:
            if isinstance(newBgtruthDict[key], list) or isinstance(newBgtruthDict[key], tuple):
                bg_methfreq = newBgtruthDict[key][0]
            else:
                bg_methfreq = newBgtruthDict[key]
            tool_methfreq = callDict[key][0]
            if abs(tool_methfreq - bg_methfreq) < 0.05:
                ret_cons_let_05.add(key)
                agre_let_05_5mc_bgtruth_list.append(bg_methfreq)
                agre_let_05_5mc_tool_list.append(tool_methfreq)

            if abs(tool_methfreq - bg_methfreq) < 0.10:
                ret_cons_let_10.add(key)
                agre_let_10_5mc_bgtruth_list.append(bg_methfreq)
                agre_let_10_5mc_tool_list.append(tool_methfreq)

            if abs(tool_methfreq - bg_methfreq) > 0.2:
                ret_diff_get_20.add(key)
                diff_get_20_5mc_bgtruth_list.append(bg_methfreq)
                diff_get_20_5mc_tool_list.append(tool_methfreq)

            if abs(tool_methfreq - bg_methfreq) > 0.4:
                ret_diff_get_40.add(key)
                diff_get_40_5mc_bgtruth_list.append(bg_methfreq)
                diff_get_40_5mc_tool_list.append(tool_methfreq)

        ## Output for scatter plotting
        df = pd.DataFrame({"BS-seq": agre_let_10_5mc_bgtruth_list, callName: agre_let_10_5mc_tool_list})
        outfn = os.path.join(pic_base_dir, f"{dsname}.outcome2.bsseq.{callName}.let10.agreement.csv.gz")
        df.to_csv(outfn, index=False)

        df = pd.DataFrame({"BS-seq": agre_let_05_5mc_bgtruth_list, callName: agre_let_05_5mc_tool_list})
        outfn = os.path.join(pic_base_dir, f"{dsname}.outcome2.bsseq.{callName}.let05.agreement.csv.gz")
        df.to_csv(outfn, index=False)

        df = pd.DataFrame({"BS-seq": diff_get_20_5mc_bgtruth_list, callName: diff_get_20_5mc_tool_list})
        outfn = os.path.join(pic_base_dir, f"{dsname}.outcome2.bsseq.{callName}.get20.discrepency.csv.gz")
        df.to_csv(outfn, index=False)

        df = pd.DataFrame({"BS-seq": diff_get_40_5mc_bgtruth_list, callName: diff_get_40_5mc_tool_list})
        outfn = os.path.join(pic_base_dir, f"{dsname}.outcome2.bsseq.{callName}.get40.discrepency.csv.gz")
        df.to_csv(outfn, index=False)

        logger.info(
            f"tool={callName}:  ret_cons_let_05={len(ret_cons_let_05)}, ret_cons_let_10={len(ret_cons_let_10)}, ret_diff_get_20={len(ret_diff_get_20)}, ret_diff_get_40={len(ret_diff_get_40)}")

        # stats_let_05, pv_let05 = wilcoxon_test_calls(newBgtruthDict, callDict, ret_cons_let_05)
        # stats_let_10, pv_let10 = wilcoxon_test_calls(newBgtruthDict, callDict, ret_cons_let_10)
        # stats_get_20, pv_get20 = wilcoxon_test_calls(newBgtruthDict, callDict, ret_diff_get_20)
        # stats_get_40, pv_get40 = wilcoxon_test_calls(newBgtruthDict, callDict, ret_diff_get_40)
        #
        # logger.info(
        #     f"tool={callName} wilcoxon test:\n <5% = ({stats_let_05}, {pv_let05});\n <10% = ({stats_let_10},{pv_let10});\n >20% = ({stats_get_20}, {pv_get20});\n >40% = ({stats_get_40}, {pv_get40})\n")

        ## Output for DeepTool plotting 5hmC
        outfn = f"{dsname}.bgtruth.{callName}.regions.cons.let.05.bed.gz"
        output_cpg_set_0base(ret_cons_let_05, os.path.join(pic_base_dir, outfn), bgtruthDict=bgtruthDict)

        outfn = f"{dsname}.bgtruth.{callName}.regions.cons.let.10.bed.gz"
        output_cpg_set_0base(ret_cons_let_10, os.path.join(pic_base_dir, outfn), bgtruthDict=bgtruthDict)

        outfn = f"{dsname}.bgtruth.{callName}.regions.diff.get.20.bed.gz"
        output_cpg_set_0base(ret_diff_get_20, os.path.join(pic_base_dir, outfn), bgtruthDict=bgtruthDict)

        outfn = f"{dsname}.bgtruth.{callName}.regions.diff.get.40.bed.gz"
        output_cpg_set_0base(ret_diff_get_40, os.path.join(pic_base_dir, outfn), bgtruthDict=bgtruthDict)

    logger.info(f"sanity_check 5hmc DONE for {dsname}: generate the agreement and discrepancy regions.")
